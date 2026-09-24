// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! CCD-reference-based hydrogen addition engine.
//!
//! Provides functions to add missing hydrogens to single residues or ligands by
//! superposing shared heavy atoms onto idealized coordinates from the Chemical
//! Component Dictionary (CCD). See `RUST_PORT_SPEC.md` §3.16 for full design details.
//!
//! # Known Limitations
//! - **Rotatable hydrogens**: Hydrogens on rotatable bonds (e.g. sidechain -OH, -SH, -NH3+)
//!   adopt the idealized conformation defined in the CCD template directly without local
//!   energy minimization or hydrogen-bond network optimization.
//! - **Protonation state**: Fixed neutral/standard tautomer states from the CCD are used;
//!   pH-dependent pKa estimation (e.g. PROPKA) is not performed.

use std::collections::{HashMap, HashSet};

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::ccd_templates::CcdBondTemplate;
use crate::error::{BridgeError, Result};
use crate::position::Position;
use crate::superposer::Superposer;

/// Minimum number of common heavy atoms required to superpose 3D structures reliably.
pub const MIN_SUPERPOSE_HEAVY_ATOMS: usize = 3;

/// Options for configuring hydrogen addition.
#[derive(Debug, Clone)]
pub struct HydrogenationOptions<'a> {
    /// Optional explicit list of heavy atom names to use for superposition.
    /// If `Some`, only matching atoms in this list are considered for superposition.
    pub fit_heavy_atoms: Option<&'a [&'a str]>,

    /// Whether to automatically exclude heavy atoms with known conformational discrepancies
    /// between free-CCD templates and polymer-bound structures (e.g. backbone carbonyl oxygen 'O'
    /// when the free carboxylate terminal 'OXT' is absent in internal peptide residues).
    ///
    /// Defaults to `true`. Even when [`fit_heavy_atoms`](Self::fit_heavy_atoms) is explicitly provided,
    /// distorted atoms are guarded against unless this flag is explicitly set to `false`.
    pub auto_exclude_distorted_atoms: bool,
}

impl<'a> Default for HydrogenationOptions<'a> {
    fn default() -> Self {
        Self {
            fit_heavy_atoms: None,
            auto_exclude_distorted_atoms: true,
        }
    }
}

/// Result summary of hydrogen addition to a component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HydrogenationReport {
    /// Number of hydrogens added.
    pub added_hydrogens: usize,
    /// Names of the added hydrogen atoms.
    pub added_atom_names: Vec<String>,
}

/// Helper: finds an atom in a component with fast O(1) direct lookup,
/// falling back to recursive pickup_atoms if needed.
fn find_component_atom(component: &AtomGroup, name: &str) -> Option<Atom> {
    component
        .get_atom(name)
        .cloned()
        .or_else(|| component.pickup_atoms(name).into_iter().next())
}

/// Helper: checks if a component contains an atom matching the given name.
fn component_has_atom(component: &AtomGroup, name: &str) -> bool {
    component.has_atom(name) || !component.pickup_atoms(name).is_empty()
}

/// Adds missing hydrogens to a single component (residue or ligand) using a CCD bond template.
///
/// Superposes the template's idealized heavy-atom coordinates onto the component's existing
/// heavy atoms using the Kabsch algorithm ([`Superposer`]), and transfers missing hydrogens
/// with the computed rigid-body transformation.
///
/// # Arguments
/// * `component` - An [`AtomGroup`] representing a single residue or ligand.
/// * `template` - The [`CcdBondTemplate`] containing idealized geometry.
///
/// # Errors
/// Returns an error if:
/// - Fewer than [`MIN_SUPERPOSE_HEAVY_ATOMS`] (3) matching non-collinear heavy atoms with known coordinates
///   are found between the component and the template.
/// - Any missing hydrogen in the template lacks idealized coordinates (`ideal_xyz` is `None`).
pub fn add_hydrogens_to_component(
    component: &AtomGroup,
    template: &CcdBondTemplate,
) -> Result<AtomGroup> {
    add_hydrogens_to_component_with_options(component, template, &HydrogenationOptions::default())
}

/// Adds missing hydrogens with custom [`HydrogenationOptions`].
pub fn add_hydrogens_to_component_with_options(
    component: &AtomGroup,
    template: &CcdBondTemplate,
    options: &HydrogenationOptions,
) -> Result<AtomGroup> {
    let mut result = component.clone();
    add_hydrogens_to_component_in_place_with_options(&mut result, template, options)?;
    Ok(result)
}

/// In-place variant of [`add_hydrogens_to_component`].
///
/// Modifies `component` directly by adding missing hydrogens and returns a [`HydrogenationReport`].
pub fn add_hydrogens_to_component_in_place(
    component: &mut AtomGroup,
    template: &CcdBondTemplate,
) -> Result<HydrogenationReport> {
    add_hydrogens_to_component_in_place_with_options(
        component,
        template,
        &HydrogenationOptions::default(),
    )
}

/// In-place variant with custom [`HydrogenationOptions`].
pub fn add_hydrogens_to_component_in_place_with_options(
    component: &mut AtomGroup,
    template: &CcdBondTemplate,
    options: &HydrogenationOptions,
) -> Result<HydrogenationReport> {
    // 1. Identify shared heavy atoms with valid coordinates in both component and template
    let mut template_heavy_group = AtomGroup::new();
    let mut actual_heavy_group = AtomGroup::new();

    // Identify distorted terminal atoms to exclude if auto-exclusion is enabled.
    // In free CCD templates (monomers), terminal functional groups (e.g. carboxylate C(=O)OXT,
    // nucleotide 5'-phosphate) adopt conformations that differ drastically (~170° psi rotation)
    // from polymer-internal backbone conformations.
    let distorted_atoms: HashSet<String> = if options.auto_exclude_distorted_atoms {
        detect_distorted_terminal_atoms(component, template)
    } else {
        HashSet::new()
    };

    for ccd_atom in &template.atoms {
        if ccd_atom.is_hydrogen() {
            continue;
        }

        // Apply explicit filter if provided
        if let Some(allowed) = options.fit_heavy_atoms {
            if !allowed.contains(&ccd_atom.name.as_str()) {
                continue;
            }
        }

        // Exclude distorted terminal atoms (guarded even when explicitly named unless auto_exclude is disabled)
        if distorted_atoms.contains(&ccd_atom.name) {
            continue;
        }

        let Some((ix, iy, iz)) = ccd_atom.ideal_xyz else {
            continue;
        };

        // Lookup atom in component using unified helper
        if let Some(actual_atom) = find_component_atom(component, &ccd_atom.name) {
            let mut ref_atom = Atom::new_with_pos(&ccd_atom.element, Position::new(ix, iy, iz))?;
            ref_atom.name = ccd_atom.name.clone();
            template_heavy_group.set_atom(&ccd_atom.name, ref_atom);
            actual_heavy_group.set_atom(&ccd_atom.name, actual_atom);
        }
    }

    let matched_count = template_heavy_group.get_number_of_atoms();
    if matched_count < MIN_SUPERPOSE_HEAVY_ATOMS {
        return Err(BridgeError::value_error(
            "common_heavy_atoms",
            format!(
                "Insufficient common heavy atoms ({matched_count} found, at least {MIN_SUPERPOSE_HEAVY_ATOMS} required) \
                 to superpose component '{}' with CCD template '{}'",
                component.name, template.comp_id
            ),
        ));
    }

    // 2. Compute rigid-body superposition from template frame to actual frame
    // Superposer::new validates non-collinearity / non-degeneracy
    let superposer = Superposer::new(&template_heavy_group, &actual_heavy_group)?;

    // 3. Identify missing hydrogens and transfer them with transformed coordinates
    let mut added_hydrogens = 0;
    let mut added_atom_names = Vec::new();

    for ccd_atom in &template.atoms {
        if !ccd_atom.is_hydrogen() {
            continue;
        }

        // Check if component already has this atom
        if component_has_atom(component, &ccd_atom.name) {
            continue;
        }

        let (ix, iy, iz) = ccd_atom.ideal_xyz.ok_or_else(|| {
            BridgeError::value_error(
                "ideal_xyz",
                format!(
                    "Missing ideal coordinates for hydrogen atom '{}' in CCD template '{}'",
                    ccd_atom.name, template.comp_id
                ),
            )
        })?;

        // Transform position using Superposer::transform_position
        let pos = superposer.transform_position(&Position::new(ix, iy, iz))?;

        let mut h_atom = Atom::new_with_pos(&ccd_atom.element, pos)?;
        h_atom.name = ccd_atom.name.clone();

        component.set_atom(&ccd_atom.name, h_atom);
        added_hydrogens += 1;
        added_atom_names.push(ccd_atom.name.clone());
    }

    Ok(HydrogenationReport {
        added_hydrogens,
        added_atom_names,
    })
}

/// Detects heavy atoms in the template whose idealized conformation is distorted
/// relative to the actual component due to polymer connectivity.
///
/// # Topology Rationale
/// When a CCD template represents a free monomer (e.g. amino acid with terminal OXT,
/// nucleotide with 5'-terminal OP3), the terminal group forms a planar or tetrahedral
/// branched carboxylate / phosphate with specific dihedral angles.
/// When the actual component is embedded inside a polymer chain, the missing terminal
/// capping atom indicates that the branching center is connected to the next residue,
/// making terminal non-bridging atoms' dihedral angles dependent on polymer conformation
/// (e.g. protein secondary structure psi angle, causing ~170° discrepancy).
///
/// Crucially, only non-bridging terminal oxygens (heavy degree == 1, bonded only to the branching center)
/// are excluded. Backbone bridging atoms (such as O5' in nucleotides, which connects both P and C5',
/// heavy degree >= 2) are not terminal and are safely retained.
fn detect_distorted_terminal_atoms(
    component: &AtomGroup,
    template: &CcdBondTemplate,
) -> HashSet<String> {
    let mut distorted = HashSet::new();

    // 1. Build adjacency map for heavy atoms in the template: atom_name -> Vec<neighbor_name>
    let mut heavy_adj: HashMap<&str, Vec<&str>> = HashMap::new();
    for (a1, a2, _order) in &template.bonds {
        let is_h1 = template.get_atom(a1).is_some_and(|a| a.is_hydrogen());
        let is_h2 = template.get_atom(a2).is_some_and(|a| a.is_hydrogen());
        if !is_h1 && !is_h2 {
            heavy_adj.entry(a1.as_str()).or_default().push(a2.as_str());
            heavy_adj.entry(a2.as_str()).or_default().push(a1.as_str());
        }
    }

    // 2. Identify missing capping heavy atoms in the component.
    // Note: We only check heavy capping atoms (e.g. OXT in amino acids, OP3 in nucleotides).
    // Hydrogens (e.g. HOP3, H3T) are not part of heavy-atom superposition and are intentionally excluded here.
    let is_known_terminal_capping_heavy_atom =
        |name: &str| -> bool { matches!(name, "OXT" | "OP3" | "O1P" | "O2P") };

    let missing_capping_heavy: Vec<&str> = template
        .atoms
        .iter()
        .filter(|a| !a.is_hydrogen() && !component_has_atom(component, &a.name))
        .map(|a| a.name.as_str())
        .filter(|&name| is_known_terminal_capping_heavy_atom(name))
        .collect();

    if missing_capping_heavy.is_empty() {
        return distorted;
    }

    // 3. For each missing capping atom, inspect its branching center atom
    for &missing_atom in &missing_capping_heavy {
        if let Some(neighbors) = heavy_adj.get(missing_atom) {
            for &center_atom in neighbors {
                if let Some(siblings) = heavy_adj.get(center_atom) {
                    for &sibling in siblings {
                        if sibling == missing_atom {
                            continue;
                        }

                        // Distorted sibling check:
                        // 1. Must be present in the component.
                        // 2. Must be an oxygen (e.g. carbonyl 'O' in peptides, non-bridging OP1/OP2 in phosphates).
                        // 3. MUST be a terminal non-bridging oxygen (heavy_degree == 1, bonded only to center_atom).
                        //    Bridging oxygens (e.g. O5' in nucleotides, bonded to both P and C5', heavy_degree >= 2)
                        //    are polymer backbone linkages and must NOT be excluded.
                        if component_has_atom(component, sibling) {
                            if let Some(ccd_sibling) = template.get_atom(sibling) {
                                if ccd_sibling.element.eq_ignore_ascii_case("O") {
                                    let heavy_degree =
                                        heavy_adj.get(sibling).map_or(0, |v| v.len());
                                    if heavy_degree == 1 {
                                        distorted.insert(sibling.to_string());
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    distorted
}
