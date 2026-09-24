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

use std::collections::HashSet;

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

    // Identify distorted atoms to exclude if auto-exclusion is enabled.
    // In free CCD templates (monomers), terminal functional groups (e.g. carboxylate C(=O)OXT,
    // nucleotide 5'-phosphate) adopt conformations that differ drastically (~170° psi rotation)
    // from polymer-internal backbone conformations.
    // We detect this generically using template topology: if a template heavy atom Y is missing
    // in the component (e.g. OXT, OP3) and shares a branching center C with atom X (e.g. O),
    // atom X reflects an unpolymerized terminal geometry and must not be used as a rigid anchor.
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

        // Exclude distorted terminal atoms (even when explicitly named, unless auto_exclude is disabled)
        if distorted_atoms.contains(&ccd_atom.name) {
            continue;
        }

        let Some((ix, iy, iz)) = ccd_atom.ideal_xyz else {
            continue;
        };

        // Lookup atom in component using fast O(1) checks (get_atom / fallback to pickup_atoms)
        if let Some(actual_atom) = component
            .get_atom(&ccd_atom.name)
            .cloned()
            .or_else(|| component.pickup_atoms(&ccd_atom.name).first().cloned())
        {
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

        // Fast O(1) check if component already has this atom
        if component.has_atom(&ccd_atom.name) || !component.pickup_atoms(&ccd_atom.name).is_empty()
        {
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

        // Transform position using Superposer::transform_position (no duplicated math)
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
/// nucleotide with 5'-terminal OP3/O3P), the terminal group forms a planar or tetrahedral
/// branched carboxylate / phosphate with specific dihedral angles.
/// When the actual component is embedded inside a polymer chain, the missing terminal
/// capping atom indicates that the branching center is connected to the next residue,
/// making the remaining terminal atom's dihedral angle dependent on polymer conformation
/// (e.g. protein secondary structure psi angle, causing ~170° discrepancy).
fn detect_distorted_terminal_atoms(
    component: &AtomGroup,
    template: &CcdBondTemplate,
) -> HashSet<String> {
    let mut distorted = HashSet::new();

    // Collect all heavy atoms present in the template but missing in the component
    let missing_template_heavy: HashSet<&str> = template
        .atoms
        .iter()
        .filter(|a| {
            !a.is_hydrogen()
                && !component.has_atom(&a.name)
                && component.pickup_atoms(&a.name).is_empty()
        })
        .map(|a| a.name.as_str())
        .collect();

    if missing_template_heavy.is_empty() {
        return distorted;
    }

    // Known terminal capping atom names across amino acids and nucleic acids
    // (OXT in amino acids, OP3/HOP3/H3T in nucleotides)
    let is_known_terminal_capping_atom =
        |name: &str| -> bool { matches!(name, "OXT" | "OP3" | "HOP3" | "H3T" | "O1P" | "O2P") };

    for &missing_atom in &missing_template_heavy {
        if !is_known_terminal_capping_atom(missing_atom) {
            continue;
        }

        // Find branching center atoms bonded to this missing capping atom
        for (a1, a2, _order) in &template.bonds {
            let center_atom = if a1 == missing_atom {
                a2.as_str()
            } else if a2 == missing_atom {
                a1.as_str()
            } else {
                continue;
            };

            // Find sibling atoms bonded to the same center atom
            for (b1, b2, _order2) in &template.bonds {
                let sibling = if b1 == center_atom && b2 != missing_atom {
                    b2.as_str()
                } else if b2 == center_atom && b1 != missing_atom {
                    b1.as_str()
                } else {
                    continue;
                };

                // If the sibling is present in the component and is a terminal/carbonyl oxygen,
                // it is subject to conformational distortion from polymer linkage
                if component.has_atom(sibling) || !component.pickup_atoms(sibling).is_empty() {
                    if let Some(ccd_sibling) = template.get_atom(sibling) {
                        // Carbonyl/terminal oxygen (e.g. 'O' in peptides)
                        if ccd_sibling.element.eq_ignore_ascii_case("O") {
                            distorted.insert(sibling.to_string());
                        }
                    }
                }
            }
        }
    }

    distorted
}
