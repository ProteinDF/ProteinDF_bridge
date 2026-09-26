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
//! - **Main-chain peptide backbone N hydrogens**: Single-component hydrogenation (PR#38) operates
//!   on individual residues without polymer chain context. Because free-amino-acid CCD templates
//!   model monomeric N-terminal geometries (carrying free amine hydrogens or, in PRO, a secondary amine
//!   N-H) rather than polymer backbone amides, all hydrogens bonded to backbone 'N' (where 'N', 'CA',
//!   and 'C' are present in the template) are excluded in PR#38. Peptide backbone amide N-H (for non-proline
//!   residues) and N-terminal capping/protonation (NH3+) are constructed with proper polymer geometry in PR#39.
//! - **Conformational distortion guard scope**: `detect_distorted_terminal_atoms` is scoped
//!   specifically to protein residues (detecting missing `OXT` and guarding backbone carbonyl `O`).
//!   Nucleic acid phosphate terminal/bridging oxygen distortion detection is not covered in PR#38
//!   (deferred to a separate task with real nucleic acid fixtures, per `RUST_PORT_SPEC.md` §3.16).
//!   Additionally, in experimental X-ray structures where a true C-terminal residue has an unresolved
//!   (unmodeled) `OXT`, the guard may conservatively exclude backbone `O` as if the residue were internal.

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
    /// between free-CCD templates and polymer-bound structures (specifically protein backbone
    /// carbonyl oxygen 'O' when the free carboxylate terminal 'OXT' is absent in internal peptide residues).
    ///
    /// This guard currently applies only to protein C-terminal carboxylate (`OXT`) and backbone
    /// carbonyl `O`. It does not apply to nucleic acids or ligands.
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

    // Build map of hydrogen -> parent heavy atom from template bonds
    let mut h_to_parent: HashMap<&str, &str> = HashMap::new();
    for (a1, a2, _order) in &template.bonds {
        let is_h1 = template.get_atom(a1).is_some_and(|a| a.is_hydrogen());
        let is_h2 = template.get_atom(a2).is_some_and(|a| a.is_hydrogen());
        if is_h1 && !is_h2 {
            h_to_parent.insert(a1.as_str(), a2.as_str());
        } else if is_h2 && !is_h1 {
            h_to_parent.insert(a2.as_str(), a1.as_str());
        }
    }

    // Determine whether the template has a standard amino acid backbone (contains N, CA, and C).
    // In standard amino acids, free monomer CCD templates represent free amines (N with H and H2,
    // or N with H in PRO). Without polymer chain context, single-component hydrogenation cannot
    // know whether a residue is an N-terminus or connected in a chain. Adding monomer N hydrogens
    // causes over-protonation in internal residues (and internal PRO has zero hydrogens on N).
    // Therefore, PR#38 unconditionally skips ALL hydrogens bonded to backbone 'N' in amino acid
    // templates. Polymer backbone amide N-H and N-terminal capping (NH3+) are handled in PR#39.
    let is_amino_acid_template = template.get_atom("N").is_some()
        && template.get_atom("CA").is_some()
        && template.get_atom("C").is_some();

    // 3. Identify missing hydrogens and resolve their positions in a staging buffer.
    // Atomicity guarantee: We compute and validate all new hydrogen atoms into a local buffer first.
    // Only after all missing hydrogens are successfully transformed and constructed do we apply
    // them to `component`. If any hydrogen fails (e.g. missing ideal coordinates), `component`
    // remains completely unmodified.
    let mut staged_hydrogens = Vec::new();

    for ccd_atom in &template.atoms {
        if !ccd_atom.is_hydrogen() {
            continue;
        }

        // Check if component already has this atom
        if component_has_atom(component, &ccd_atom.name) {
            continue;
        }

        // Parent heavy atom check:
        // A hydrogen must only be added if its parent heavy atom is known from template bonds
        // AND actually exists in the component.
        // - If bond information is missing for a hydrogen, skip it safely.
        // - For example, HXT is bonded to OXT; if OXT is absent (as in internal peptide residues),
        //   HXT must NOT be added.
        let Some(&parent_heavy) = h_to_parent.get(ccd_atom.name.as_str()) else {
            continue;
        };
        if !component_has_atom(component, parent_heavy) {
            continue;
        }

        // Bond-topology exclusion rule for peptide backbone N:
        // If the template is an amino acid backbone (has N, CA, C), skip ALL hydrogens bonded
        // to heavy atom "N" regardless of their name (e.g. 'H', 'H2', 'H3' in standard amino acids,
        // or 'H' in PRO). This avoids misfiring on ligands that happen to contain an atom named "N".
        if is_amino_acid_template && parent_heavy == "N" {
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

        staged_hydrogens.push((ccd_atom.name.clone(), h_atom));
    }

    let added_atom_names: Vec<String> = staged_hydrogens
        .iter()
        .map(|(name, _)| name.clone())
        .collect();

    // Apply all validated hydrogens atomically
    let added_hydrogens = staged_hydrogens.len();
    for (name, atom) in staged_hydrogens {
        component.set_atom(&name, atom);
    }

    Ok(HydrogenationReport {
        added_hydrogens,
        added_atom_names,
    })
}

/// Detects heavy atoms in the template whose idealized conformation is distorted
/// relative to the actual component due to polymer connectivity.
///
/// # Topology Rationale (Protein-Specific)
/// When a CCD template represents a free amino acid (with terminal `OXT`), the terminal
/// carboxylate C(=O)OXT adopts a planar geometry with a specific dihedral angle.
/// When the actual component is embedded inside a polypeptide chain (missing `OXT`), the
/// carbonyl oxygen's (`O`) orientation is governed by the backbone dihedral angle (psi),
/// which typically deviates by ~170° from the free CCD template.
///
/// To prevent corrupting the rigid-body superposition, carbonyl `O` is automatically
/// identified and excluded from the superposition fitting set.
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
    // Scoped strictly to protein C-terminal OXT (per RUST_PORT_SPEC.md §3.16 scope decision).
    let is_known_terminal_capping_heavy_atom = |name: &str| -> bool { matches!(name, "OXT") };

    let missing_capping_heavy: Vec<&str> = template
        .atoms
        .iter()
        .filter(|a| !a.is_hydrogen())
        .map(|a| a.name.as_str())
        .filter(|&name| is_known_terminal_capping_heavy_atom(name))
        .filter(|&name| !component_has_atom(component, name))
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
                        // 2. Must be an oxygen (backbone carbonyl 'O' in peptides).
                        // 3. Must be a terminal non-bridging oxygen (heavy_degree == 1, bonded only to center_atom 'C').
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
