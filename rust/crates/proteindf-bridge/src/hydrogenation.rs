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

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::ccd_templates::CcdBondTemplate;
use crate::error::{BridgeError, Result};
use crate::position::Position;
use crate::superposer::Superposer;

/// Minimum number of common heavy atoms required to superpose 3D structures reliably.
pub const MIN_SUPERPOSE_HEAVY_ATOMS: usize = 3;

/// Options for configuring hydrogen addition.
#[derive(Debug, Clone, Default)]
pub struct HydrogenationOptions<'a> {
    /// Optional explicit list of heavy atom names to use for superposition.
    /// If `None`, all shared heavy atoms (excluding peptide backbone carbonyl oxygen 'O'
    /// in internal residues) are used.
    pub fit_heavy_atoms: Option<&'a [&'a str]>,
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
/// - Fewer than [`MIN_SUPERPOSE_HEAVY_ATOMS`] (3) matching heavy atoms with known coordinates
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

    // When a template has a terminal carboxylate oxygen (OXT) but the component does not,
    // the component's 'O' is a backbone peptide carbonyl oxygen rather than a free carboxylate oxygen.
    // In free amino acids (CCD entries), the carboxylate 'O' conformation differs by ~170° around CA-C from
    // the peptide backbone carbonyl 'O' (dependent on secondary structure psi angle).
    // Including it in rigid-body superposition severely distorts the fit (RUST_PORT_SPEC.md §3.16).
    let is_internal_peptide_residue = template.get_atom("OXT").is_some()
        && component.pickup_atoms("OXT").is_empty()
        && !component.pickup_atoms("N").is_empty()
        && !component.pickup_atoms("CA").is_empty()
        && !component.pickup_atoms("C").is_empty();

    for ccd_atom in &template.atoms {
        if ccd_atom.is_hydrogen() {
            continue;
        }

        if let Some(allowed) = options.fit_heavy_atoms {
            if !allowed.contains(&ccd_atom.name.as_str()) {
                continue;
            }
        } else if is_internal_peptide_residue && ccd_atom.name == "O" {
            // Skip backbone carbonyl O to avoid distortion from psi angle discrepancy
            continue;
        }

        let Some((ix, iy, iz)) = ccd_atom.ideal_xyz else {
            continue;
        };

        // Check if component has this atom
        let actual_atoms = component.pickup_atoms(&ccd_atom.name);
        if let Some(actual_atom) = actual_atoms.first() {
            let mut ref_atom = Atom::new_with_pos(&ccd_atom.element, Position::new(ix, iy, iz))?;
            ref_atom.name = ccd_atom.name.clone();
            template_heavy_group.set_atom(&ccd_atom.name, ref_atom);
            actual_heavy_group.set_atom(&ccd_atom.name, (*actual_atom).clone());
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
    let superposer = Superposer::new(&template_heavy_group, &actual_heavy_group)?;

    // 3. Identify missing hydrogens and transfer them with transformed coordinates
    let mut added_hydrogens = 0;
    let mut added_atom_names = Vec::new();

    for ccd_atom in &template.atoms {
        if !ccd_atom.is_hydrogen() {
            continue;
        }

        // Skip if component already has this atom
        if !component.pickup_atoms(&ccd_atom.name).is_empty() {
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

        // Apply rigid transformation: pos' = R * (pos - center1) + center2
        let mut pos = Position::new(ix, iy, iz);
        pos -= superposer.center1();
        pos.rotate(superposer.rotation_mat())?;
        pos += superposer.center2();

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
