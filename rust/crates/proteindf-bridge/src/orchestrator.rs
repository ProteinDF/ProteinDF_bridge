// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Top-level hydrogenation orchestrator for macromolecular structures.
//!
//! Provides functions to add missing hydrogens to entire macromolecular assemblies,
//! protein chains, or individual residues by coordinating:
//! 1. **Backbone amide and N-terminal hydrogen construction** ([`crate::backbone_hydrogen`]):
//!    Determines sequence-level polymer connectivity and trans-peptide geometry,
//!    adding backbone amide `H` to internal residues and tetrahedral ammonium
//!    `H1`, `H2`, `H3` (or `H1`, `H2` for proline) to N-termini.
//! 2. **CCD-reference-based sidechain and general hydrogen addition** ([`crate::hydrogenation`]):
//!    Superposes idealized heavy atom templates from the Chemical Component Dictionary
//!    (CCD) to generate sidechain and alpha-carbon hydrogens while excluding main-chain N hydrogens.
//!
//! See `RUST_PORT_SPEC.md` §3.16 for full design details and literature sources.

use std::collections::HashMap;

use crate::atom_group::AtomGroup;
use crate::backbone_hydrogen::add_backbone_hydrogens_to_residue_in_place;
use crate::brd::sort_nicely;
use crate::ccd_templates::CcdTemplateDb;
use crate::error::Result;
use crate::hydrogenation::{
    add_hydrogens_to_component_in_place_with_options, HydrogenationOptions, HydrogenationReport,
};

/// Summary report of the overall hydrogenation process on a structure.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct OverallHydrogenationReport {
    /// Total number of hydrogens added across all residues/components.
    pub total_added_hydrogens: usize,
    /// Total number of hydrogens (or spurious atoms) removed across all residues/components.
    pub total_removed_hydrogens: usize,
    /// Number of residues/components that were modified (at least one hydrogen added or removed).
    pub hydrogenated_residues: usize,
    /// Residues/components where NO changes were made (completely skipped due to missing CCD
    /// template, insufficient heavy atoms such as HOH water, or geometric failure), with reason: `(path, reason)`.
    pub skipped_residues: Vec<(String, String)>,
    /// Errors or warnings encountered during Step 1 (backbone) or Step 2 (sidechain/CCD) processing,
    /// with context: `(path, error_message)`.
    ///
    /// Even if a residue was partially modified (e.g. backbone added but sidechain failed),
    /// the step failure is guaranteed to be recorded here and never silenced.
    pub step_errors: Vec<(String, String)>,
    /// Detailed per-residue hydrogenation reports, keyed by residue path.
    pub residue_reports: HashMap<String, HydrogenationReport>,
}

impl OverallHydrogenationReport {
    /// Creates an empty overall hydrogenation report.
    pub fn new() -> Self {
        Self::default()
    }

    /// Merges a per-residue report into this overall report.
    pub fn record_residue(&mut self, path: String, report: HydrogenationReport) {
        self.total_added_hydrogens += report.added_hydrogens;
        self.total_removed_hydrogens += report.removed_hydrogens;
        if report.added_hydrogens > 0 || report.removed_hydrogens > 0 {
            self.hydrogenated_residues += 1;
        }
        self.residue_reports.insert(path, report);
    }

    /// Records a skipped residue or component where no modifications occurred.
    pub fn record_skipped(&mut self, path: String, reason: String) {
        self.skipped_residues.push((path, reason));
    }

    /// Records a failure or warning during a hydrogenation step.
    pub fn record_error(&mut self, path: String, error: String) {
        self.step_errors.push((path, error));
    }
}

/// Helper: checks if a group can be treated as an amino acid polymer component.
/// Standard amino acids in a polypeptide chain possess backbone 'N' and 'CA' atoms.
fn is_amino_acid_component(group: &AtomGroup) -> bool {
    group.has_atom("N") && group.has_atom("CA")
}

/// Helper: checks if two adjacent residues in sequence order have a plausible peptide bond.
/// Plausible $C_{prev} - N$ distance is between 0.8 and 2.5 Å.
fn is_plausible_peptide_bond(prev_res: &AtomGroup, curr_res: &AtomGroup) -> bool {
    if let (Some(c_prev), Some(n_curr)) = (prev_res.get_atom("C"), curr_res.get_atom("N")) {
        let dist = (c_prev.xyz - n_curr.xyz).length();
        (0.8..=2.5).contains(&dist)
    } else {
        false
    }
}

/// Adds missing hydrogens to a single residue or component.
///
/// If `is_amino_acid` is true, performs:
/// 1. Step 1: Backbone amide / N-terminal hydrogen addition.
/// 2. Step 2: Sidechain / general hydrogen addition via CCD template.
///
/// If not an amino acid (e.g. ligand, water, cofactor), performs only Step 2.
///
/// Failures on individual residues (such as lack of CCD template, insufficient common
/// heavy atoms for rigid superposition like HOH water molecules, or coordinate degeneracy)
/// do not abort structure-wide processing; instead, they are recorded in `report` and
/// traversal continues.
fn hydrogenate_single_residue(
    residue: &mut AtomGroup,
    prev_residue: Option<&AtomGroup>,
    db: &CcdTemplateDb,
    report: &mut OverallHydrogenationReport,
) -> Result<()> {
    let res_path = residue.path().to_string();
    let is_aa = is_amino_acid_component(residue);
    let mut combined_report = HydrogenationReport {
        added_hydrogens: 0,
        added_atom_names: Vec::new(),
        removed_hydrogens: 0,
        removed_atom_names: Vec::new(),
    };

    let mut step1_err = None;
    let mut step2_err = None;

    // Step 1: Backbone hydrogenation for amino acid residues
    if is_aa {
        match add_backbone_hydrogens_to_residue_in_place(residue, prev_residue) {
            Ok(bb_report) => {
                combined_report.added_hydrogens += bb_report.added_hydrogens;
                combined_report
                    .added_atom_names
                    .extend(bb_report.added_atom_names);
                combined_report.removed_hydrogens += bb_report.removed_hydrogens;
                combined_report
                    .removed_atom_names
                    .extend(bb_report.removed_atom_names);
            }
            Err(e) => {
                let err_msg = format!("Backbone hydrogenation error: {e}");
                report.record_error(res_path.clone(), err_msg.clone());
                step1_err = Some(err_msg);
            }
        }
    }

    // Step 2: CCD-reference-based sidechain and general hydrogen addition
    if let Some(template) = db.lookup(&residue.name) {
        let options = HydrogenationOptions::default();
        match add_hydrogens_to_component_in_place_with_options(residue, template, &options) {
            Ok(sc_report) => {
                combined_report.added_hydrogens += sc_report.added_hydrogens;
                combined_report
                    .added_atom_names
                    .extend(sc_report.added_atom_names);
                combined_report.removed_hydrogens += sc_report.removed_hydrogens;
                combined_report
                    .removed_atom_names
                    .extend(sc_report.removed_atom_names);
            }
            Err(e) => {
                // E.g., HOH water has only 1 heavy atom ('O'), failing MIN_SUPERPOSE_HEAVY_ATOMS (3),
                // or degenerate coordinates. Record step failure and continue.
                let err_msg = format!("Sidechain/general hydrogenation error: {e}");
                report.record_error(res_path.clone(), err_msg.clone());
                step2_err = Some(err_msg);
            }
        }
    } else {
        let err_msg = format!("No CCD template found for residue '{}'", residue.name);
        report.record_error(res_path.clone(), err_msg.clone());
        step2_err = Some(err_msg);
    }

    let has_errors = step1_err.is_some() || step2_err.is_some();
    let has_modifications =
        combined_report.added_hydrogens > 0 || combined_report.removed_hydrogens > 0;

    if !has_errors {
        // Both steps succeeded without any errors.
        // Even if 0 hydrogens were added/removed (already fully hydrogenated), record into residue_reports.
        report.record_residue(res_path, combined_report);
    } else if has_modifications {
        // At least one step encountered an error, but the other step succeeded and modified the residue.
        // Record the modifications into residue_reports (the error is already recorded in step_errors).
        report.record_residue(res_path, combined_report);
    } else {
        // No modifications occurred AND errors/missing templates occurred: completely skipped.
        let reason = step2_err
            .or(step1_err)
            .unwrap_or_else(|| "No hydrogens added or removed".to_string());
        report.record_skipped(res_path, reason);
    }

    Ok(())
}

/// Hydrogenates a polymer chain group (where child subgroups are residues/components).
fn hydrogenate_chain(
    chain: &mut AtomGroup,
    db: &CcdTemplateDb,
    report: &mut OverallHydrogenationReport,
) -> Result<()> {
    let mut res_keys = chain.get_group_list();
    sort_nicely(&mut res_keys);

    // Snapshot of previous residue for peptide bond tracking
    let mut prev_residue_snapshot: Option<AtomGroup> = None;

    for key in &res_keys {
        if let Some(residue) = chain.get_group_mut(key) {
            let is_aa = is_amino_acid_component(residue);

            let prev_ref = if is_aa {
                prev_residue_snapshot
                    .as_ref()
                    .filter(|prev| is_plausible_peptide_bond(prev, residue))
            } else {
                None
            };

            hydrogenate_single_residue(residue, prev_ref, db, report)?;

            // Update prev_residue_snapshot for the next residue
            if is_aa && residue.has_atom("C") {
                prev_residue_snapshot = Some(residue.clone());
            } else {
                // Non-amino acid (or missing C atom) breaks polypeptide continuity
                prev_residue_snapshot = None;
            }
        }
    }

    Ok(())
}

/// Recursively traverses an [`AtomGroup`] hierarchy and applies hydrogen addition.
///
/// Automatically identifies:
/// - Chain-like groups: Groups whose children contain atoms (leaf residue/component groups).
/// - Single component groups: Groups containing atoms directly.
/// - Model/Root container groups: Groups containing subgroups that are further traversed.
pub fn hydrogenate_atomgroup(
    group: &mut AtomGroup,
    db: &CcdTemplateDb,
) -> Result<OverallHydrogenationReport> {
    let mut report = OverallHydrogenationReport::new();
    traverse_and_hydrogenate(group, db, &mut report)?;
    Ok(report)
}

fn traverse_and_hydrogenate(
    group: &mut AtomGroup,
    db: &CcdTemplateDb,
    report: &mut OverallHydrogenationReport,
) -> Result<()> {
    if group.get_number_of_groups() == 0 {
        // Leaf group with direct atoms: single residue or component
        if group.get_number_of_atoms() > 0 {
            hydrogenate_single_residue(group, None, db, report)?;
        }
        return Ok(());
    }

    // Check if the children of this group are leaf groups (i.e. this group acts as a chain)
    let all_children_are_leaves = group
        .groups()
        .all(|(_, child)| child.get_number_of_groups() == 0);

    if all_children_are_leaves {
        hydrogenate_chain(group, db, report)?;
    } else {
        // Container group (e.g. Root, Model): sort child keys and recurse
        let mut child_keys = group.get_group_list();
        sort_nicely(&mut child_keys);

        for key in &child_keys {
            if let Some(child) = group.get_group_mut(key) {
                traverse_and_hydrogenate(child, db, report)?;
            }
        }
    }

    Ok(())
}

impl AtomGroup {
    /// Adds missing hydrogens across this entire structure using standard polymer
    /// connectivity and CCD idealized templates.
    ///
    /// Delegates to [`crate::orchestrator::hydrogenate_atomgroup`].
    pub fn add_missing_hydrogens(
        &mut self,
        db: &CcdTemplateDb,
    ) -> Result<OverallHydrogenationReport> {
        hydrogenate_atomgroup(self, db)
    }
}
