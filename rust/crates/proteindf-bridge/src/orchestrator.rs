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
use crate::error::{BridgeError, Result};
use crate::his_state::{his_states, is_histidine_comp_id, HisState, HisStateOptions};
use crate::hydrogenation::{
    add_hydrogens_to_component_in_place_excluding, HydrogenationOptions, HydrogenationReport,
};

/// Default histidine protonation mode (`"auto"`).
pub const DEFAULT_HIS_MODE: &str = "auto";

/// Histidine ring-nitrogen hydrogenation mode during structure traversal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum HisHydrogenationMode {
    /// Select `HD1`/`HE2` per residue according to precomputed [`his_states`] results
    /// keyed by residue path.
    Auto(HashMap<String, HisState>),
    /// Attach both `HD1` and `HE2` to histidine residues (traditional full protonation behavior).
    Hip,
}

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
    /// Processing failures, geometric errors, or missing sidechain templates on partially modified residues,
    /// with context: `(path, error_message)`.
    ///
    /// For completely unmodified components (such as uncataloged ligands or waters), missing templates
    /// are recorded in [`skipped_residues`](Self::skipped_residues) to avoid noise. However, if a residue
    /// was partially modified (e.g. backbone added but sidechain template missing), the issue is guaranteed
    /// to be recorded here so partial failures are never silenced.
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

/// Outcome of a single hydrogenation step (Step 1 backbone or Step 2 sidechain/general).
#[derive(Debug, Clone, PartialEq, Eq)]
enum StepOutcome {
    /// Step was not applicable (e.g. Step 1 on a non-amino-acid component).
    NotApplicable,
    /// Step completed successfully.
    Success,
    /// Step encountered a genuine processing or geometric failure.
    Failed(String),
    /// CCD template was not found (expected for uncataloged ligands/waters/ions).
    TemplateMissing(String),
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
    his_mode: &HisHydrogenationMode,
    report: &mut OverallHydrogenationReport,
) -> Result<()> {
    let res_path = residue.path().to_string();
    let is_aa = is_amino_acid_component(residue);
    let is_his = is_histidine_comp_id(&residue.name);
    let excluded_his_hydrogens: &[&str] = if is_his && residue.get_number_of_atoms() > 0 {
        match his_mode {
            HisHydrogenationMode::Auto(ref sites_map) => {
                let state = sites_map.get(&res_path).ok_or_else(|| {
                    BridgeError::value_error(
                        &res_path,
                        format!("histidine state not found for residue path: {res_path}"),
                    )
                })?;
                match state {
                    HisState::Hid => &["HE2"],
                    HisState::Hie => &["HD1"],
                    HisState::Hip => &[],
                    HisState::MetalBridge => &["HD1", "HE2"],
                }
            }
            HisHydrogenationMode::Hip => &[],
        }
    } else {
        &[]
    };

    let mut combined_report = HydrogenationReport {
        added_hydrogens: 0,
        added_atom_names: Vec::new(),
        removed_hydrogens: 0,
        removed_atom_names: Vec::new(),
    };

    // Step 1: Backbone hydrogenation for amino acid residues
    let step1_outcome = if is_aa {
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
                StepOutcome::Success
            }
            Err(e) => StepOutcome::Failed(format!("Backbone hydrogenation error: {e}")),
        }
    } else {
        StepOutcome::NotApplicable
    };

    // Step 2: CCD-reference-based sidechain and general hydrogen addition
    let template_opt =
        db.lookup(&residue.name)
            .or_else(|| if is_his { db.lookup("HIS") } else { None });
    let step2_outcome = if let Some(template) = template_opt {
        let options = HydrogenationOptions::default();
        match add_hydrogens_to_component_in_place_excluding(
            residue,
            template,
            &options,
            excluded_his_hydrogens,
        ) {
            Ok(sc_report) => {
                combined_report.added_hydrogens += sc_report.added_hydrogens;
                combined_report
                    .added_atom_names
                    .extend(sc_report.added_atom_names);
                combined_report.removed_hydrogens += sc_report.removed_hydrogens;
                combined_report
                    .removed_atom_names
                    .extend(sc_report.removed_atom_names);
                StepOutcome::Success
            }
            Err(e) => {
                // E.g., HOH water has only 1 heavy atom ('O'), failing MIN_SUPERPOSE_HEAVY_ATOMS (3),
                // or degenerate coordinates.
                StepOutcome::Failed(format!("Sidechain/general hydrogenation error: {e}"))
            }
        }
    } else {
        // Missing CCD template is expected for uncataloged ligands, waters, or metal ions.
        StepOutcome::TemplateMissing(format!(
            "No CCD template found for residue '{}'",
            residue.name
        ))
    };

    let has_errors = matches!(step1_outcome, StepOutcome::Failed(_))
        || matches!(
            step2_outcome,
            StepOutcome::Failed(_) | StepOutcome::TemplateMissing(_)
        );
    let has_modifications =
        combined_report.added_hydrogens > 0 || combined_report.removed_hydrogens > 0;

    if !has_errors {
        // Both steps succeeded without any errors (or Step 1 was not applicable and Step 2 succeeded).
        // Even if 0 hydrogens were added/removed (already fully hydrogenated), record into residue_reports.
        report.record_residue(res_path, combined_report);
    } else if has_modifications {
        // At least one step encountered an error or missing template, but another step succeeded and modified the residue.
        // 1. Record the modifications into residue_reports:
        report.record_residue(res_path.clone(), combined_report);

        // 2. Record all encountered errors into step_errors so partial modifications are never silenced.
        // Missing template on a modified residue is also recorded here (unlike unmodified ligands).
        // Since we record at most once per step here, no linear scan or deduplication is needed (O(1)).
        if let StepOutcome::Failed(err) = step1_outcome {
            report.record_error(res_path.clone(), err);
        }
        match step2_outcome {
            StepOutcome::Failed(err) | StepOutcome::TemplateMissing(err) => {
                report.record_error(res_path, err);
            }
            _ => {}
        }
    } else {
        // No modifications occurred AND errors/missing templates occurred: completely skipped.
        // 1. Record genuine processing failures in step_errors.
        // Missing templates on completely unmodified components are intentionally excluded
        // from step_errors to avoid noise for standard uncataloged ligands/waters.
        if let StepOutcome::Failed(ref err) = step1_outcome {
            report.record_error(res_path.clone(), err.clone());
        }
        if let StepOutcome::Failed(ref err) = step2_outcome {
            report.record_error(res_path.clone(), err.clone());
        }

        // 2. Record in skipped_residues with the primary skip/failure reason (Step 2 preferred).
        let reason = match (&step2_outcome, &step1_outcome) {
            (StepOutcome::Failed(err), _) | (StepOutcome::TemplateMissing(err), _) => err.clone(),
            (_, StepOutcome::Failed(err)) => err.clone(),
            _ => "No hydrogens added or removed".to_string(),
        };
        report.record_skipped(res_path, reason);
    }

    Ok(())
}

/// Hydrogenates a polymer chain group (where child subgroups are residues/components).
fn hydrogenate_chain(
    chain: &mut AtomGroup,
    db: &CcdTemplateDb,
    his_mode: &HisHydrogenationMode,
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

            hydrogenate_single_residue(residue, prev_ref, db, his_mode, report)?;

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

/// Recursively traverses an [`AtomGroup`] hierarchy and applies hydrogen addition
/// using the default histidine protonation mode (`"auto"`).
///
/// Automatically identifies:
/// - Chain-like groups: Groups whose children contain atoms (leaf residue/component groups).
/// - Single component groups: Groups containing atoms directly.
/// - Model/Root container groups: Groups containing subgroups that are further traversed.
pub fn hydrogenate_atomgroup(
    group: &mut AtomGroup,
    db: &CcdTemplateDb,
) -> Result<OverallHydrogenationReport> {
    hydrogenate_atomgroup_with_his(group, db, DEFAULT_HIS_MODE)
}

/// Recursively traverses an [`AtomGroup`] hierarchy and applies hydrogen addition
/// with an explicit histidine mode (`his`).
///
/// Supported modes:
/// - `"auto"` (default): selects `HD1`/`HE2` using [`his_states`] with [`HisStateOptions::default`]:
///   - [`HisState::Hid`]: `HD1` only
///   - [`HisState::Hie`]: `HE2` only
///   - [`HisState::Hip`]: both `HD1` and `HE2`
///   - [`HisState::MetalBridge`]: neither `HD1` nor `HE2`
/// - `"hip"`: adds both `HD1` and `HE2` to histidine residues (traditional full protonation behavior)
///
/// Any other value for `his` returns a [`BridgeError::ValueError`]. Existing hydrogens
/// in the input are never removed, and residue names remain unchanged.
pub fn hydrogenate_atomgroup_with_his(
    group: &mut AtomGroup,
    db: &CcdTemplateDb,
    his: &str,
) -> Result<OverallHydrogenationReport> {
    let his_mode = match his {
        "auto" => {
            let his_report = his_states(group, &HisStateOptions::default());
            let sites_map = his_report
                .sites
                .into_iter()
                .map(|s| (s.path, s.state))
                .collect::<HashMap<_, _>>();
            HisHydrogenationMode::Auto(sites_map)
        }
        "hip" => HisHydrogenationMode::Hip,
        other => {
            return Err(BridgeError::value_error(
                "his",
                format!("unsupported his mode '{other}'; expected 'auto' or 'hip'"),
            ));
        }
    };

    let mut report = OverallHydrogenationReport::new();
    traverse_and_hydrogenate(group, db, &his_mode, &mut report)?;
    Ok(report)
}

fn traverse_and_hydrogenate(
    group: &mut AtomGroup,
    db: &CcdTemplateDb,
    his_mode: &HisHydrogenationMode,
    report: &mut OverallHydrogenationReport,
) -> Result<()> {
    if group.get_number_of_groups() == 0 {
        // Leaf group with direct atoms: single residue or component
        if group.get_number_of_atoms() > 0 {
            hydrogenate_single_residue(group, None, db, his_mode, report)?;
        }
        return Ok(());
    }

    // Check if the children of this group are leaf groups (i.e. this group acts as a chain)
    let all_children_are_leaves = group
        .groups()
        .all(|(_, child)| child.get_number_of_groups() == 0);

    if all_children_are_leaves {
        hydrogenate_chain(group, db, his_mode, report)?;
    } else {
        // Container group (e.g. Root, Model): sort child keys and recurse
        let mut child_keys = group.get_group_list();
        sort_nicely(&mut child_keys);

        for key in &child_keys {
            if let Some(child) = group.get_group_mut(key) {
                traverse_and_hydrogenate(child, db, his_mode, report)?;
            }
        }
    }

    Ok(())
}

impl AtomGroup {
    /// Adds missing hydrogens across this entire structure using standard polymer
    /// connectivity, CCD idealized templates, and the default histidine mode (`"auto"`).
    ///
    /// Delegates to [`crate::orchestrator::hydrogenate_atomgroup`].
    pub fn add_missing_hydrogens(
        &mut self,
        db: &CcdTemplateDb,
    ) -> Result<OverallHydrogenationReport> {
        self.add_missing_hydrogens_with_his(db, DEFAULT_HIS_MODE)
    }

    /// Adds missing hydrogens across this entire structure with an explicit histidine mode (`his`).
    ///
    /// - `"auto"` (default): selects `HD1`/`HE2` using [`his_states`] with [`HisStateOptions::default`].
    /// - `"hip"`: adds both `HD1` and `HE2` to histidine residues (traditional full protonation behavior).
    ///
    /// Passing any other value returns [`BridgeError::ValueError`].
    ///
    /// Delegates to [`crate::orchestrator::hydrogenate_atomgroup_with_his`].
    pub fn add_missing_hydrogens_with_his(
        &mut self,
        db: &CcdTemplateDb,
        his: &str,
    ) -> Result<OverallHydrogenationReport> {
        hydrogenate_atomgroup_with_his(self, db, his)
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;
    use crate::atom::Atom;
    use crate::format::pdb::Pdb;
    use crate::his_state::HisDecidedBy;
    use crate::position::Position;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
    }

    fn make_atom(symbol: &str, name: &str, xyz: Position) -> Atom {
        let mut atom = Atom::new_with_pos(symbol, xyz).expect("valid element symbol");
        atom.name = name.to_string();
        atom
    }

    /// Builds a single histidine residue containing all non-OXT heavy atoms from the
    /// embedded `"HIS"` CCD template at their idealized coordinates.
    fn make_his_heavy_from_ccd(comp_id: &str) -> AtomGroup {
        let db = CcdTemplateDb::global();
        let tpl = db.lookup("HIS").expect("HIS template must exist");
        let mut res = AtomGroup::with_name(comp_id);
        for ccd_atom in &tpl.atoms {
            if !ccd_atom.is_hydrogen() && ccd_atom.name != "OXT" {
                let (x, y, z) = ccd_atom.ideal_xyz.expect("ideal_xyz must exist");
                res.set_atom(
                    &ccd_atom.name,
                    make_atom(&ccd_atom.element, &ccd_atom.name, Position::new(x, y, z)),
                );
            }
        }
        res
    }

    fn strip_all_hydrogens(group: &mut AtomGroup) {
        let h_keys: Vec<String> = group
            .atoms()
            .filter(|(_, a)| a.atomic_number() == 1 || a.symbol() == Ok("H"))
            .map(|(k, _)| k.clone())
            .collect();
        for k in h_keys {
            group.remove_atom(&k);
        }
        for (_, sub) in group.groups_mut() {
            strip_all_hydrogens(sub);
        }
    }

    fn count_all_hydrogens(group: &AtomGroup) -> usize {
        let direct_h = group
            .atoms()
            .filter(|(_, a)| a.atomic_number() == 1 || a.symbol() == Ok("H"))
            .count();
        let sub_h: usize = group
            .groups()
            .map(|(_, sub)| count_all_hydrogens(sub))
            .sum();
        direct_h + sub_h
    }

    // ① HID に HD1 だけ付く
    #[test]
    fn test_01_hid_attaches_only_hd1() {
        let db = CcdTemplateDb::global();
        let opts = HisStateOptions::default();

        // Case A: Isolated HIS -> Default -> Hid
        let mut his = make_his_heavy_from_ccd("HIS");
        let state_rep = his.his_states(&opts);
        assert_eq!(state_rep.sites.len(), 1);
        assert_eq!(state_rep.sites[0].state, HisState::Hid);
        assert_eq!(state_rep.sites[0].decided_by, HisDecidedBy::Default);

        let report = his.add_missing_hydrogens(db).expect("hydrogenate Hid");
        assert!(his.has_atom("HD1"), "HID must have HD1 attached");
        assert!(!his.has_atom("HE2"), "HID must NOT have HE2 attached");
        assert_eq!(his.name, "HIS", "residue name must remain HIS");
        // N-term H1, H2, H3 (3) + HA, HB2, HB3, HD1, HD2, HE1 (6) = 9
        assert_eq!(report.total_added_hydrogens, 9);

        // Case B: ZN coordinated to NE2 -> Metal -> Hid
        let his2 = make_his_heavy_from_ccd("HIS");
        let nd1_pos = his2.get_atom("ND1").unwrap().xyz;
        let ne2_pos = his2.get_atom("NE2").unwrap().xyz;
        let dir_ne2_out = (ne2_pos - nd1_pos) * (1.0 / (ne2_pos - nd1_pos).length());
        let zn_pos = ne2_pos + dir_ne2_out * 2.0;

        let mut zn_res = AtomGroup::with_name("ZN");
        zn_res.set_atom("ZN", make_atom("Zn", "ZN", zn_pos));

        let mut chain = AtomGroup::with_name("A");
        chain.set_group("1", his2);
        chain.set_group("2", zn_res);

        let state_rep2 = chain.his_states(&opts);
        assert_eq!(state_rep2.sites[0].state, HisState::Hid);
        assert_eq!(state_rep2.sites[0].decided_by, HisDecidedBy::Metal);

        let _ = chain
            .add_missing_hydrogens(db)
            .expect("hydrogenate Hid with metal");
        let res1 = chain.get_group("1").unwrap();
        assert!(res1.has_atom("HD1"));
        assert!(!res1.has_atom("HE2"));
        assert_eq!(res1.name, "HIS");
    }

    // ② HIE に HE2 だけ付く
    #[test]
    fn test_02_hie_attaches_only_he2() {
        let db = CcdTemplateDb::global();
        let opts = HisStateOptions::default();

        let his = make_his_heavy_from_ccd("HIS");
        let nd1_pos = his.get_atom("ND1").unwrap().xyz;
        let ne2_pos = his.get_atom("NE2").unwrap().xyz;
        let dir_nd1_out = (nd1_pos - ne2_pos) * (1.0 / (nd1_pos - ne2_pos).length());
        // Place ZN 2.0 Å from ND1 away from NE2 (> 2.5 Å from NE2) -> Metal -> Hie
        let zn_pos = nd1_pos + dir_nd1_out * 2.0;

        let mut zn_res = AtomGroup::with_name("ZN");
        zn_res.set_atom("ZN", make_atom("Zn", "ZN", zn_pos));

        let mut chain = AtomGroup::with_name("A");
        chain.set_group("1", his);
        chain.set_group("2", zn_res);

        let state_rep = chain.his_states(&opts);
        assert_eq!(state_rep.sites.len(), 1);
        assert_eq!(state_rep.sites[0].state, HisState::Hie);
        assert_eq!(state_rep.sites[0].decided_by, HisDecidedBy::Metal);

        let report = chain.add_missing_hydrogens(db).expect("hydrogenate Hie");
        let res1 = chain.get_group("1").unwrap();
        assert!(!res1.has_atom("HD1"), "HIE must NOT have HD1 attached");
        assert!(res1.has_atom("HE2"), "HIE must have HE2 attached");
        assert_eq!(res1.name, "HIS", "residue name must remain HIS");
        // N-term H1, H2, H3 (3) + HA, HB2, HB3, HD2, HE1, HE2 (6) = 9
        assert_eq!(report.total_added_hydrogens, 9);
    }

    // ③ HIP に両方付く
    #[test]
    fn test_03_hip_attaches_both_hd1_and_he2() {
        let db = CcdTemplateDb::global();
        let opts = HisStateOptions::default();

        let his = make_his_heavy_from_ccd("HIS");
        let nd1_pos = his.get_atom("ND1").unwrap().xyz;
        let ne2_pos = his.get_atom("NE2").unwrap().xyz;
        let dir_nd1_out = (nd1_pos - ne2_pos) * (1.0 / (nd1_pos - ne2_pos).length());
        // Place ASP OD1 at 3.0 Å from ND1 (<= 3.5 Å salt bridge) -> SaltBridge -> Hip
        let od1_pos = nd1_pos + dir_nd1_out * 3.0;

        let mut asp_res = AtomGroup::with_name("ASP");
        asp_res.set_atom("OD1", make_atom("O", "OD1", od1_pos));

        let mut chain = AtomGroup::with_name("A");
        chain.set_group("1", his);
        chain.set_group("2", asp_res);

        let state_rep = chain.his_states(&opts);
        assert_eq!(state_rep.sites.len(), 1);
        assert_eq!(state_rep.sites[0].state, HisState::Hip);
        assert_eq!(state_rep.sites[0].decided_by, HisDecidedBy::SaltBridge);

        let report = chain.add_missing_hydrogens(db).expect("hydrogenate Hip");
        let res1 = chain.get_group("1").unwrap();
        assert!(res1.has_atom("HD1"), "HIP must have HD1 attached");
        assert!(res1.has_atom("HE2"), "HIP must have HE2 attached");
        assert_eq!(res1.name, "HIS", "residue name must remain HIS");
        // N-term H1, H2, H3 (3) + HA, HB2, HB3, HD1, HD2, HE1, HE2 (7) = 10
        assert_eq!(report.residue_reports["/1/"].added_hydrogens, 10);
    }

    // ④ MetalBridge は水素なし
    #[test]
    fn test_04_metal_bridge_attaches_no_ring_nitrogen_hydrogens() {
        let db = CcdTemplateDb::global();
        let opts = HisStateOptions::default();

        let his = make_his_heavy_from_ccd("HIS");
        let nd1_pos = his.get_atom("ND1").unwrap().xyz;
        let ne2_pos = his.get_atom("NE2").unwrap().xyz;
        let dir_nd1_out = (nd1_pos - ne2_pos) * (1.0 / (nd1_pos - ne2_pos).length());
        let dir_ne2_out = (ne2_pos - nd1_pos) * (1.0 / (ne2_pos - nd1_pos).length());

        let mut zn1 = AtomGroup::with_name("ZN");
        zn1.set_atom("ZN", make_atom("Zn", "ZN", nd1_pos + dir_nd1_out * 2.0));
        let mut zn2 = AtomGroup::with_name("ZN");
        zn2.set_atom("ZN", make_atom("Zn", "ZN", ne2_pos + dir_ne2_out * 2.0));

        let mut chain = AtomGroup::with_name("A");
        chain.set_group("1", his);
        chain.set_group("2", zn1);
        chain.set_group("3", zn2);

        let state_rep = chain.his_states(&opts);
        assert_eq!(state_rep.sites.len(), 1);
        assert_eq!(state_rep.sites[0].state, HisState::MetalBridge);
        assert_eq!(state_rep.sites[0].decided_by, HisDecidedBy::Metal);

        let report = chain
            .add_missing_hydrogens(db)
            .expect("hydrogenate MetalBridge");
        let res1 = chain.get_group("1").unwrap();
        assert!(
            !res1.has_atom("HD1"),
            "MetalBridge must NOT have HD1 attached"
        );
        assert!(
            !res1.has_atom("HE2"),
            "MetalBridge must NOT have HE2 attached"
        );
        assert!(res1.has_atom("HA"));
        assert!(res1.has_atom("HD2"));
        assert!(res1.has_atom("HE1"));
        assert_eq!(res1.name, "HIS", "residue name must remain HIS");
        // N-term H1, H2, H3 (3) + HA, HB2, HB3, HD2, HE1 (5) = 8
        assert_eq!(report.total_added_hydrogens, 8);
    }

    // ⑤ 入力の HE2 を消さない
    #[test]
    fn test_05_preserves_existing_input_he2() {
        let db = CcdTemplateDb::global();
        let opts = HisStateOptions::default();

        let mut his = make_his_heavy_from_ccd("HIS");
        let custom_he2_pos = Position::new(12.345, -6.789, 0.543);
        his.set_atom("HE2", make_atom("H", "HE2", custom_he2_pos));

        // Existing HE2 makes his_states() decide Hie (by Hydrogens)
        let state_rep = his.his_states(&opts);
        assert_eq!(state_rep.sites.len(), 1);
        assert_eq!(state_rep.sites[0].state, HisState::Hie);
        assert_eq!(state_rep.sites[0].decided_by, HisDecidedBy::Hydrogens);

        let report = his
            .add_missing_hydrogens(db)
            .expect("hydrogenate HIS with existing HE2");
        assert!(his.has_atom("HE2"), "existing HE2 must not be removed");
        assert_eq!(
            his.get_atom("HE2").unwrap().xyz,
            custom_he2_pos,
            "existing HE2 coordinates must be preserved"
        );
        assert!(!his.has_atom("HD1"), "HIE must not add HD1");
        assert_eq!(report.total_removed_hydrogens, 0);
        // N-term H1, H2, H3 (3) + HA, HB2, HB3, HD2, HE1 (5) = 8 newly added
        assert_eq!(report.total_added_hydrogens, 8);
        assert_eq!(his.name, "HIS");

        // Also verify with his="hip" that existing HE2 is preserved and HD1 is added
        let mut his_hip = make_his_heavy_from_ccd("HIS");
        his_hip.set_atom("HE2", make_atom("H", "HE2", custom_he2_pos));
        let rep_hip = his_hip
            .add_missing_hydrogens_with_his(db, "hip")
            .expect("hydrogenate with his=hip");
        assert!(his_hip.has_atom("HE2"), "existing HE2 must not be removed");
        assert_eq!(
            his_hip.get_atom("HE2").unwrap().xyz,
            custom_he2_pos,
            "existing HE2 coordinates must be preserved"
        );
        assert!(his_hip.has_atom("HD1"), "HIP must attach HD1");
        assert_eq!(rep_hip.total_removed_hydrogens, 0);
        // N-term H1, H2, H3 (3) + HA, HB2, HB3, HD1, HD2, HE1 (6) = 9
        assert_eq!(rep_hip.total_added_hydrogens, 9);
    }

    // ⑥ his="hip" で両方(HD1・HE2)を付加(従来の動作)、および "auto"・"hip" 以外の値でエラー
    #[test]
    fn test_06_his_hip_and_invalid_modes() {
        let db = CcdTemplateDb::global();

        // 1. his="hip" attaches both HD1 and HE2
        let mut his_hip = make_his_heavy_from_ccd("HIS");
        let report_hip = his_hip
            .add_missing_hydrogens_with_his(db, "hip")
            .expect("hydrogenate with his=hip");
        assert!(his_hip.has_atom("HD1"), "his=\"hip\" must attach HD1");
        assert!(his_hip.has_atom("HE2"), "his=\"hip\" must attach HE2");
        assert!(his_hip.has_atom("HA"));
        assert!(his_hip.has_atom("HD2"));
        assert!(his_hip.has_atom("HE1"));
        assert_eq!(his_hip.name, "HIS");
        // N-term H1, H2, H3 (3) + HA, HB2, HB3, HD1, HD2, HE1, HE2 (7) = 10
        assert_eq!(report_hip.total_added_hydrogens, 10);

        // 2. Invalid modes ("none", "off", "", etc.) must return an error
        let mut his_invalid = make_his_heavy_from_ccd("HIS");
        for invalid_mode in ["none", "off", "hid", "hie", ""] {
            let res = his_invalid.add_missing_hydrogens_with_his(db, invalid_mode);
            assert!(
                matches!(res, Err(BridgeError::ValueError { .. })),
                "mode '{invalid_mode}' must return BridgeError::ValueError"
            );
        }

        // 3. If a histidine residue has no entry in the Auto sites map, an error is returned
        let mut his_unmapped = make_his_heavy_from_ccd("HIS");
        let mut report = OverallHydrogenationReport::new();
        let empty_auto_mode = HisHydrogenationMode::Auto(HashMap::new());
        let err = traverse_and_hydrogenate(&mut his_unmapped, db, &empty_auto_mode, &mut report);
        assert!(
            matches!(err, Err(BridgeError::ValueError { .. })),
            "missing path in Auto sites map must return BridgeError::ValueError"
        );
    }

    // ⑦ 1hls で、付いた水素の数と残基の状態を報告する
    #[test]
    fn test_07_1hls_hydrogen_count_and_his_states() {
        let pdb_path = test_data_dir().join("1hls.pdb");
        let pdb = Pdb::from_file(&pdb_path, None).expect("load 1hls.pdb");
        let db = CcdTemplateDb::global();
        let his_opts = HisStateOptions::default();

        // 1. Native 1hls.pdb (starts with 379 hydrogens; B/5, B/10, B/16 already have HD1)
        let mut ag_native = pdb.get_atomgroup(None, None).expect("1hls atomgroup");
        assert_eq!(count_all_hydrogens(&ag_native), 379);

        let his_rep_native = ag_native.his_states(&his_opts);
        assert_eq!(his_rep_native.sites.len(), 3);
        for site in &his_rep_native.sites {
            assert_eq!(site.comp_id, "HIS");
            assert_eq!(site.state, HisState::Hid);
            assert_eq!(site.decided_by, HisDecidedBy::Hydrogens);
        }

        let rep_native_auto = ag_native
            .add_missing_hydrogens(db)
            .expect("hydrogenate native 1hls (auto)");
        println!(
            "1hls native (his=auto): added={}, hydrogenated_residues={}, final_h={}",
            rep_native_auto.total_added_hydrogens,
            rep_native_auto.hydrogenated_residues,
            count_all_hydrogens(&ag_native)
        );
        assert_eq!(rep_native_auto.total_added_hydrogens, 12);
        assert_eq!(rep_native_auto.hydrogenated_residues, 12);
        assert_eq!(count_all_hydrogens(&ag_native), 391);

        // Compare native 1hls with his="hip" (adds missing HE2 to the 3 HIS residues: 12 + 3 = 15 added, 394 total)
        let mut ag_native_hip = pdb.get_atomgroup(None, None).expect("1hls atomgroup");
        let rep_native_hip = ag_native_hip
            .add_missing_hydrogens_with_his(db, "hip")
            .expect("hydrogenate native 1hls (hip)");
        println!(
            "1hls native (his=hip): added={}, hydrogenated_residues={}, final_h={}",
            rep_native_hip.total_added_hydrogens,
            rep_native_hip.hydrogenated_residues,
            count_all_hydrogens(&ag_native_hip)
        );
        assert_eq!(rep_native_hip.total_added_hydrogens, 15);
        assert_eq!(rep_native_hip.hydrogenated_residues, 15);
        assert_eq!(count_all_hydrogens(&ag_native_hip), 394);

        // 2. Heavy-atom-only 1hls.pdb (all hydrogens stripped first -> 0 initial hydrogens)
        let mut ag_stripped_auto = pdb.get_atomgroup(None, None).expect("1hls atomgroup");
        strip_all_hydrogens(&mut ag_stripped_auto);
        assert_eq!(count_all_hydrogens(&ag_stripped_auto), 0);

        let his_rep_stripped = ag_stripped_auto.his_states(&his_opts);
        assert_eq!(his_rep_stripped.sites.len(), 3);
        for site in &his_rep_stripped.sites {
            println!(
                "1hls stripped site: path={}, comp_id={}, state={:?}, decided_by={:?}",
                site.path, site.comp_id, site.state, site.decided_by
            );
            assert_eq!(site.comp_id, "HIS");
            assert_eq!(site.state, HisState::Hid);
            assert_eq!(site.decided_by, HisDecidedBy::Burial);
        }

        let rep_stripped_auto = ag_stripped_auto
            .add_missing_hydrogens(db)
            .expect("hydrogenate stripped 1hls (auto)");
        println!(
            "1hls stripped (his=auto): added={}, hydrogenated_residues={}, final_h={}",
            rep_stripped_auto.total_added_hydrogens,
            rep_stripped_auto.hydrogenated_residues,
            count_all_hydrogens(&ag_stripped_auto)
        );
        assert_eq!(rep_stripped_auto.total_added_hydrogens, 391);
        assert_eq!(rep_stripped_auto.hydrogenated_residues, 51);
        assert_eq!(count_all_hydrogens(&ag_stripped_auto), 391);

        let chain_b_auto = ag_stripped_auto
            .get_group("model_1")
            .unwrap()
            .get_group("B")
            .unwrap();
        for res_key in ["5", "10", "16"] {
            let his_res = chain_b_auto.get_group(res_key).unwrap();
            assert_eq!(his_res.name, "HIS");
            assert!(his_res.has_atom("HD1"));
            assert!(!his_res.has_atom("HE2"));
            let path = format!("/model_1/B/{res_key}/");
            assert_eq!(rep_stripped_auto.residue_reports[&path].added_hydrogens, 7);
        }

        // 3. Compare with his="hip" on stripped 1hls.pdb (conventional behavior -> 394 total)
        let mut ag_stripped_hip = pdb.get_atomgroup(None, None).expect("1hls atomgroup");
        strip_all_hydrogens(&mut ag_stripped_hip);
        let rep_stripped_hip = ag_stripped_hip
            .add_missing_hydrogens_with_his(db, "hip")
            .expect("hydrogenate stripped 1hls (hip)");
        println!(
            "1hls stripped (his=hip): added={}, hydrogenated_residues={}, final_h={}",
            rep_stripped_hip.total_added_hydrogens,
            rep_stripped_hip.hydrogenated_residues,
            count_all_hydrogens(&ag_stripped_hip)
        );
        assert_eq!(rep_stripped_hip.total_added_hydrogens, 394);
        assert_eq!(rep_stripped_hip.hydrogenated_residues, 51);
        assert_eq!(count_all_hydrogens(&ag_stripped_hip), 394);

        let chain_b_hip = ag_stripped_hip
            .get_group("model_1")
            .unwrap()
            .get_group("B")
            .unwrap();
        for res_key in ["5", "10", "16"] {
            let his_res = chain_b_hip.get_group(res_key).unwrap();
            assert_eq!(his_res.name, "HIS");
            assert!(his_res.has_atom("HD1"));
            assert!(his_res.has_atom("HE2"));
            let path = format!("/model_1/B/{res_key}/");
            assert_eq!(rep_stripped_hip.residue_reports[&path].added_hydrogens, 8);
        }
    }
}
