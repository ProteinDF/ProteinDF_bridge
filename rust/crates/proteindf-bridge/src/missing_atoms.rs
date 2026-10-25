// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Missing-atom inspection for macromolecular structures using CCD templates.
//!
//! Compares actual atoms in each residue/component against expected atoms derived
//! from the Chemical Component Dictionary ([`CcdTemplateDb`]), accounting for
//! polymer terminal and protonation-dependent optional atoms.

use std::collections::{BTreeMap, HashSet};

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::brd::sort_nicely;
use crate::ccd_templates::{CcdAtom, CcdBondTemplate, CcdTemplateDb};
use crate::position::Position;

/// Minimum plausible $O3'_{prev} - P$ bond distance in Angstroms for nucleic acid backbone continuity.
const MIN_PHOSPHODIESTER_O3_P_DISTANCE: f64 = 1.3;

/// Maximum plausible $O3'_{prev} - P$ bond distance in Angstroms for nucleic acid backbone continuity.
const MAX_PHOSPHODIESTER_O3_P_DISTANCE: f64 = 2.0;

/// Maximum $SG - SG$ distance in Angstroms to treat two cysteine residues as disulfide-bonded.
pub(crate) const MAX_DISULFIDE_SG_SG_DISTANCE: f64 = 2.5;

/// Options for missing-atom inspection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MissingAtomsOptions {
    /// Whether to include hydrogen atoms in the check (default: `false`).
    pub include_hydrogens: bool,
}

/// Completion status of a single residue/component against its CCD template.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResidueStatus {
    /// All expected atoms are present (even if unexpected atoms exist).
    Complete,
    /// One or more expected atoms are missing.
    Incomplete,
    /// No CCD template was found in the database for this component.
    NoTemplate,
}

/// Per-residue result of missing-atom inspection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResidueAtomCheck {
    /// Hierarchical path of the residue group (e.g. `"/model_1/A/1/"`).
    pub path: String,
    /// Chemical component identifier (residue group's `name`, e.g. `"ALA"`).
    pub comp_id: String,
    /// Inspection status (`Complete`, `Incomplete`, or `NoTemplate`).
    pub status: ResidueStatus,
    /// Number of expected atoms (after filtering by element and removing optional atoms).
    pub expected: usize,
    /// Number of expected atoms that are present in the residue.
    pub present: usize,
    /// Missing expected atom names, in template order.
    pub missing: Vec<String>,
    /// Atom names present in the residue that are neither expected nor optional,
    /// in tree order with duplicates deduplicated (first appearance kept).
    pub unexpected: Vec<String>,
}

/// Structure-wide report of missing-atom inspection across residues.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct MissingAtomsReport {
    /// Inspection results for all visited residues in traversal order.
    pub residues: Vec<ResidueAtomCheck>,
}

impl MissingAtomsReport {
    /// Returns all residues with status [`ResidueStatus::Incomplete`].
    pub fn incomplete(&self) -> Vec<&ResidueAtomCheck> {
        self.residues
            .iter()
            .filter(|r| r.status == ResidueStatus::Incomplete)
            .collect()
    }

    /// Returns all residues with status [`ResidueStatus::NoTemplate`].
    pub fn no_template(&self) -> Vec<&ResidueAtomCheck> {
        self.residues
            .iter()
            .filter(|r| r.status == ResidueStatus::NoTemplate)
            .collect()
    }

    /// Returns the total number of missing atoms across all residues.
    pub fn total_missing(&self) -> usize {
        self.residues.iter().map(|r| r.missing.len()).sum()
    }

    /// Returns `true` if there are no [`ResidueStatus::Incomplete`] residues.
    pub fn is_complete(&self) -> bool {
        !self
            .residues
            .iter()
            .any(|r| r.status == ResidueStatus::Incomplete)
    }

    /// Returns a map from component ID to residue count for all [`ResidueStatus::NoTemplate`] residues.
    pub fn no_template_components(&self) -> BTreeMap<String, usize> {
        let mut counts = BTreeMap::new();
        for r in &self.residues {
            if r.status == ResidueStatus::NoTemplate {
                *counts.entry(r.comp_id.clone()).or_insert(0) += 1;
            }
        }
        counts
    }
}

/// Finds an atom in `residue` strictly by `Atom.name` (ignoring storage keys).
fn find_atom_by_name<'a>(residue: &'a AtomGroup, name: &str) -> Option<&'a Atom> {
    residue.atoms().map(|(_, a)| a).find(|a| a.name == name)
}

/// Checks if `residue` contains an atom with the given `Atom.name` (ignoring storage keys).
fn has_atom_by_name(residue: &AtomGroup, name: &str) -> bool {
    find_atom_by_name(residue, name).is_some()
}

/// Checks if a residue/template represents a protein component (possesses `N` and `CA`).
fn is_protein_component(residue: &AtomGroup, template: &CcdBondTemplate) -> bool {
    (template.get_atom("N").is_some() && template.get_atom("CA").is_some())
        || (has_atom_by_name(residue, "N") && has_atom_by_name(residue, "CA"))
}

/// Checks if a template represents a nucleic acid component (possesses `O3'` and `C1'`).
fn is_nucleic_template(template: &CcdBondTemplate) -> bool {
    template.get_atom("O3'").is_some() && template.get_atom("C1'").is_some()
}

/// Checks if a template represents a cysteine-like thiol component (possesses `SG` and `HG`).
fn is_cysteine_template(template: &CcdBondTemplate) -> bool {
    template.get_atom("SG").is_some() && template.get_atom("HG").is_some()
}

/// Given `(residue_id, sg_position)` entries collected within a single scope (model or chain),
/// returns the set of `residue_id`s whose `SG` atom lies within [`MAX_DISULFIDE_SG_SG_DISTANCE`]
/// (2.5 Å) of another residue's `SG`.
pub(crate) fn find_disulfide_bonded_ids(sgs: &[(usize, Position)]) -> HashSet<usize> {
    let mut bonded = HashSet::new();
    for (i, (id1, pos1)) in sgs.iter().enumerate() {
        for (id2, pos2) in sgs.iter().skip(i + 1) {
            if id1 != id2 && pos1.distance_from(pos2) <= MAX_DISULFIDE_SG_SG_DISTANCE {
                bonded.insert(*id1);
                bonded.insert(*id2);
            }
        }
    }
    bonded
}

/// Determines whether a nucleic acid residue is at the 5' terminus.
///
/// A residue is considered 5'-terminal if:
/// - There is no preceding residue in the same chain (`prev_residue` is `None`), or
/// - The distance between the preceding residue's `O3'` and this residue's `P` is outside
///   `1.3..=2.0` Å.
///
/// If `prev_residue` is `Some` and this residue lacks a `P` atom, the distance cannot be
/// computed and the residue is treated as non-terminal (`false`).
pub(crate) fn is_five_prime_terminal(
    residue: &AtomGroup,
    prev_residue: Option<&AtomGroup>,
) -> bool {
    match prev_residue {
        None => true,
        Some(prev) => {
            let Some(p_curr) = find_atom_by_name(residue, "P") else {
                return false;
            };
            let Some(o3_prev) = find_atom_by_name(prev, "O3'") else {
                return true;
            };
            let dist = o3_prev.xyz.distance_from(&p_curr.xyz);
            !(MIN_PHOSPHODIESTER_O3_P_DISTANCE..=MAX_PHOSPHODIESTER_O3_P_DISTANCE).contains(&dist)
        }
    }
}

/// Returns the set of optional atom names for `residue` under `template` and `options`.
///
/// Optional atoms are neither required (their absence is not reported in `missing`) nor
/// rejected when present (their presence is not reported in `unexpected`).
pub(crate) fn optional_atom_names(
    residue: &AtomGroup,
    prev_residue: Option<&AtomGroup>,
    is_disulfide_bonded: bool,
    template: &CcdBondTemplate,
    options: &MissingAtomsOptions,
) -> HashSet<&'static str> {
    let mut optional = HashSet::new();

    if is_protein_component(residue, template) {
        optional.insert("OXT");
        if options.include_hydrogens {
            optional.extend(["H", "H1", "H2", "H3", "HXT"]);
        }
    }

    if options.include_hydrogens
        && is_disulfide_bonded
        && is_cysteine_template(template)
        && has_atom_by_name(residue, "SG")
    {
        optional.insert("HG");
    }

    if is_nucleic_template(template) {
        optional.insert("OP3");
        if options.include_hydrogens {
            optional.extend(["HOP2", "HOP3", "HO3'"]);
        }

        if is_five_prime_terminal(residue, prev_residue) {
            let phosphate_atoms = ["P", "OP1", "OP2", "OP3"];
            if phosphate_atoms
                .iter()
                .all(|&name| !has_atom_by_name(residue, name))
            {
                optional.extend(phosphate_atoms);
            }
        }
    }

    optional
}

/// Returns the expected [`CcdAtom`] entries from `template` for `residue` in template order.
///
/// Filters out hydrogens when `options.include_hydrogens` is `false` and removes all optional
/// atoms determined by [`optional_atom_names`].
pub(crate) fn expected_atoms<'a>(
    residue: &AtomGroup,
    prev_residue: Option<&AtomGroup>,
    is_disulfide_bonded: bool,
    template: &'a CcdBondTemplate,
    options: &MissingAtomsOptions,
) -> Vec<&'a CcdAtom> {
    let optional = optional_atom_names(
        residue,
        prev_residue,
        is_disulfide_bonded,
        template,
        options,
    );
    template
        .atoms
        .iter()
        .filter(|a| options.include_hydrogens || !a.is_hydrogen())
        .filter(|a| !optional.contains(a.name.as_str()))
        .collect()
}

/// Checks a single leaf residue against `db`.
fn check_single_residue(
    residue: &AtomGroup,
    prev_residue: Option<&AtomGroup>,
    is_disulfide_bonded: bool,
    db: &CcdTemplateDb,
    options: &MissingAtomsOptions,
) -> ResidueAtomCheck {
    let path = residue.path().to_string();
    let comp_id = residue.name.clone();

    let Some(template) = db.lookup(&comp_id) else {
        return ResidueAtomCheck {
            path,
            comp_id,
            status: ResidueStatus::NoTemplate,
            expected: 0,
            present: 0,
            missing: Vec::new(),
            unexpected: Vec::new(),
        };
    };

    let optional = optional_atom_names(
        residue,
        prev_residue,
        is_disulfide_bonded,
        template,
        options,
    );
    let expected_list = expected_atoms(
        residue,
        prev_residue,
        is_disulfide_bonded,
        template,
        options,
    );
    let expected_set: HashSet<&str> = expected_list.iter().map(|a| a.name.as_str()).collect();

    // Collect atom names from the residue for the relevant elements (preserving tree order)
    let mut tree_atom_names: Vec<&str> = Vec::new();
    let mut present_set: HashSet<&str> = HashSet::new();
    for (_, atom) in residue.atoms() {
        let is_h = atom.atomic_number() == 1;
        if !options.include_hydrogens && is_h {
            continue;
        }
        tree_atom_names.push(atom.name.as_str());
        present_set.insert(atom.name.as_str());
    }

    let mut missing = Vec::new();
    let mut present = 0;
    for exp_atom in &expected_list {
        if present_set.contains(exp_atom.name.as_str()) {
            present += 1;
        } else {
            missing.push(exp_atom.name.clone());
        }
    }

    let mut unexpected = Vec::new();
    let mut seen_unexpected: HashSet<&str> = HashSet::new();
    for name in tree_atom_names {
        if !expected_set.contains(name) && !optional.contains(name) && seen_unexpected.insert(name)
        {
            unexpected.push(name.to_string());
        }
    }

    let status = if missing.is_empty() {
        ResidueStatus::Complete
    } else {
        ResidueStatus::Incomplete
    };

    ResidueAtomCheck {
        path,
        comp_id,
        status,
        expected: expected_list.len(),
        present,
        missing,
        unexpected,
    }
}

/// Entry for a leaf residue collected during a single tree traversal.
struct CollectedResidue<'a> {
    residue: &'a AtomGroup,
    prev_residue: Option<&'a AtomGroup>,
}

/// State accumulated during a single pass over the [`AtomGroup`] tree.
#[derive(Default)]
struct TraversalState<'a> {
    residues: Vec<CollectedResidue<'a>>,
    sgs_by_scope: Vec<Vec<(usize, Position)>>,
}

impl<'a> TraversalState<'a> {
    fn new_scope(&mut self) -> usize {
        let id = self.sgs_by_scope.len();
        self.sgs_by_scope.push(Vec::new());
        id
    }

    fn push_residue(
        &mut self,
        residue: &'a AtomGroup,
        prev_residue: Option<&'a AtomGroup>,
        scope_id: usize,
        collect_sgs: bool,
    ) {
        let res_idx = self.residues.len();
        if collect_sgs {
            if let Some(sg) = find_atom_by_name(residue, "SG") {
                self.sgs_by_scope[scope_id].push((res_idx, sg.xyz));
            }
        }
        self.residues.push(CollectedResidue {
            residue,
            prev_residue,
        });
    }
}

/// Collects all residues in a chain group (where direct children are leaf residue groups).
fn collect_from_chain<'a>(
    chain: &'a AtomGroup,
    scope_id: usize,
    collect_sgs: bool,
    state: &mut TraversalState<'a>,
) {
    let mut res_keys = chain.get_group_list();
    sort_nicely(&mut res_keys);

    let mut prev_residue: Option<&'a AtomGroup> = None;

    for key in &res_keys {
        if let Some(residue) = chain.get_group(key) {
            if residue.get_number_of_atoms() > 0 {
                state.push_residue(residue, prev_residue, scope_id, collect_sgs);
                prev_residue = Some(residue);
            }
        }
    }
}

/// Recursively traverses an [`AtomGroup`] hierarchy in a single pass using the same
/// leaf/chain rules as `orchestrator::traverse_and_hydrogenate`, collecting residues
/// and per-scope `SG` coordinates.
///
/// A disulfide search scope corresponds to a group that contains chains (e.g. a model
/// `/model_N/` in a full PDB/mmCIF tree), or the root of the call itself when called on
/// a single chain or residue.
fn collect_from_group<'a>(
    group: &'a AtomGroup,
    inherited_scope: Option<usize>,
    collect_sgs: bool,
    state: &mut TraversalState<'a>,
) {
    if group.get_number_of_groups() == 0 {
        if group.get_number_of_atoms() > 0 {
            let scope = inherited_scope.unwrap_or_else(|| state.new_scope());
            state.push_residue(group, None, scope, collect_sgs);
        }
        return;
    }

    let all_children_are_leaves = group
        .groups()
        .all(|(_, child)| child.get_number_of_groups() == 0);

    if all_children_are_leaves {
        let scope = inherited_scope.unwrap_or_else(|| state.new_scope());
        collect_from_chain(group, scope, collect_sgs, state);
    } else {
        let mut child_keys = group.get_group_list();
        sort_nicely(&mut child_keys);

        let mut shared_chain_scope = inherited_scope;

        for key in &child_keys {
            if let Some(child) = group.get_group(key) {
                if child.get_number_of_groups() == 0 {
                    if child.get_number_of_atoms() > 0 {
                        let scope = *shared_chain_scope.get_or_insert_with(|| state.new_scope());
                        state.push_residue(child, None, scope, collect_sgs);
                    }
                } else if child.groups().all(|(_, gc)| gc.get_number_of_groups() == 0) {
                    let scope = *shared_chain_scope.get_or_insert_with(|| state.new_scope());
                    collect_from_chain(child, scope, collect_sgs, state);
                } else {
                    collect_from_group(child, inherited_scope, collect_sgs, state);
                }
            }
        }
    }
}

/// Inspects an [`AtomGroup`] hierarchy for missing and unexpected atoms against `db`.
pub fn check_missing_atoms(
    group: &AtomGroup,
    db: &CcdTemplateDb,
    options: &MissingAtomsOptions,
) -> MissingAtomsReport {
    let mut state = TraversalState::default();
    collect_from_group(group, None, options.include_hydrogens, &mut state);

    let mut disulfide_bonded = HashSet::new();
    if options.include_hydrogens {
        for sgs in &state.sgs_by_scope {
            disulfide_bonded.extend(find_disulfide_bonded_ids(sgs));
        }
    }

    let residues = state
        .residues
        .iter()
        .enumerate()
        .map(|(idx, item)| {
            check_single_residue(
                item.residue,
                item.prev_residue,
                disulfide_bonded.contains(&idx),
                db,
                options,
            )
        })
        .collect();

    MissingAtomsReport { residues }
}

impl AtomGroup {
    /// Inspects each residue in this hierarchy against CCD templates in `db`
    /// and reports missing and unexpected atoms without modifying the tree.
    pub fn check_missing_atoms(
        &self,
        db: &CcdTemplateDb,
        options: &MissingAtomsOptions,
    ) -> MissingAtomsReport {
        check_missing_atoms(self, db, options)
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;
    use crate::format::mmcif::SimpleMmcif;
    use crate::position::Position;

    fn make_atom(symbol: &str, name: &str, xyz: Position) -> Atom {
        let mut atom = Atom::new_with_pos(symbol, xyz).expect("valid element symbol");
        atom.name = name.to_string();
        atom
    }

    fn build_residue_from_template(
        comp_id: &str,
        include_hydrogens: bool,
        exclude_names: &[&str],
    ) -> AtomGroup {
        let db = CcdTemplateDb::global();
        let template = db.lookup(comp_id).expect("template should exist");
        let mut res = AtomGroup::with_name(comp_id);
        for (idx, ccd_atom) in template.atoms.iter().enumerate() {
            if !include_hydrogens && ccd_atom.is_hydrogen() {
                continue;
            }
            if exclude_names.contains(&ccd_atom.name.as_str()) {
                continue;
            }
            let (x, y, z) = ccd_atom.ideal_xyz.unwrap_or((0.0, 0.0, 0.0));
            let atom = make_atom(&ccd_atom.element, &ccd_atom.name, Position::new(x, y, z));
            let key = format!("{}_{}", idx + 1, ccd_atom.name);
            res.set_atom(&key, atom);
        }
        res
    }

    #[test]
    fn test_ala_complete_without_oxt_missing_cb_and_hydrogens() {
        let db = CcdTemplateDb::global();
        let default_opts = MissingAtomsOptions::default();
        assert!(!default_opts.include_hydrogens);

        // 1. ALA without OXT is complete in heavy-atom mode
        let ala_no_oxt = build_residue_from_template("ALA", false, &["OXT"]);
        let report = ala_no_oxt.check_missing_atoms(db, &default_opts);
        assert_eq!(report.residues.len(), 1);
        assert!(report.is_complete());
        assert_eq!(report.total_missing(), 0);
        assert!(report.incomplete().is_empty());
        let r = &report.residues[0];
        assert_eq!(r.comp_id, "ALA");
        assert_eq!(r.status, ResidueStatus::Complete);
        assert_eq!(r.expected, 5);
        assert_eq!(r.present, 5);
        assert!(r.missing.is_empty());
        assert!(r.unexpected.is_empty());

        // 2. ALA missing CB is incomplete
        let ala_no_cb = build_residue_from_template("ALA", false, &["CB", "OXT"]);
        let report_cb = ala_no_cb.check_missing_atoms(db, &default_opts);
        assert!(!report_cb.is_complete());
        assert_eq!(report_cb.total_missing(), 1);
        assert_eq!(report_cb.incomplete().len(), 1);
        let r_cb = &report_cb.residues[0];
        assert_eq!(r_cb.status, ResidueStatus::Incomplete);
        assert_eq!(r_cb.expected, 5);
        assert_eq!(r_cb.present, 4);
        assert_eq!(r_cb.missing, vec!["CB"]);
        assert!(r_cb.unexpected.is_empty());

        // 3. Hydrogen mode: HA, HB1..HB3 are missing, while H, H2, HXT are optional and not in missing
        let h_opts = MissingAtomsOptions {
            include_hydrogens: true,
        };
        let report_h = ala_no_oxt.check_missing_atoms(db, &h_opts);
        assert!(!report_h.is_complete());
        let r_h = &report_h.residues[0];
        assert_eq!(r_h.status, ResidueStatus::Incomplete);
        assert_eq!(r_h.expected, 9); // N, CA, C, O, CB, HA, HB1, HB2, HB3
        assert_eq!(r_h.present, 5);
        assert_eq!(r_h.missing, vec!["HA", "HB1", "HB2", "HB3"]);
        assert!(r_h.unexpected.is_empty());

        // Optional hydrogens (H, H1, H2, H3, HXT) and OXT are not reported as unexpected when present
        let mut ala_nterm = build_residue_from_template("ALA", true, &["H", "HXT"]);
        ala_nterm.set_atom("90_H1", make_atom("H", "H1", Position::new(0.0, 0.0, 0.0)));
        ala_nterm.set_atom("91_H3", make_atom("H", "H3", Position::new(0.0, 0.0, 0.0)));
        let report_nterm = ala_nterm.check_missing_atoms(db, &h_opts);
        assert!(report_nterm.is_complete());
        assert!(report_nterm.residues[0].unexpected.is_empty());
    }

    #[test]
    fn test_hoh_heavy_and_hydrogens() {
        let db = CcdTemplateDb::global();
        let hoh = build_residue_from_template("HOH", false, &[]);

        // Heavy-atom mode: O is present -> Complete
        let rep_heavy = hoh.check_missing_atoms(db, &MissingAtomsOptions::default());
        assert!(rep_heavy.is_complete());
        assert_eq!(rep_heavy.residues.len(), 1);
        let r_heavy = &rep_heavy.residues[0];
        assert_eq!(r_heavy.status, ResidueStatus::Complete);
        assert_eq!(r_heavy.expected, 1);
        assert_eq!(r_heavy.present, 1);
        assert!(r_heavy.missing.is_empty());

        // Hydrogen mode: H1 and H2 are missing -> Incomplete
        let rep_h = hoh.check_missing_atoms(
            db,
            &MissingAtomsOptions {
                include_hydrogens: true,
            },
        );
        assert!(!rep_h.is_complete());
        let r_h = &rep_h.residues[0];
        assert_eq!(r_h.status, ResidueStatus::Incomplete);
        assert_eq!(r_h.expected, 3);
        assert_eq!(r_h.present, 1);
        assert_eq!(r_h.missing, vec!["H1", "H2"]);
    }

    #[test]
    fn test_unknown_component_and_custom_ccd_db() {
        let db = CcdTemplateDb::global();
        let opts = MissingAtomsOptions::default();

        let mut lig = AtomGroup::with_name("LIG");
        lig.set_atom("1_C1", make_atom("C", "C1", Position::new(0.0, 0.0, 0.0)));
        lig.set_atom("2_O1", make_atom("O", "O1", Position::new(1.2, 0.0, 0.0)));

        // Without template in global DB -> NoTemplate
        let rep_unknown = lig.check_missing_atoms(db, &opts);
        assert!(rep_unknown.is_complete()); // Incomplete is 0
        assert_eq!(rep_unknown.total_missing(), 0);
        assert_eq!(rep_unknown.no_template().len(), 1);
        assert_eq!(
            rep_unknown.no_template_components(),
            BTreeMap::from([("LIG".to_string(), 1)])
        );
        let r_unk = &rep_unknown.residues[0];
        assert_eq!(r_unk.status, ResidueStatus::NoTemplate);
        assert_eq!(r_unk.expected, 0);
        assert_eq!(r_unk.present, 0);
        assert!(r_unk.missing.is_empty());
        assert!(r_unk.unexpected.is_empty());

        // Add minimal CCD block to a custom DB -> LIG is now inspected
        let minimal_cif = "\
data_comp_LIG
_chem_comp.id LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
LIG C1 C
LIG O1 O
LIG C2 C
LIG H1 H
";
        let cif = SimpleMmcif::from_str(minimal_cif).expect("parse minimal CCD cif");
        let block = cif.get_data_block("comp_LIG").expect("comp_LIG block");
        let template =
            CcdBondTemplate::from_mmcif_block(block, "LIG").expect("build template from block");

        let mut custom_db = db.clone();
        custom_db.insert(template);

        let rep_checked = lig.check_missing_atoms(&custom_db, &opts);
        assert!(!rep_checked.is_complete());
        assert!(rep_checked.no_template().is_empty());
        assert!(rep_checked.no_template_components().is_empty());
        let r_chk = &rep_checked.residues[0];
        assert_eq!(r_chk.status, ResidueStatus::Incomplete);
        assert_eq!(r_chk.expected, 3);
        assert_eq!(r_chk.present, 2);
        assert_eq!(r_chk.missing, vec!["C2"]);
    }

    #[test]
    fn test_da_missing_n9_and_five_prime_terminus() {
        let db = CcdTemplateDb::global();
        let opts = MissingAtomsOptions::default();

        // 1. DA missing N9
        let da_no_n9 = build_residue_from_template("DA", false, &["OP3", "N9"]);
        let rep_n9 = da_no_n9.check_missing_atoms(db, &opts);
        assert!(!rep_n9.is_complete());
        let r_n9 = &rep_n9.residues[0];
        assert_eq!(r_n9.status, ResidueStatus::Incomplete);
        assert_eq!(r_n9.missing, vec!["N9"]);

        // 2. 5'-terminal DA without any of P, OP1, OP2, OP3 -> Complete
        let da_5term_no_phos =
            build_residue_from_template("DA", false, &["P", "OP1", "OP2", "OP3"]);
        let rep_5term = da_5term_no_phos.check_missing_atoms(db, &opts);
        assert!(rep_5term.is_complete());
        let r_5term = &rep_5term.residues[0];
        assert_eq!(r_5term.status, ResidueStatus::Complete);
        assert_eq!(r_5term.expected, 18);
        assert_eq!(r_5term.present, 18);
        assert!(r_5term.missing.is_empty());

        // 3. 5'-terminal DA missing only P (OP1, OP2 present) -> P is reported missing
        let da_5term_no_p = build_residue_from_template("DA", false, &["P", "OP3"]);
        let rep_no_p = da_5term_no_p.check_missing_atoms(db, &opts);
        assert!(!rep_no_p.is_complete());
        let r_no_p = &rep_no_p.residues[0];
        assert_eq!(r_no_p.status, ResidueStatus::Incomplete);
        assert_eq!(r_no_p.missing, vec!["P"]);

        // 4. Internal residue (preceded by another residue in the same chain) lacking P, OP1, OP2, OP3:
        // Since P is absent and a previous residue exists, it is NOT 5'-terminal -> P, OP1, OP2 are missing
        let mut chain = AtomGroup::with_name("A");
        chain.set_group("1", da_5term_no_phos.clone());
        chain.set_group("2", da_5term_no_phos);
        let rep_chain = chain.check_missing_atoms(db, &opts);
        assert_eq!(rep_chain.residues.len(), 2);
        assert_eq!(rep_chain.residues[0].status, ResidueStatus::Complete);
        assert_eq!(rep_chain.residues[1].status, ResidueStatus::Incomplete);
        assert_eq!(rep_chain.residues[1].missing, vec!["P", "OP1", "OP2"]);

        // 5. Distance check for 5'-terminus when P is present:
        // Within 1.3..=2.0 A -> not 5'-terminal; outside 1.3..=2.0 A -> 5'-terminal
        let mut prev_res = AtomGroup::with_name("DA");
        prev_res.set_atom("1_O3'", make_atom("O", "O3'", Position::new(0.0, 0.0, 0.0)));
        let mut curr_connected = AtomGroup::with_name("DA");
        curr_connected.set_atom("1_P", make_atom("P", "P", Position::new(1.6, 0.0, 0.0)));
        assert!(!is_five_prime_terminal(&curr_connected, Some(&prev_res)));

        let mut curr_far = AtomGroup::with_name("DA");
        curr_far.set_atom("1_P", make_atom("P", "P", Position::new(3.5, 0.0, 0.0)));
        assert!(is_five_prime_terminal(&curr_far, Some(&prev_res)));
    }

    #[test]
    fn test_cys_disulfide_hg_optional() {
        let db = CcdTemplateDb::global();
        let h_opts = MissingAtomsOptions {
            include_hydrogens: true,
        };

        let make_cys_with_sg_at = |sg_pos: Position| -> AtomGroup {
            let mut cys =
                build_residue_from_template("CYS", true, &["OXT", "H", "H2", "HXT", "HG"]);
            for (_, atom) in cys.atoms_mut() {
                if atom.name == "SG" {
                    atom.xyz = sg_pos;
                }
            }
            cys
        };

        // 1. Two CYS residues with SG 2.05 A apart in the same model -> HG is optional (not missing)
        let mut model_close = AtomGroup::with_name("model_1");
        let mut chain_a = AtomGroup::with_name("A");
        chain_a.set_group("1", make_cys_with_sg_at(Position::new(0.0, 0.0, 0.0)));
        let mut chain_b = AtomGroup::with_name("B");
        chain_b.set_group("1", make_cys_with_sg_at(Position::new(2.05, 0.0, 0.0)));
        model_close.set_group("A", chain_a);
        model_close.set_group("B", chain_b);

        let rep_close = model_close.check_missing_atoms(db, &h_opts);
        assert_eq!(rep_close.residues.len(), 2);
        assert!(rep_close.is_complete());
        assert_eq!(rep_close.total_missing(), 0);
        for r in &rep_close.residues {
            assert_eq!(r.status, ResidueStatus::Complete);
            assert!(r.missing.is_empty());
            assert!(r.unexpected.is_empty());
        }

        // Same-chain pair 2.05 A apart also treats HG as optional, and present HG is not unexpected
        let mut chain_close = AtomGroup::with_name("A");
        chain_close.set_group("1", make_cys_with_sg_at(Position::new(0.0, 0.0, 0.0)));
        let mut cys_with_hg = make_cys_with_sg_at(Position::new(2.05, 0.0, 0.0));
        cys_with_hg.set_atom("99_HG", make_atom("H", "HG", Position::new(2.5, 0.5, 0.0)));
        chain_close.set_group("2", cys_with_hg);
        let rep_chain_close = chain_close.check_missing_atoms(db, &h_opts);
        assert!(rep_chain_close.is_complete());
        assert!(rep_chain_close.residues[0].missing.is_empty());
        assert!(rep_chain_close.residues[1].unexpected.is_empty());

        // 2. Two CYS residues with SG 5.0 A apart -> HG is missing on both
        let mut chain_far = AtomGroup::with_name("A");
        chain_far.set_group("1", make_cys_with_sg_at(Position::new(0.0, 0.0, 0.0)));
        chain_far.set_group("2", make_cys_with_sg_at(Position::new(5.0, 0.0, 0.0)));
        let rep_far = chain_far.check_missing_atoms(db, &h_opts);
        assert!(!rep_far.is_complete());
        assert_eq!(rep_far.total_missing(), 2);
        for r in &rep_far.residues {
            assert_eq!(r.status, ResidueStatus::Incomplete);
            assert_eq!(r.missing, vec!["HG"]);
        }

        // 3. Nearby SG (2.05 A) in a different model is NOT counted -> HG is missing on both
        let mut root_two_models = AtomGroup::new();
        let mut m1 = AtomGroup::with_name("model_1");
        let mut m1_chain = AtomGroup::with_name("A");
        m1_chain.set_group("1", make_cys_with_sg_at(Position::new(0.0, 0.0, 0.0)));
        m1.set_group("A", m1_chain);

        let mut m2 = AtomGroup::with_name("model_2");
        let mut m2_chain = AtomGroup::with_name("A");
        m2_chain.set_group("1", make_cys_with_sg_at(Position::new(2.05, 0.0, 0.0)));
        m2.set_group("A", m2_chain);

        root_two_models.set_group("model_1", m1);
        root_two_models.set_group("model_2", m2);

        let rep_two_models = root_two_models.check_missing_atoms(db, &h_opts);
        assert_eq!(rep_two_models.residues.len(), 2);
        assert!(!rep_two_models.is_complete());
        assert_eq!(rep_two_models.total_missing(), 2);
        for r in &rep_two_models.residues {
            assert_eq!(r.status, ResidueStatus::Incomplete);
            assert_eq!(r.missing, vec!["HG"]);
        }

        // 4. Calling check_missing_atoms on a single residue alone cannot see its partner -> HG is missing
        let single_cys = model_close.get_group("A").unwrap().get_group("1").unwrap();
        let rep_single = single_cys.check_missing_atoms(db, &h_opts);
        assert_eq!(rep_single.residues.len(), 1);
        assert!(!rep_single.is_complete());
        assert_eq!(rep_single.residues[0].status, ResidueStatus::Incomplete);
        assert_eq!(rep_single.residues[0].missing, vec!["HG"]);
    }

    #[test]
    fn test_unexpected_atoms() {
        let db = CcdTemplateDb::global();
        let opts = MissingAtomsOptions::default();

        let mut ala = build_residue_from_template("ALA", false, &["OXT"]);
        // Optional atom OXT is not unexpected when present
        ala.set_atom(
            "80_OXT",
            make_atom("O", "OXT", Position::new(0.6, 0.4, -1.7)),
        );
        // Atoms neither expected nor optional are reported in tree order
        ala.set_atom("91_ZN", make_atom("Zn", "ZN", Position::new(1.0, 1.0, 1.0)));
        ala.set_atom(
            "92_FOO",
            make_atom("C", "FOO", Position::new(2.0, 2.0, 2.0)),
        );

        let report = ala.check_missing_atoms(db, &opts);
        assert!(report.is_complete());
        let r = &report.residues[0];
        assert_eq!(r.status, ResidueStatus::Complete);
        assert_eq!(r.expected, 5);
        assert_eq!(r.present, 5);
        assert!(r.missing.is_empty());
        assert_eq!(r.unexpected, vec!["ZN", "FOO"]);
    }

    #[test]
    fn test_duplicate_atom_names() {
        let db = CcdTemplateDb::global();
        let opts = MissingAtomsOptions::default();

        let mut ala = build_residue_from_template("ALA", false, &["OXT"]);
        // Duplicate expected atom CA under a second key
        ala.set_atom(
            "90_CA",
            make_atom("C", "CA", Position::new(0.25, 0.41, 0.69)),
        );
        // Duplicate unexpected atom "XX"
        ala.set_atom("91_XX", make_atom("C", "XX", Position::new(1.0, 1.0, 1.0)));
        ala.set_atom("92_YY", make_atom("O", "YY", Position::new(2.0, 2.0, 2.0)));
        ala.set_atom("93_XX", make_atom("C", "XX", Position::new(3.0, 3.0, 3.0)));

        let report = ala.check_missing_atoms(db, &opts);
        assert!(report.is_complete());
        let r = &report.residues[0];
        assert_eq!(r.status, ResidueStatus::Complete);
        assert_eq!(r.expected, 5);
        assert_eq!(r.present, 5); // CA is counted once
        assert!(r.missing.is_empty());
        assert_eq!(r.unexpected, vec!["XX", "YY"]);
    }

    #[test]
    fn test_empty_tree() {
        let db = CcdTemplateDb::global();
        let opts = MissingAtomsOptions::default();

        let empty = AtomGroup::new();
        let rep_empty = empty.check_missing_atoms(db, &opts);
        assert!(rep_empty.residues.is_empty());
        assert!(rep_empty.is_complete());
        assert_eq!(rep_empty.total_missing(), 0);
        assert!(rep_empty.incomplete().is_empty());
        assert!(rep_empty.no_template().is_empty());
        assert!(rep_empty.no_template_components().is_empty());
    }

    #[test]
    fn test_tree_unchanged_and_callable_across_hierarchy_levels() {
        let db = CcdTemplateDb::global();
        let opts = MissingAtomsOptions::default();

        // Hierarchy: root -> model_1 -> A -> residues "10", "2", "1" (checks sort_nicely and immutability)
        let mut root = AtomGroup::new();
        let mut model = AtomGroup::with_name("model_1");
        let mut chain = AtomGroup::with_name("A");
        chain.set_group(
            "10",
            build_residue_from_template("ALA", false, &["CB", "OXT"]),
        );
        chain.set_group("2", build_residue_from_template("HOH", false, &[]));
        chain.set_group("1", build_residue_from_template("ALA", false, &["OXT"]));
        model.set_group("A", chain);
        root.set_group("model_1", model);

        let root_before = root.clone();
        let rep_root = root.check_missing_atoms(db, &opts);
        assert_eq!(
            root, root_before,
            "tree must remain unchanged after check_missing_atoms"
        );

        // Verify sort_nicely order: /model_1/A/1/, /model_1/A/2/, /model_1/A/10/
        let paths: Vec<&str> = rep_root.residues.iter().map(|r| r.path.as_str()).collect();
        assert_eq!(
            paths,
            vec!["/model_1/A/1/", "/model_1/A/2/", "/model_1/A/10/"]
        );
        assert_eq!(rep_root.total_missing(), 1);
        assert_eq!(rep_root.incomplete().len(), 1);
        assert_eq!(rep_root.incomplete()[0].path, "/model_1/A/10/");

        // Callable at model, chain, and single-residue levels
        let model_ref = root.get_group("model_1").unwrap();
        assert_eq!(model_ref.check_missing_atoms(db, &opts), rep_root);

        let chain_ref = model_ref.get_group("A").unwrap();
        assert_eq!(chain_ref.check_missing_atoms(db, &opts), rep_root);

        let res1_ref = chain_ref.get_group("1").unwrap();
        let rep_res1 = res1_ref.check_missing_atoms(db, &opts);
        assert_eq!(rep_res1.residues.len(), 1);
        assert_eq!(rep_res1.residues[0].path, "/model_1/A/1/");
        assert_eq!(rep_res1.residues[0].status, ResidueStatus::Complete);
    }

    #[test]
    fn test_real_2fb4_and_2mgo() {
        let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let db = CcdTemplateDb::global();
        let default_opts = MissingAtomsOptions::default();
        let h_opts = MissingAtomsOptions {
            include_hydrogens: true,
        };

        // 2FB4 (X-ray structure)
        let cif_2fb4 = SimpleMmcif::from_file(data_dir.join("2FB4.cif")).expect("load 2FB4.cif");
        let ag_2fb4 = cif_2fb4
            .get_structure_atomgroup(None, None)
            .expect("2FB4 atomgroup");
        let ag_2fb4_before = ag_2fb4.clone();
        let rep_2fb4 = ag_2fb4.check_missing_atoms(db, &default_opts);
        assert_eq!(ag_2fb4, ag_2fb4_before);
        println!(
            "2FB4 (include_hydrogens=false): residues={}, incomplete={}, total_missing={}, no_template={:?}",
            rep_2fb4.residues.len(),
            rep_2fb4.incomplete().len(),
            rep_2fb4.total_missing(),
            rep_2fb4.no_template_components()
        );
        assert_eq!(rep_2fb4.residues.len(), 550);
        assert!(rep_2fb4.is_complete());
        assert_eq!(rep_2fb4.incomplete().len(), 0);
        assert_eq!(rep_2fb4.total_missing(), 0);
        assert!(rep_2fb4.no_template_components().is_empty());

        let rep_2fb4_h = ag_2fb4.check_missing_atoms(db, &h_opts);
        println!(
            "2FB4 (include_hydrogens=true): residues={}, incomplete={}, total_missing={}, no_template={:?}",
            rep_2fb4_h.residues.len(),
            rep_2fb4_h.incomplete().len(),
            rep_2fb4_h.total_missing(),
            rep_2fb4_h.no_template_components()
        );
        assert_eq!(rep_2fb4_h.residues.len(), 550);
        assert!(!rep_2fb4_h.is_complete());
        assert_eq!(rep_2fb4_h.incomplete().len(), 550);
        assert_eq!(rep_2fb4_h.total_missing(), 3031);
        assert!(rep_2fb4_h.no_template_components().is_empty());

        // 2MGO (NMR structure with 20 models and hydrogens)
        let cif_2mgo = SimpleMmcif::from_file(data_dir.join("2MGO.cif")).expect("load 2MGO.cif");
        let ag_2mgo = cif_2mgo
            .get_structure_atomgroup(None, None)
            .expect("2MGO atomgroup");
        let ag_2mgo_before = ag_2mgo.clone();
        let rep_2mgo = ag_2mgo.check_missing_atoms(db, &default_opts);
        assert_eq!(ag_2mgo, ag_2mgo_before);
        println!(
            "2MGO (include_hydrogens=false): residues={}, incomplete={}, total_missing={}, no_template={:?}",
            rep_2mgo.residues.len(),
            rep_2mgo.incomplete().len(),
            rep_2mgo.total_missing(),
            rep_2mgo.no_template_components()
        );
        assert_eq!(rep_2mgo.residues.len(), 180);
        assert!(rep_2mgo.is_complete());
        assert_eq!(rep_2mgo.incomplete().len(), 0);
        assert_eq!(rep_2mgo.total_missing(), 0);
        assert!(rep_2mgo.no_template_components().is_empty());

        let rep_2mgo_h = ag_2mgo.check_missing_atoms(db, &h_opts);
        println!(
            "2MGO (include_hydrogens=true): residues={}, incomplete={}, total_missing={}, no_template={:?}",
            rep_2mgo_h.residues.len(),
            rep_2mgo_h.incomplete().len(),
            rep_2mgo_h.total_missing(),
            rep_2mgo_h.no_template_components()
        );
        assert_eq!(rep_2mgo_h.residues.len(), 180);
        assert!(rep_2mgo_h.is_complete());
        assert_eq!(rep_2mgo_h.incomplete().len(), 0);
        assert_eq!(rep_2mgo_h.total_missing(), 0);
        assert!(rep_2mgo_h.no_template_components().is_empty());
    }
}
