// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Missing heavy-atom completion for macromolecular structures using CCD templates.
//!
//! Reconstructs missing expected heavy atoms in each residue/component by superposing
//! shared anchor heavy atoms from the Chemical Component Dictionary ([`CcdTemplateDb`])
//! onto the existing residue geometry via [`Superposer`] and transforming idealized
//! coordinates (`ideal_xyz`).

use std::collections::{HashMap, HashSet};

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::brd::sort_nicely;
use crate::ccd_templates::{CcdAtom, CcdTemplateDb};
use crate::hydrogenation::MIN_SUPERPOSE_HEAVY_ATOMS;
use crate::missing_atoms::{
    collect_from_group, expected_atoms, find_disulfide_bonded_ids, MissingAtomsOptions,
    TraversalState,
};
use crate::position::Position;
use crate::superposer::Superposer;

/// Default maximum allowed anchor superposition RMSD in Angstroms.
pub const DEFAULT_MAX_COMPLETE_RMSD: f64 = 1.0;

/// Options for missing heavy-atom completion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CompleteOptions {
    /// Maximum allowed heavy-atom anchor superposition RMSD in Angstroms (default: `1.0`).
    pub max_rmsd: f64,
}

impl Default for CompleteOptions {
    fn default() -> Self {
        Self {
            max_rmsd: DEFAULT_MAX_COMPLETE_RMSD,
        }
    }
}

/// Completion outcome status of a single residue/component.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CompletionStatus {
    /// Missing expected heavy atoms were reconstructed and added to the residue.
    Completed,
    /// All expected heavy atoms were already present; the residue was not modified.
    #[default]
    Unchanged,
    /// No CCD template was found in the database for this component; the residue was not modified.
    NoTemplate,
    /// Fewer than 3 anchor heavy atoms or collinear/degenerate anchors; the residue was not modified.
    InsufficientAnchors,
    /// Anchor superposition RMSD exceeded `max_rmsd` or a missing atom lacked `ideal_xyz`;
    /// the residue was not modified.
    PoorFit,
}

/// Per-residue result of missing heavy-atom completion.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ResidueCompletion {
    /// Hierarchical path of the residue group (e.g. `"/model_1/A/1/"`).
    pub path: String,
    /// Chemical component identifier (residue group's `name`, e.g. `"ALA"`).
    pub comp_id: String,
    /// Completion outcome (`Completed`, `Unchanged`, `NoTemplate`, `InsufficientAnchors`, or `PoorFit`).
    pub status: CompletionStatus,
    /// Names of heavy atoms added to the residue, in template order.
    pub added: Vec<String>,
    /// Number of anchor heavy atoms used for superposition.
    pub anchors: usize,
    /// Anchor superposition RMSD in Angstroms, if superposition was computed.
    pub rmsd: Option<f64>,
}

/// Structure-wide report of missing heavy-atom completion across residues.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CompletionReport {
    /// Completion results for all visited residues in traversal order.
    pub residues: Vec<ResidueCompletion>,
}

impl CompletionReport {
    /// Returns all residues with status [`CompletionStatus::Completed`].
    pub fn completed(&self) -> Vec<&ResidueCompletion> {
        self.residues
            .iter()
            .filter(|r| r.status == CompletionStatus::Completed)
            .collect()
    }

    /// Returns all skipped residues ([`CompletionStatus::NoTemplate`],
    /// [`CompletionStatus::InsufficientAnchors`], or [`CompletionStatus::PoorFit`]).
    pub fn skipped(&self) -> Vec<&ResidueCompletion> {
        self.residues
            .iter()
            .filter(|r| {
                matches!(
                    r.status,
                    CompletionStatus::NoTemplate
                        | CompletionStatus::InsufficientAnchors
                        | CompletionStatus::PoorFit
                )
            })
            .collect()
    }

    /// Returns the total number of heavy atoms added across all residues.
    pub fn total_added(&self) -> usize {
        self.residues.iter().map(|r| r.added.len()).sum()
    }
}

/// Plans heavy-atom completion for a single leaf residue without mutating it.
fn plan_single_residue(
    residue: &AtomGroup,
    prev_residue: Option<&AtomGroup>,
    is_disulfide_bonded: bool,
    db: &CcdTemplateDb,
    options: &CompleteOptions,
) -> (ResidueCompletion, Vec<(String, Atom)>) {
    let path = residue.path().to_string();
    let comp_id = residue.name.clone();

    let Some(template) = db.lookup(&comp_id) else {
        return (
            ResidueCompletion {
                path,
                comp_id,
                status: CompletionStatus::NoTemplate,
                added: Vec::new(),
                anchors: 0,
                rmsd: None,
            },
            Vec::new(),
        );
    };

    let missing_opts = MissingAtomsOptions {
        include_hydrogens: false,
    };
    let expected_list = expected_atoms(
        residue,
        prev_residue,
        is_disulfide_bonded,
        template,
        &missing_opts,
    );

    // Collect heavy atoms present in the residue in tree order (keeping the first appearance of each name)
    let mut present_heavy: HashMap<&str, &Atom> = HashMap::new();
    for (_, atom) in residue.atoms() {
        if atom.atomic_number() == 1 {
            continue;
        }
        present_heavy.entry(atom.name.as_str()).or_insert(atom);
    }

    let missing_atoms: Vec<&CcdAtom> = expected_list
        .into_iter()
        .filter(|exp_atom| !present_heavy.contains_key(exp_atom.name.as_str()))
        .collect();

    if missing_atoms.is_empty() {
        return (
            ResidueCompletion {
                path,
                comp_id,
                status: CompletionStatus::Unchanged,
                added: Vec::new(),
                anchors: 0,
                rmsd: None,
            },
            Vec::new(),
        );
    }

    // Build anchor groups: heavy atoms present in both residue and template (in template order)
    let mut template_heavy_group = AtomGroup::new();
    let mut actual_heavy_group = AtomGroup::new();

    for ccd_atom in &template.atoms {
        if ccd_atom.is_hydrogen() {
            continue;
        }
        let Some(&actual_atom) = present_heavy.get(ccd_atom.name.as_str()) else {
            continue;
        };
        let Some((ix, iy, iz)) = ccd_atom.ideal_xyz else {
            continue;
        };
        let Ok(mut ref_atom) = Atom::new_with_pos(&ccd_atom.element, Position::new(ix, iy, iz))
        else {
            continue;
        };
        ref_atom.name = ccd_atom.name.clone();
        template_heavy_group.set_atom(&ccd_atom.name, ref_atom);
        actual_heavy_group.set_atom(&ccd_atom.name, actual_atom.clone());
    }

    let anchors = template_heavy_group.get_number_of_atoms();
    if anchors < MIN_SUPERPOSE_HEAVY_ATOMS {
        return (
            ResidueCompletion {
                path,
                comp_id,
                status: CompletionStatus::InsufficientAnchors,
                added: Vec::new(),
                anchors,
                rmsd: None,
            },
            Vec::new(),
        );
    }

    let Ok(superposer) = Superposer::new(&template_heavy_group, &actual_heavy_group) else {
        return (
            ResidueCompletion {
                path,
                comp_id,
                status: CompletionStatus::InsufficientAnchors,
                added: Vec::new(),
                anchors,
                rmsd: None,
            },
            Vec::new(),
        );
    };

    let rmsd = superposer.rmsd();
    if !rmsd.is_finite() || rmsd > options.max_rmsd {
        return (
            ResidueCompletion {
                path,
                comp_id,
                status: CompletionStatus::PoorFit,
                added: Vec::new(),
                anchors,
                rmsd: Some(rmsd),
            },
            Vec::new(),
        );
    }

    let mut staged_atoms = Vec::with_capacity(missing_atoms.len());
    for ccd_atom in missing_atoms {
        let Some((ix, iy, iz)) = ccd_atom.ideal_xyz else {
            return (
                ResidueCompletion {
                    path,
                    comp_id,
                    status: CompletionStatus::PoorFit,
                    added: Vec::new(),
                    anchors,
                    rmsd: Some(rmsd),
                },
                Vec::new(),
            );
        };
        let Ok(pos) = superposer.transform_position(&Position::new(ix, iy, iz)) else {
            return (
                ResidueCompletion {
                    path,
                    comp_id,
                    status: CompletionStatus::PoorFit,
                    added: Vec::new(),
                    anchors,
                    rmsd: Some(rmsd),
                },
                Vec::new(),
            );
        };
        let Ok(mut new_atom) = Atom::new_with_pos(&ccd_atom.element, pos) else {
            return (
                ResidueCompletion {
                    path,
                    comp_id,
                    status: CompletionStatus::PoorFit,
                    added: Vec::new(),
                    anchors,
                    rmsd: Some(rmsd),
                },
                Vec::new(),
            );
        };
        new_atom.name = ccd_atom.name.clone();
        staged_atoms.push((ccd_atom.name.clone(), new_atom));
    }

    let added = staged_atoms.iter().map(|(name, _)| name.clone()).collect();
    (
        ResidueCompletion {
            path,
            comp_id,
            status: CompletionStatus::Completed,
            added,
            anchors,
            rmsd: Some(rmsd),
        },
        staged_atoms,
    )
}

/// Visits mutable leaf residue groups in the exact same traversal order as [`collect_from_group`].
fn visit_residues_mut(group: &mut AtomGroup, f: &mut impl FnMut(&mut AtomGroup)) {
    if group.get_number_of_groups() == 0 {
        if group.get_number_of_atoms() > 0 {
            f(group);
        }
        return;
    }

    let mut child_keys = group.get_group_list();
    sort_nicely(&mut child_keys);

    for key in &child_keys {
        if let Some(child) = group.get_group_mut(key) {
            visit_residues_mut(child, f);
        }
    }
}

/// Reconstructs and adds missing expected heavy atoms across an [`AtomGroup`] hierarchy
/// using CCD templates in `db`.
pub fn complete_missing_atoms(
    group: &mut AtomGroup,
    db: &CcdTemplateDb,
    options: &CompleteOptions,
) -> CompletionReport {
    let mut state = TraversalState::default();
    collect_from_group(group, None, true, &mut state);

    let mut disulfide_bonded = HashSet::new();
    for sgs in &state.sgs_by_scope {
        disulfide_bonded.extend(find_disulfide_bonded_ids(sgs));
    }

    let mut residues = Vec::with_capacity(state.residues.len());
    let mut staged_per_residue: Vec<Vec<(String, Atom)>> = Vec::with_capacity(state.residues.len());
    let mut any_added = false;

    for (idx, item) in state.residues.iter().enumerate() {
        let (res_comp, staged) = plan_single_residue(
            item.residue,
            item.prev_residue,
            disulfide_bonded.contains(&idx),
            db,
            options,
        );
        if !staged.is_empty() {
            any_added = true;
        }
        residues.push(res_comp);
        staged_per_residue.push(staged);
    }

    if any_added {
        let mut iter = staged_per_residue.into_iter();
        visit_residues_mut(group, &mut |residue| {
            if let Some(staged) = iter.next() {
                for (name, atom) in staged {
                    residue.set_atom(&name, atom);
                }
            }
        });
    }

    CompletionReport { residues }
}

impl AtomGroup {
    /// Reconstructs and adds missing expected heavy atoms in this hierarchy against
    /// CCD templates in `db`.
    pub fn complete_missing_atoms(
        &mut self,
        db: &CcdTemplateDb,
        options: &CompleteOptions,
    ) -> CompletionReport {
        complete_missing_atoms(self, db, options)
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;
    use crate::ccd_templates::CcdBondTemplate;
    use crate::format::mmcif::SimpleMmcif;
    use crate::missing_atoms::ResidueStatus;

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
    fn test_ala_missing_cb_completion() {
        let db = CcdTemplateDb::global();
        let opts = CompleteOptions::default();
        let check_opts = MissingAtomsOptions::default();

        // Build ALA without CB and optional OXT, then apply a rigid-body shift to verify non-trivial frame transform
        let mut ala = build_residue_from_template("ALA", false, &["CB", "OXT"]);
        ala.shift_by(Position::new(3.5, -2.0, 1.25));

        let orig_atoms: Vec<(String, Position)> =
            ala.atoms().map(|(k, a)| (k.clone(), a.xyz)).collect();

        let rep = ala.complete_missing_atoms(db, &opts);
        assert_eq!(rep.residues.len(), 1);
        assert_eq!(rep.completed().len(), 1);
        assert!(rep.skipped().is_empty());
        assert_eq!(rep.total_added(), 1);

        let r = &rep.residues[0];
        assert_eq!(r.comp_id, "ALA");
        assert_eq!(r.status, CompletionStatus::Completed);
        assert_eq!(r.added, vec!["CB"]);
        assert_eq!(r.anchors, 4);
        assert!(r.rmsd.expect("rmsd present") < 1e-6);

        // Check that check_missing_atoms reports Complete and optional OXT was not added
        let check_rep = ala.check_missing_atoms(db, &check_opts);
        assert!(check_rep.is_complete());
        assert_eq!(check_rep.residues[0].status, ResidueStatus::Complete);
        assert!(!ala.has_atom("OXT"));

        // CB-CA bond distance ~ 1.52 Å
        let ca = ala.get_atom("CA").expect("CA present");
        let cb = ala.get_atom("CB").expect("CB added");
        let dist = cb.xyz.distance_from(&ca.xyz);
        assert!(
            (dist - 1.524).abs() < 0.02,
            "expected CB-CA distance ~1.52 Å, got {dist:.4}"
        );

        // Original atoms' coordinates must be completely unchanged
        for (key, orig_xyz) in &orig_atoms {
            let current = ala.get_atom(key).expect("original atom still present");
            assert_eq!(current.xyz, *orig_xyz);
        }

        // Duplicate anchor name: first atom in tree order is used as the anchor
        let mut ala_dup = build_residue_from_template("ALA", false, &["CB", "OXT"]);
        ala_dup.set_atom(
            "99_N_dup",
            make_atom("N", "N", Position::new(50.0, 50.0, 50.0)),
        );
        let rep_dup = ala_dup.complete_missing_atoms(db, &opts);
        assert_eq!(rep_dup.residues[0].status, CompletionStatus::Completed);
        assert_eq!(rep_dup.residues[0].added, vec!["CB"]);
        assert_eq!(rep_dup.residues[0].anchors, 4);
        assert!(rep_dup.residues[0].rmsd.unwrap() < 1e-6);
    }

    #[test]
    fn test_lys_missing_cd_ce_nz() {
        let db = CcdTemplateDb::global();
        let opts = CompleteOptions::default();
        let check_opts = MissingAtomsOptions::default();

        let full_lys = build_residue_from_template("LYS", false, &["OXT"]);
        let mut lys = build_residue_from_template("LYS", false, &["CD", "CE", "NZ", "OXT"]);

        let rep = lys.complete_missing_atoms(db, &opts);
        assert_eq!(rep.residues.len(), 1);
        assert_eq!(rep.completed().len(), 1);
        assert!(rep.skipped().is_empty());
        assert_eq!(rep.total_added(), 3);

        let r = &rep.residues[0];
        assert_eq!(r.comp_id, "LYS");
        assert_eq!(r.status, CompletionStatus::Completed);
        assert_eq!(r.added, vec!["CD", "CE", "NZ"]);
        assert_eq!(r.anchors, 6); // N, CA, C, O, CB, CG
        assert!(r.rmsd.expect("rmsd present") < 1e-6);

        let check_rep = lys.check_missing_atoms(db, &check_opts);
        assert!(check_rep.is_complete());

        for name in ["CD", "CE", "NZ"] {
            let added_atom = lys.get_atom(name).expect("added atom present");
            let expected_atom = full_lys.get_atom(name).expect("template atom present");
            assert!(added_atom.xyz.distance_from(&expected_atom.xyz) < 1e-5);
        }
    }

    #[test]
    fn test_insufficient_anchors_fewer_than_three_and_collinear() {
        let db = CcdTemplateDb::global();
        let opts = CompleteOptions::default();

        // 1. ALA with only N and CA (2 anchors < 3)
        let mut ala_two = build_residue_from_template("ALA", false, &["C", "O", "CB", "OXT"]);
        let ala_two_before = ala_two.clone();

        let rep = ala_two.complete_missing_atoms(db, &opts);
        assert_eq!(rep.residues.len(), 1);
        assert!(rep.completed().is_empty());
        assert_eq!(rep.skipped().len(), 1);
        assert_eq!(rep.total_added(), 0);
        assert_eq!(ala_two, ala_two_before);

        let r = &rep.residues[0];
        assert_eq!(r.status, CompletionStatus::InsufficientAnchors);
        assert!(r.added.is_empty());
        assert_eq!(r.anchors, 2);
        assert_eq!(r.rmsd, None);

        // 2. ALA with 3 collinear anchors (N, CA, C placed on a line, missing O and CB)
        let mut ala_collinear = AtomGroup::with_name("ALA");
        ala_collinear.set_atom("N", make_atom("N", "N", Position::new(0.0, 0.0, 0.0)));
        ala_collinear.set_atom("CA", make_atom("C", "CA", Position::new(1.5, 0.0, 0.0)));
        ala_collinear.set_atom("C", make_atom("C", "C", Position::new(3.0, 0.0, 0.0)));
        let ala_collinear_before = ala_collinear.clone();

        let rep_col = ala_collinear.complete_missing_atoms(db, &opts);
        assert_eq!(rep_col.residues.len(), 1);
        assert!(rep_col.completed().is_empty());
        assert_eq!(rep_col.skipped().len(), 1);
        assert_eq!(ala_collinear, ala_collinear_before);

        let r_col = &rep_col.residues[0];
        assert_eq!(r_col.status, CompletionStatus::InsufficientAnchors);
        assert!(r_col.added.is_empty());
        assert_eq!(r_col.anchors, 3);
        assert_eq!(r_col.rmsd, None);
    }

    #[test]
    fn test_poor_fit_swapped_names_and_relaxed_max_rmsd() {
        let db = CcdTemplateDb::global();

        // Build ALA missing CB and OXT, then swap the coordinates of N and O to corrupt anchor naming
        let mut ala = build_residue_from_template("ALA", false, &["CB", "OXT"]);
        let n_pos = ala.get_atom("N").unwrap().xyz;
        let o_pos = ala.get_atom("O").unwrap().xyz;
        ala.get_atom_mut("1_N").unwrap().xyz = o_pos;
        ala.get_atom_mut("4_O").unwrap().xyz = n_pos;
        let ala_before = ala.clone();

        // With default max_rmsd (1.0 Å), should fail with PoorFit and leave tree unmodified
        let rep_strict = ala.complete_missing_atoms(db, &CompleteOptions::default());
        assert_eq!(rep_strict.residues.len(), 1);
        assert!(rep_strict.completed().is_empty());
        assert_eq!(rep_strict.skipped().len(), 1);
        assert_eq!(rep_strict.total_added(), 0);
        assert_eq!(ala, ala_before);

        let r_strict = &rep_strict.residues[0];
        assert_eq!(r_strict.status, CompletionStatus::PoorFit);
        assert!(r_strict.added.is_empty());
        assert_eq!(r_strict.anchors, 4);
        let fit_rmsd = r_strict.rmsd.expect("rmsd should be computed");
        assert!(
            fit_rmsd > 1.0,
            "corrupted naming RMSD should exceed 1.0 Å, got {fit_rmsd:.4}"
        );

        // With relaxed max_rmsd (5.0 Å), should complete
        let relaxed_opts = CompleteOptions { max_rmsd: 5.0 };
        let rep_relaxed = ala.complete_missing_atoms(db, &relaxed_opts);
        assert_eq!(rep_relaxed.completed().len(), 1);
        assert!(rep_relaxed.skipped().is_empty());
        assert_eq!(rep_relaxed.total_added(), 1);

        let r_relaxed = &rep_relaxed.residues[0];
        assert_eq!(r_relaxed.status, CompletionStatus::Completed);
        assert_eq!(r_relaxed.added, vec!["CB"]);
        assert_eq!(r_relaxed.anchors, 4);
        assert_eq!(r_relaxed.rmsd, Some(fit_rmsd));
        assert!(ala.has_atom("CB"));
    }

    #[test]
    fn test_da_missing_n9_and_optional_phosphates_not_added() {
        let db = CcdTemplateDb::global();
        let opts = CompleteOptions::default();
        let check_opts = MissingAtomsOptions::default();

        // 5'-terminal DA lacking optional 5'-phosphate (P, OP1, OP2, OP3) and missing required N9
        let full_da = build_residue_from_template("DA", false, &["P", "OP1", "OP2", "OP3"]);
        let mut da = build_residue_from_template("DA", false, &["P", "OP1", "OP2", "OP3", "N9"]);

        let rep = da.complete_missing_atoms(db, &opts);
        assert_eq!(rep.residues.len(), 1);
        assert_eq!(rep.completed().len(), 1);
        assert_eq!(rep.total_added(), 1);

        let r = &rep.residues[0];
        assert_eq!(r.comp_id, "DA");
        assert_eq!(r.status, CompletionStatus::Completed);
        assert_eq!(r.added, vec!["N9"]);
        assert!(r.rmsd.expect("rmsd present") < 1e-6);

        // Optional phosphate atoms must not be added
        for opt_name in ["P", "OP1", "OP2", "OP3"] {
            assert!(!da.has_atom(opt_name));
        }

        let check_rep = da.check_missing_atoms(db, &check_opts);
        assert!(check_rep.is_complete());

        let n9 = da.get_atom("N9").expect("N9 added");
        let expected_n9 = full_da.get_atom("N9").unwrap();
        assert!(n9.xyz.distance_from(&expected_n9.xyz) < 1e-5);
    }

    #[test]
    fn test_unknown_component_no_template_and_custom_db_completion() {
        let db = CcdTemplateDb::global();
        let opts = CompleteOptions::default();
        let check_opts = MissingAtomsOptions::default();

        let mut lig = AtomGroup::with_name("LIG");
        lig.set_atom("C1", make_atom("C", "C1", Position::new(0.0, 0.0, 0.0)));
        lig.set_atom("C2", make_atom("C", "C2", Position::new(1.5, 0.0, 0.0)));
        lig.set_atom("O1", make_atom("O", "O1", Position::new(1.5, 1.4, 0.0)));
        lig.set_atom("C3", make_atom("C", "C3", Position::new(0.0, 0.8, 1.2)));
        let lig_before = lig.clone();

        // 1. Without template in global DB -> NoTemplate
        let rep_no_tpl = lig.complete_missing_atoms(db, &opts);
        assert_eq!(rep_no_tpl.residues.len(), 1);
        assert!(rep_no_tpl.completed().is_empty());
        assert_eq!(rep_no_tpl.skipped().len(), 1);
        assert_eq!(rep_no_tpl.total_added(), 0);
        assert_eq!(lig, lig_before);

        let r_no = &rep_no_tpl.residues[0];
        assert_eq!(r_no.comp_id, "LIG");
        assert_eq!(r_no.status, CompletionStatus::NoTemplate);
        assert!(r_no.added.is_empty());
        assert_eq!(r_no.anchors, 0);
        assert_eq!(r_no.rmsd, None);

        // 2. Insert custom LIG template with an additional heavy atom N1
        let mut custom_db = db.clone();
        custom_db.insert(CcdBondTemplate {
            comp_id: "LIG".to_string(),
            atoms: vec![
                CcdAtom {
                    name: "C1".to_string(),
                    element: "C".to_string(),
                    ideal_xyz: Some((0.0, 0.0, 0.0)),
                },
                CcdAtom {
                    name: "C2".to_string(),
                    element: "C".to_string(),
                    ideal_xyz: Some((1.5, 0.0, 0.0)),
                },
                CcdAtom {
                    name: "O1".to_string(),
                    element: "O".to_string(),
                    ideal_xyz: Some((1.5, 1.4, 0.0)),
                },
                CcdAtom {
                    name: "C3".to_string(),
                    element: "C".to_string(),
                    ideal_xyz: Some((0.0, 0.8, 1.2)),
                },
                CcdAtom {
                    name: "N1".to_string(),
                    element: "N".to_string(),
                    ideal_xyz: Some((0.0, -1.4, 0.5)),
                },
            ],
            bonds: Vec::new(),
        });

        let rep_custom = lig.complete_missing_atoms(&custom_db, &opts);
        assert_eq!(rep_custom.completed().len(), 1);
        assert!(rep_custom.skipped().is_empty());
        assert_eq!(rep_custom.total_added(), 1);

        let r_ok = &rep_custom.residues[0];
        assert_eq!(r_ok.status, CompletionStatus::Completed);
        assert_eq!(r_ok.added, vec!["N1"]);
        assert_eq!(r_ok.anchors, 4);
        assert!(r_ok.rmsd.unwrap() < 1e-6);
        assert!(lig
            .check_missing_atoms(&custom_db, &check_opts)
            .is_complete());
        let n1 = lig.get_atom("N1").expect("N1 added");
        assert!(n1.xyz.distance_from(&Position::new(0.0, -1.4, 0.5)) < 1e-5);

        // 3. Missing atom without ideal_xyz -> PoorFit and residue untouched
        let mut no_xyz_db = db.clone();
        no_xyz_db.insert(CcdBondTemplate {
            comp_id: "LIG".to_string(),
            atoms: vec![
                CcdAtom {
                    name: "C1".to_string(),
                    element: "C".to_string(),
                    ideal_xyz: Some((0.0, 0.0, 0.0)),
                },
                CcdAtom {
                    name: "C2".to_string(),
                    element: "C".to_string(),
                    ideal_xyz: Some((1.5, 0.0, 0.0)),
                },
                CcdAtom {
                    name: "O1".to_string(),
                    element: "O".to_string(),
                    ideal_xyz: Some((1.5, 1.4, 0.0)),
                },
                CcdAtom {
                    name: "C3".to_string(),
                    element: "C".to_string(),
                    ideal_xyz: Some((0.0, 0.8, 1.2)),
                },
                CcdAtom {
                    name: "N1".to_string(),
                    element: "N".to_string(),
                    ideal_xyz: None,
                },
            ],
            bonds: Vec::new(),
        });
        let mut lig2 = lig_before.clone();
        let rep_no_xyz = lig2.complete_missing_atoms(&no_xyz_db, &opts);
        assert_eq!(rep_no_xyz.residues[0].status, CompletionStatus::PoorFit);
        assert!(rep_no_xyz.residues[0].added.is_empty());
        assert_eq!(lig2, lig_before);
    }

    #[test]
    fn test_complete_structure_unchanged_and_idempotent_second_call() {
        let db = CcdTemplateDb::global();
        let opts = CompleteOptions::default();

        // 1. Already complete structure is unchanged
        let mut complete_ala = build_residue_from_template("ALA", false, &["OXT"]);
        let before = complete_ala.clone();
        let rep1 = complete_ala.complete_missing_atoms(db, &opts);
        assert_eq!(rep1.residues.len(), 1);
        assert_eq!(rep1.residues[0].status, CompletionStatus::Unchanged);
        assert!(rep1.residues[0].added.is_empty());
        assert_eq!(rep1.residues[0].anchors, 0);
        assert_eq!(rep1.residues[0].rmsd, None);
        assert!(rep1.completed().is_empty());
        assert!(rep1.skipped().is_empty());
        assert_eq!(rep1.total_added(), 0);
        assert_eq!(complete_ala, before);

        // 2. Completing an incomplete residue and calling a second time adds nothing on the second call
        let mut inc_ala = build_residue_from_template("ALA", false, &["CB", "OXT"]);
        let rep_first = inc_ala.complete_missing_atoms(db, &opts);
        assert_eq!(rep_first.total_added(), 1);
        assert_eq!(rep_first.residues[0].status, CompletionStatus::Completed);

        let after_first = inc_ala.clone();
        let rep_second = inc_ala.complete_missing_atoms(db, &opts);
        assert_eq!(rep_second.total_added(), 0);
        assert_eq!(rep_second.residues[0].status, CompletionStatus::Unchanged);
        assert_eq!(inc_ala, after_first);
    }

    #[test]
    fn test_empty_tree() {
        let db = CcdTemplateDb::global();
        let opts = CompleteOptions::default();

        let mut empty = AtomGroup::new();
        let rep = empty.complete_missing_atoms(db, &opts);
        assert!(rep.residues.is_empty());
        assert!(rep.completed().is_empty());
        assert!(rep.skipped().is_empty());
        assert_eq!(rep.total_added(), 0);
    }

    #[test]
    fn test_callable_across_hierarchy_levels() {
        let db = CcdTemplateDb::global();
        let opts = CompleteOptions::default();

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

        let mut root_copy = root.clone();
        let mut model_copy = root.get_group("model_1").unwrap().clone();
        let mut chain_copy = root
            .get_group("model_1")
            .unwrap()
            .get_group("A")
            .unwrap()
            .clone();
        let mut res10_copy = chain_copy.get_group("10").unwrap().clone();

        let rep_root = root_copy.complete_missing_atoms(db, &opts);
        let rep_model = model_copy.complete_missing_atoms(db, &opts);
        let rep_chain = chain_copy.complete_missing_atoms(db, &opts);
        let rep_res10 = res10_copy.complete_missing_atoms(db, &opts);

        let paths: Vec<&str> = rep_root.residues.iter().map(|r| r.path.as_str()).collect();
        assert_eq!(
            paths,
            vec!["/model_1/A/1/", "/model_1/A/2/", "/model_1/A/10/"]
        );
        assert_eq!(rep_root.total_added(), 1);
        assert_eq!(rep_root.completed().len(), 1);
        assert_eq!(rep_root.completed()[0].path, "/model_1/A/10/");

        assert_eq!(rep_model, rep_root);
        assert_eq!(rep_chain, rep_root);
        assert_eq!(rep_res10.residues.len(), 1);
        assert_eq!(rep_res10.residues[0], rep_root.residues[2]);
        assert_eq!(res10_copy, *chain_copy.get_group("10").unwrap());
    }

    #[test]
    fn test_real_2fb4_and_2mgo() {
        let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let db = CcdTemplateDb::global();
        let opts = CompleteOptions::default();
        let check_opts = MissingAtomsOptions::default();

        // 1. 2FB4 (X-ray structure, 550 residues)
        // Remove 7 heavy atoms across 5 residues and reconstruct them
        let cif_2fb4 = SimpleMmcif::from_file(data_dir.join("2FB4.cif")).expect("load 2FB4.cif");
        let mut ag_2fb4 = cif_2fb4
            .get_structure_atomgroup(None, None)
            .expect("2FB4 atomgroup");

        let targets_2fb4: &[(&str, &str, &[&str])] = &[
            ("H", "14", &["CB", "CG"]), // PRO
            ("H", "33", &["CB"]),       // ALA
            ("H", "40", &["CB"]),       // ALA
            ("H", "41", &["CB", "CG"]), // PRO
            ("L", "13", &["CB"]),       // THR
        ];
        let mut orig_2fb4 = Vec::new();
        for &(chain_id, res_id, atom_names) in targets_2fb4 {
            let res = ag_2fb4
                .get_group_mut("model_1")
                .unwrap()
                .get_group_mut(chain_id)
                .unwrap()
                .get_group_mut(res_id)
                .unwrap();
            let mut res_orig = Vec::new();
            for &aname in atom_names {
                let key = res
                    .atoms()
                    .find(|(_, a)| a.name == aname)
                    .map(|(k, _)| k.clone())
                    .expect("target atom exists");
                let removed = res.remove_atom(&key).unwrap();
                res_orig.push((aname, removed.xyz));
            }
            orig_2fb4.push((chain_id, res_id, res_orig));
        }

        let removed_count_2fb4: usize = orig_2fb4.iter().map(|(_, _, v)| v.len()).sum();
        assert_eq!(removed_count_2fb4, 7);
        assert_eq!(
            ag_2fb4.check_missing_atoms(db, &check_opts).total_missing(),
            7
        );

        let rep_2fb4 = ag_2fb4.complete_missing_atoms(db, &opts);
        assert_eq!(rep_2fb4.residues.len(), 550);
        assert_eq!(rep_2fb4.completed().len(), 5);
        assert!(rep_2fb4.skipped().is_empty());
        assert_eq!(rep_2fb4.total_added(), 7);
        assert!(ag_2fb4.check_missing_atoms(db, &check_opts).is_complete());

        let max_anchor_rmsd_2fb4 = rep_2fb4
            .completed()
            .iter()
            .filter_map(|r| r.rmsd)
            .fold(0.0_f64, f64::max);

        let mut max_res_reconstructed_rmsd_2fb4 = 0.0_f64;
        let mut max_atom_dist_2fb4 = 0.0_f64;
        let mut total_sq_err_2fb4 = 0.0_f64;
        for (chain_id, res_id, atoms) in &orig_2fb4 {
            let res = ag_2fb4
                .get_group("model_1")
                .unwrap()
                .get_group(chain_id)
                .unwrap()
                .get_group(res_id)
                .unwrap();
            let mut res_sq = 0.0_f64;
            for &(aname, orig_xyz) in atoms {
                let new_xyz = res.get_atom(aname).expect("reconstructed").xyz;
                let dist = orig_xyz.distance_from(&new_xyz);
                max_atom_dist_2fb4 = max_atom_dist_2fb4.max(dist);
                res_sq += dist * dist;
                total_sq_err_2fb4 += dist * dist;
            }
            let res_rmsd = (res_sq / (atoms.len() as f64)).sqrt();
            max_res_reconstructed_rmsd_2fb4 = max_res_reconstructed_rmsd_2fb4.max(res_rmsd);
        }
        let overall_reconstructed_rmsd_2fb4 =
            (total_sq_err_2fb4 / (removed_count_2fb4 as f64)).sqrt();
        println!(
            "2FB4 completion: removed={}, added={}, completed_residues={}, max_anchor_rmsd={:.4} Å, max_residue_reconstructed_rmsd={:.4} Å, max_atom_dist={:.4} Å, overall_reconstructed_rmsd={:.4} Å",
            removed_count_2fb4,
            rep_2fb4.total_added(),
            rep_2fb4.completed().len(),
            max_anchor_rmsd_2fb4,
            max_res_reconstructed_rmsd_2fb4,
            max_atom_dist_2fb4,
            overall_reconstructed_rmsd_2fb4
        );
        assert!(max_anchor_rmsd_2fb4 < 0.35);
        assert!(max_res_reconstructed_rmsd_2fb4 < 0.40);
        assert!(max_atom_dist_2fb4 < 0.45);

        // 2. 2MGO (NMR structure, 20 models, 180 residues)
        // Remove 7 heavy atoms across 5 residues and reconstruct them
        let cif_2mgo = SimpleMmcif::from_file(data_dir.join("2MGO.cif")).expect("load 2MGO.cif");
        let mut ag_2mgo = cif_2mgo
            .get_structure_atomgroup(None, None)
            .expect("2MGO atomgroup");
        let targets_2mgo: &[(&str, &str, &str, &[&str])] = &[
            ("model_1", "A", "1", &["CB"]),       // CYS
            ("model_1", "A", "3", &["CB"]),       // ILE
            ("model_1", "A", "7", &["CB", "CG"]), // PRO
            ("model_2", "A", "7", &["CB", "CG"]), // PRO
            ("model_3", "A", "7", &["CB"]),       // PRO
        ];
        let mut orig_2mgo = Vec::new();
        for &(model_id, chain_id, res_id, atom_names) in targets_2mgo {
            let res = ag_2mgo
                .get_group_mut(model_id)
                .unwrap()
                .get_group_mut(chain_id)
                .unwrap()
                .get_group_mut(res_id)
                .unwrap();
            let mut res_orig = Vec::new();
            for &aname in atom_names {
                let key = res
                    .atoms()
                    .find(|(_, a)| a.name == aname)
                    .map(|(k, _)| k.clone())
                    .expect("target atom exists");
                let removed = res.remove_atom(&key).unwrap();
                res_orig.push((aname, removed.xyz));
            }
            orig_2mgo.push((model_id, chain_id, res_id, res_orig));
        }

        let removed_count_2mgo: usize = orig_2mgo.iter().map(|(_, _, _, v)| v.len()).sum();
        assert_eq!(removed_count_2mgo, 7);
        assert_eq!(
            ag_2mgo.check_missing_atoms(db, &check_opts).total_missing(),
            7
        );

        let rep_2mgo = ag_2mgo.complete_missing_atoms(db, &opts);
        assert_eq!(rep_2mgo.residues.len(), 180);
        assert_eq!(rep_2mgo.completed().len(), 5);
        assert!(rep_2mgo.skipped().is_empty());
        assert_eq!(rep_2mgo.total_added(), 7);
        assert!(ag_2mgo.check_missing_atoms(db, &check_opts).is_complete());

        let max_anchor_rmsd_2mgo = rep_2mgo
            .completed()
            .iter()
            .filter_map(|r| r.rmsd)
            .fold(0.0_f64, f64::max);

        let mut max_res_reconstructed_rmsd_2mgo = 0.0_f64;
        let mut max_atom_dist_2mgo = 0.0_f64;
        let mut total_sq_err_2mgo = 0.0_f64;
        for (model_id, chain_id, res_id, atoms) in &orig_2mgo {
            let res = ag_2mgo
                .get_group(model_id)
                .unwrap()
                .get_group(chain_id)
                .unwrap()
                .get_group(res_id)
                .unwrap();
            let mut res_sq = 0.0_f64;
            for &(aname, orig_xyz) in atoms {
                let new_xyz = res.get_atom(aname).expect("reconstructed").xyz;
                let dist = orig_xyz.distance_from(&new_xyz);
                max_atom_dist_2mgo = max_atom_dist_2mgo.max(dist);
                res_sq += dist * dist;
                total_sq_err_2mgo += dist * dist;
            }
            let res_rmsd = (res_sq / (atoms.len() as f64)).sqrt();
            max_res_reconstructed_rmsd_2mgo = max_res_reconstructed_rmsd_2mgo.max(res_rmsd);
        }
        let overall_reconstructed_rmsd_2mgo =
            (total_sq_err_2mgo / (removed_count_2mgo as f64)).sqrt();
        println!(
            "2MGO completion: removed={}, added={}, completed_residues={}, max_anchor_rmsd={:.4} Å, max_residue_reconstructed_rmsd={:.4} Å, max_atom_dist={:.4} Å, overall_reconstructed_rmsd={:.4} Å",
            removed_count_2mgo,
            rep_2mgo.total_added(),
            rep_2mgo.completed().len(),
            max_anchor_rmsd_2mgo,
            max_res_reconstructed_rmsd_2mgo,
            max_atom_dist_2mgo,
            overall_reconstructed_rmsd_2mgo
        );
        assert!(max_anchor_rmsd_2mgo < 0.60);
        assert!(max_res_reconstructed_rmsd_2mgo < 0.55);
        assert!(max_atom_dist_2mgo < 0.70);
    }
}
