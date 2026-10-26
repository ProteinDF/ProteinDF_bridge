// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Solvent-accessible surface area (SASA) calculation using the Shrake-Rupley algorithm.
//!
//! Computes per-atom and per-residue solvent-accessible surface areas on heavy atoms
//! using a Fibonacci sphere point distribution and Bondi (1964) van der Waals radii.

use std::f64::consts::PI;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::brd::sort_nicely;
use crate::position::Position;

/// Default solvent probe radius in Angstroms (water molecule, 1.4 Å).
pub const DEFAULT_PROBE_RADIUS: f64 = 1.4;

/// Default number of test points on each expanded atomic sphere.
pub const DEFAULT_N_POINTS: usize = 92;

/// Default van der Waals radius in Angstroms for elements not in Bondi's table.
pub const DEFAULT_BONDI_RADIUS: f64 = 1.8;

/// Bondi (1964) van der Waals radii in Angstroms, keyed by element symbol.
///
/// Reference: A. Bondi, "van der Waals Volumes and Radii",
/// *J. Phys. Chem.* 1964, 68, 3, 441–451.
pub const BONDI_RADII: &[(&str, f64)] = &[
    ("H", 1.20),
    ("He", 1.40),
    ("Li", 1.82),
    ("C", 1.70),
    ("N", 1.55),
    ("O", 1.52),
    ("F", 1.47),
    ("Ne", 1.54),
    ("Na", 2.27),
    ("Mg", 1.73),
    ("Si", 2.10),
    ("P", 1.80),
    ("S", 1.80),
    ("Cl", 1.75),
    ("Ar", 1.88),
    ("K", 2.75),
    ("Ni", 1.63),
    ("Cu", 1.40),
    ("Zn", 1.39),
    ("Ga", 1.87),
    ("As", 1.85),
    ("Se", 1.90),
    ("Br", 1.85),
    ("Kr", 2.02),
    ("Pd", 1.63),
    ("Ag", 1.72),
    ("Cd", 1.58),
    ("In", 1.93),
    ("Sn", 2.17),
    ("Te", 2.06),
    ("I", 1.98),
    ("Xe", 2.16),
    ("Pt", 1.75),
    ("Au", 1.66),
    ("Hg", 1.55),
    ("Tl", 1.96),
    ("Pb", 2.02),
    ("U", 1.86),
];

/// Returns the Bondi van der Waals radius (in Å) for the given element symbol.
///
/// Elements not present in Bondi's table fall back to [`DEFAULT_BONDI_RADIUS`] (1.8 Å).
pub fn bondi_radius(symbol: &str) -> f64 {
    let trimmed = symbol.trim();
    for &(sym, radius) in BONDI_RADII {
        if sym.eq_ignore_ascii_case(trimmed) {
            return radius;
        }
    }
    DEFAULT_BONDI_RADIUS
}

/// Options for solvent-accessible surface area (SASA) calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct SasaOptions {
    /// Solvent probe radius in Angstroms (default: `1.4`).
    pub probe_radius: f64,
    /// Number of test points placed on each sphere via Fibonacci lattice (default: `92`).
    pub n_points: usize,
}

impl Default for SasaOptions {
    fn default() -> Self {
        Self {
            probe_radius: DEFAULT_PROBE_RADIUS,
            n_points: DEFAULT_N_POINTS,
        }
    }
}

/// Solvent-accessible surface area for a single heavy atom.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct AtomSasa {
    /// Hierarchical path of the atom (e.g. `"/model_1/A/1/CA"`).
    pub path: String,
    /// Atom name (e.g. `"CA"`).
    pub name: String,
    /// Solvent-accessible surface area in Ų.
    pub area: f64,
}

/// Solvent-accessible surface area aggregated for a single residue/component.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ResidueSasa {
    /// Hierarchical path of the residue group (e.g. `"/model_1/A/1/"`).
    pub path: String,
    /// Chemical component identifier (residue group's `name`, e.g. `"ALA"`).
    pub comp_id: String,
    /// Sum of heavy-atom solvent-accessible surface areas in this residue in Ų.
    pub area: f64,
}

/// Structure-wide report of per-atom and per-residue solvent-accessible surface areas.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SasaReport {
    /// Per-atom SASA results for all heavy atoms in traversal order.
    pub atoms: Vec<AtomSasa>,
    /// Per-residue SASA results for all leaf residue groups in traversal order.
    pub residues: Vec<ResidueSasa>,
}

/// Generates `n_points` quasi-uniformly distributed points on the unit sphere
/// using a Fibonacci spiral lattice.
fn fibonacci_sphere(n_points: usize) -> Vec<Position> {
    if n_points == 0 {
        return Vec::new();
    }
    let golden_angle = PI * (3.0 - 5.0_f64.sqrt());
    let n = n_points as f64;
    let mut points = Vec::with_capacity(n_points);
    for i in 0..n_points {
        let y = 1.0 - (2.0 * (i as f64 + 0.5)) / n;
        let radius = (1.0 - y * y).max(0.0).sqrt();
        let theta = golden_angle * (i as f64);
        let x = radius * theta.cos();
        let z = radius * theta.sin();
        points.push(Position::new(x, y, z));
    }
    points
}

/// Returns true if `atom` is hydrogen.
fn is_hydrogen(atom: &Atom) -> bool {
    atom.atomic_number() == 1 || atom.symbol().is_ok_and(|s| s.eq_ignore_ascii_case("H"))
}

/// Returns the Bondi van der Waals radius for `atom`, falling back to [`DEFAULT_BONDI_RADIUS`].
fn atom_bondi_radius(atom: &Atom) -> f64 {
    match atom.symbol() {
        Ok(sym) => bondi_radius(sym),
        Err(_) => DEFAULT_BONDI_RADIUS,
    }
}

/// Collects leaf residue groups (groups with no subgroups and at least one direct atom)
/// in the same traversal order as `orchestrator::traverse_and_hydrogenate`.
fn collect_residues<'a>(group: &'a AtomGroup, residues: &mut Vec<&'a AtomGroup>) {
    if group.get_number_of_groups() == 0 {
        if group.get_number_of_atoms() > 0 {
            residues.push(group);
        }
        return;
    }

    let mut child_keys = group.get_group_list();
    sort_nicely(&mut child_keys);

    for key in &child_keys {
        if let Some(child) = group.get_group(key) {
            collect_residues(child, residues);
        }
    }
}

struct HeavyAtomEntry {
    residue_idx: usize,
    path: String,
    name: String,
    xyz: Position,
    expanded_radius: f64,
}

impl AtomGroup {
    /// Calculates per-atom and per-residue solvent-accessible surface areas (SASA)
    /// for all heavy atoms within this group hierarchy using the Shrake-Rupley method.
    ///
    /// - Hydrogen atoms are ignored.
    /// - Occlusion is evaluated exclusively among heavy atoms within `self`.
    /// - The hierarchy is not modified.
    pub fn sasa(&self, options: &SasaOptions) -> SasaReport {
        let mut residues = Vec::new();
        collect_residues(self, &mut residues);

        let mut heavy_atoms = Vec::new();
        for (res_idx, residue) in residues.iter().enumerate() {
            for (key, atom) in residue.atoms() {
                if is_hydrogen(atom) {
                    continue;
                }
                let path = if atom.path.is_empty() {
                    format!("{}{}", residue.path(), key)
                } else {
                    atom.path.clone()
                };
                let vdw = atom_bondi_radius(atom);
                let expanded_radius = (vdw + options.probe_radius).max(0.0);
                heavy_atoms.push(HeavyAtomEntry {
                    residue_idx: res_idx,
                    path,
                    name: atom.name.clone(),
                    xyz: atom.xyz,
                    expanded_radius,
                });
            }
        }

        let unit_points = fibonacci_sphere(options.n_points);
        let n_heavy = heavy_atoms.len();
        let mut atom_reports = Vec::with_capacity(n_heavy);
        let mut residue_areas = vec![0.0; residues.len()];

        // Precompute squared expanded radii for fast inside-sphere tests
        let expanded_radii_sq: Vec<f64> = heavy_atoms
            .iter()
            .map(|a| a.expanded_radius * a.expanded_radius)
            .collect();

        let mut neighbors = Vec::new();

        for (i, atom_i) in heavy_atoms.iter().enumerate() {
            let r_i = atom_i.expanded_radius;
            let area = if options.n_points == 0 || r_i <= 0.0 {
                0.0
            } else {
                neighbors.clear();
                for (j, atom_j) in heavy_atoms.iter().enumerate() {
                    if i == j {
                        continue;
                    }
                    let r_j = atom_j.expanded_radius;
                    if r_j <= 0.0 {
                        continue;
                    }
                    let max_dist = r_i + r_j;
                    if atom_i.xyz.square_distance_from(&atom_j.xyz) < max_dist * max_dist {
                        neighbors.push((atom_j.xyz, expanded_radii_sq[j]));
                    }
                }

                let mut exposed = 0usize;
                for u in &unit_points {
                    let p = atom_i.xyz + (*u * r_i);
                    let buried = neighbors
                        .iter()
                        .any(|&(center_j, r_j_sq)| p.square_distance_from(&center_j) < r_j_sq);
                    if !buried {
                        exposed += 1;
                    }
                }

                4.0 * PI * r_i * r_i * (exposed as f64) / (options.n_points as f64)
            };

            atom_reports.push(AtomSasa {
                path: atom_i.path.clone(),
                name: atom_i.name.clone(),
                area,
            });
            residue_areas[atom_i.residue_idx] += area;
        }

        let residue_reports = residues
            .iter()
            .enumerate()
            .map(|(idx, residue)| ResidueSasa {
                path: residue.path().to_string(),
                comp_id: residue.name.clone(),
                area: residue_areas[idx],
            })
            .collect();

        SasaReport {
            atoms: atom_reports,
            residues: residue_reports,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;
    use crate::format::mmcif::SimpleMmcif;
    use crate::format::pdb::Pdb;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
    }

    fn make_atom(symbol: &str, name: &str, xyz: Position) -> Atom {
        let mut atom = Atom::new_with_pos(symbol, xyz).expect("valid element symbol");
        atom.name = name.to_string();
        atom
    }

    #[test]
    fn test_single_atom_sasa() {
        let mut res = AtomGroup::with_name("ALA");
        res.set_path("/model_1/A/1/".to_string());
        res.set_atom("CA", make_atom("C", "CA", Position::new(0.0, 0.0, 0.0)));

        let opts = SasaOptions::default();
        let report = res.sasa(&opts);

        assert_eq!(report.atoms.len(), 1);
        assert_eq!(report.residues.len(), 1);

        let r_eff = 1.70 + 1.4; // C Bondi radius 1.70 + probe 1.4 = 3.10 Å
        let expected_area = 4.0 * PI * r_eff * r_eff;
        let actual_area = report.atoms[0].area;
        let abs_error = (actual_area - expected_area).abs();
        println!(
            "single_atom (C): actual={:.6}, expected_4pi_r2={:.6}, abs_error={:.6e}",
            actual_area, expected_area, abs_error
        );

        assert!(abs_error < 1e-10);
        assert_eq!(report.atoms[0].path, "/model_1/A/1/CA");
        assert_eq!(report.atoms[0].name, "CA");
        assert_eq!(report.residues[0].path, "/model_1/A/1/");
        assert_eq!(report.residues[0].comp_id, "ALA");
        assert!((report.residues[0].area - expected_area).abs() < 1e-10);

        // Element not in Bondi table falls back to 1.8 Å without error
        let mut fe_res = AtomGroup::with_name("HEM");
        fe_res.set_path("/model_1/A/2/".to_string());
        fe_res.set_atom("FE", make_atom("Fe", "FE", Position::new(0.0, 0.0, 0.0)));
        let fe_report = fe_res.sasa(&opts);
        let expected_fe_area = 4.0 * PI * (1.8 + 1.4_f64).powi(2);
        assert_eq!(fe_report.atoms.len(), 1);
        assert!((fe_report.atoms[0].area - expected_fe_area).abs() < 1e-10);
    }

    #[test]
    fn test_two_atoms_overlap_and_separated() {
        let opts = SasaOptions::default();
        let r_eff = 1.70 + 1.4; // 3.10 Å
        let single_area = 4.0 * PI * r_eff * r_eff;

        // Overlapping: two C atoms separated by d = 2.0 Å (< 2 * 3.10 = 6.20 Å)
        let d_overlap = 2.0;
        let mut overlapping = AtomGroup::with_name("GLY");
        overlapping.set_path("/model_1/A/1/".to_string());
        overlapping.set_atom("C1", make_atom("C", "C1", Position::new(0.0, 0.0, 0.0)));
        overlapping.set_atom(
            "C2",
            make_atom("C", "C2", Position::new(d_overlap, 0.0, 0.0)),
        );

        let overlap_report = overlapping.sasa(&opts);
        assert_eq!(overlap_report.atoms.len(), 2);
        for atom_sasa in &overlap_report.atoms {
            assert!(atom_sasa.area > 0.0);
            assert!(
                atom_sasa.area < single_area,
                "overlapping atom area ({}) must be smaller than isolated area ({})",
                atom_sasa.area,
                single_area
            );
        }
        // Analytical exposed area for each sphere: 4πR² - 2πRh where h = R - d/2
        let h = r_eff - d_overlap / 2.0;
        let analytical_each = 4.0 * PI * r_eff * r_eff - 2.0 * PI * r_eff * h;
        let point_step_area = single_area / (opts.n_points as f64);
        for atom_sasa in &overlap_report.atoms {
            assert!(
                (atom_sasa.area - analytical_each).abs() <= 2.5 * point_step_area,
                "overlapping atom area ({}) should be close to analytical ({})",
                atom_sasa.area,
                analytical_each
            );
        }

        // Separated: two C atoms separated by d = 10.0 Å (> 2 * 3.10 = 6.20 Å)
        let mut separated = AtomGroup::with_name("GLY");
        separated.set_path("/model_1/A/1/".to_string());
        separated.set_atom("C1", make_atom("C", "C1", Position::new(0.0, 0.0, 0.0)));
        separated.set_atom("C2", make_atom("C", "C2", Position::new(10.0, 0.0, 0.0)));

        let sep_report = separated.sasa(&opts);
        assert_eq!(sep_report.atoms.len(), 2);
        assert!((sep_report.atoms[0].area - single_area).abs() < 1e-10);
        assert!((sep_report.atoms[1].area - single_area).abs() < 1e-10);
        assert!((sep_report.residues[0].area - 2.0 * single_area).abs() < 1e-10);
    }

    #[test]
    fn test_hydrogens_ignored() {
        let opts = SasaOptions::default();

        let mut heavy_only = AtomGroup::with_name("GLY");
        heavy_only.set_path("/model_1/A/1/".to_string());
        heavy_only.set_atom("N", make_atom("N", "N", Position::new(0.0, 0.0, 0.0)));
        heavy_only.set_atom("CA", make_atom("C", "CA", Position::new(1.46, 0.0, 0.0)));
        heavy_only.set_atom("C", make_atom("C", "C", Position::new(2.0, 1.4, 0.0)));
        heavy_only.set_atom("O", make_atom("O", "O", Position::new(1.3, 2.4, 0.0)));

        let mut with_hydrogens = heavy_only.clone();
        with_hydrogens.set_atom("H", make_atom("H", "H", Position::new(-0.5, -0.8, 0.0)));
        with_hydrogens.set_atom("HA2", make_atom("H", "HA2", Position::new(1.8, -0.5, 0.9)));
        with_hydrogens.set_atom("HA3", make_atom("H", "HA3", Position::new(1.8, -0.5, -0.9)));

        let report_heavy = heavy_only.sasa(&opts);
        let report_with_h = with_hydrogens.sasa(&opts);

        assert_eq!(report_with_h.atoms.len(), 4);
        assert_eq!(report_with_h, report_heavy);
    }

    #[test]
    fn test_empty_and_hierarchy_levels_and_immutability() {
        let opts = SasaOptions::default();

        // Empty tree
        let empty = AtomGroup::new();
        let empty_before = empty.clone();
        let empty_report = empty.sasa(&opts);
        assert!(empty_report.atoms.is_empty());
        assert!(empty_report.residues.is_empty());
        assert_eq!(empty, empty_before);

        // Build Root -> Model -> Chain -> Residue hierarchy
        let mut residue = AtomGroup::with_name("ALA");
        residue.set_atom("N", make_atom("N", "N", Position::new(0.0, 0.0, 0.0)));
        residue.set_atom("CA", make_atom("C", "CA", Position::new(1.45, 0.0, 0.0)));
        residue.set_atom("C", make_atom("C", "C", Position::new(2.0, 1.42, 0.0)));
        residue.set_atom("O", make_atom("O", "O", Position::new(1.25, 2.39, 0.0)));
        residue.set_atom("CB", make_atom("C", "CB", Position::new(1.98, -0.77, 1.21)));

        let mut chain = AtomGroup::with_name("A");
        chain.set_group("1", residue);

        let mut model = AtomGroup::with_name("model_1");
        model.set_group("A", chain);

        let mut root = AtomGroup::new();
        root.set_group("model_1", model);

        let root_snapshot = root.clone();

        let rep_root = root.sasa(&opts);
        let rep_model = root.get_group("model_1").unwrap().sasa(&opts);
        let rep_chain = root
            .get_group("model_1")
            .unwrap()
            .get_group("A")
            .unwrap()
            .sasa(&opts);
        let rep_res = root
            .get_group("model_1")
            .unwrap()
            .get_group("A")
            .unwrap()
            .get_group("1")
            .unwrap()
            .sasa(&opts);

        assert_eq!(rep_root.atoms.len(), 5);
        assert_eq!(rep_root.residues.len(), 1);
        assert_eq!(rep_root, rep_model);
        assert_eq!(rep_root, rep_chain);
        assert_eq!(rep_root, rep_res);
        assert_eq!(root, root_snapshot, "tree must not be modified by sasa()");
    }

    #[test]
    fn test_real_structures_1hls_and_2fb4() {
        let opts = SasaOptions::default();

        // 1hls.pdb (NMR model 1, 51 residues, 403 heavy atoms + hydrogens)
        let pdb_path = test_data_dir().join("1hls.pdb");
        let pdb = Pdb::from_file(&pdb_path, None).expect("load 1hls.pdb");
        let ag_1hls = pdb.get_atomgroup(None, None).expect("1hls atomgroup");
        let ag_1hls_before = ag_1hls.clone();

        let report_1hls = ag_1hls.sasa(&opts);
        assert_eq!(ag_1hls, ag_1hls_before);
        assert_eq!(report_1hls.atoms.len(), 403);
        assert_eq!(report_1hls.residues.len(), 51);

        let sum_atoms_1hls: f64 = report_1hls.atoms.iter().map(|a| a.area).sum();
        let sum_res_1hls: f64 = report_1hls.residues.iter().map(|r| r.area).sum();
        println!(
            "1hls.pdb: atoms={}, residues={}, total_area={:.6} (res_sum={:.6})",
            report_1hls.atoms.len(),
            report_1hls.residues.len(),
            sum_atoms_1hls,
            sum_res_1hls
        );
        assert!((sum_atoms_1hls - sum_res_1hls).abs() < 1e-9);
        assert!((sum_atoms_1hls - 3992.152817).abs() < 1e-4);

        // Calling on model_1 directly yields the exact same result as root
        let model_1 = ag_1hls.get_group("model_1").expect("model_1");
        assert_eq!(model_1.sasa(&opts), report_1hls);

        // 2FB4.cif (X-ray structure, 550 residues, 3417 heavy atoms)
        let cif_path = test_data_dir().join("2FB4.cif");
        let cif = SimpleMmcif::from_file(&cif_path).expect("load 2FB4.cif");
        let ag_2fb4 = cif.get_atomgroup("2FB4").expect("2FB4 atomgroup");
        let report_2fb4 = ag_2fb4.sasa(&opts);
        assert_eq!(report_2fb4.atoms.len(), 3417);
        assert_eq!(report_2fb4.residues.len(), 550);

        let sum_atoms_2fb4: f64 = report_2fb4.atoms.iter().map(|a| a.area).sum();
        let sum_res_2fb4: f64 = report_2fb4.residues.iter().map(|r| r.area).sum();
        assert!((sum_atoms_2fb4 - sum_res_2fb4).abs() < 1e-9);
        assert!((sum_atoms_2fb4 - 20832.689166).abs() < 1e-4);
    }
}
