// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Geometry inspection (bond lengths, bond angles, aromatic ring planarity).
//!
//! Ideal bond lengths and angles are derived from the idealized coordinates
//! (`ideal_xyz`) of the CCD templates in [`CcdTemplateDb`]. Bond length and bond angle
//! thresholds are ratios of the ideal value; the planarity threshold is in Angstrom.
//! Thresholds must always be given explicitly at this level. Only violations are reported.

use crate::atom_group::AtomGroup;
use crate::brd::sort_nicely;
use crate::ccd_templates::{CcdBondTemplate, CcdTemplateDb};
use crate::ch_pi::{calc_ring_center_and_normal, AROMATIC_RINGS};
use crate::error::{BridgeError, Result};
use crate::position::Position;

/// Kind of geometry violation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeometryViolationKind {
    /// Bond length deviates from the template ideal value.
    BondLength,
    /// Bond angle deviates from the template ideal value.
    BondAngle,
    /// Aromatic ring atom deviates from the best-fit ring plane.
    Planarity,
}

/// A single geometry violation.
#[derive(Debug, Clone, PartialEq)]
pub struct GeometryViolation {
    /// Violation kind.
    pub kind: GeometryViolationKind,
    /// Hierarchical path of the residue group.
    pub path: String,
    /// Component identifier of the residue (e.g. `"ALA"`).
    pub comp_id: String,
    /// Atom names involved: 2 for a bond length, 3 for a bond angle (vertex in the
    /// middle), and the ring atom names for planarity.
    pub atoms: Vec<String>,
    /// Ideal value (Angstrom for lengths, degrees for angles, `0.0` for planarity).
    pub ideal: f64,
    /// Observed value (same unit as `ideal`; for planarity, the maximum out-of-plane distance).
    pub observed: f64,
    /// Absolute deviation `|observed - ideal|`.
    pub deviation: f64,
}

/// List of geometry violations. Atoms without violations are not included.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GeometryReport {
    /// Violations in traversal order.
    pub violations: Vec<GeometryViolation>,
}

impl GeometryReport {
    /// Returns the number of violations.
    pub fn len(&self) -> usize {
        self.violations.len()
    }

    /// Returns `true` if there are no violations.
    pub fn is_empty(&self) -> bool {
        self.violations.is_empty()
    }
}

fn validate_threshold(name: &str, value: f64) -> Result<()> {
    if value.is_nan() || value < 0.0 {
        return Err(BridgeError::value_error(
            name,
            format!("threshold must be a non-negative number, got {value}"),
        ));
    }
    Ok(())
}

fn collect_leaf_residues<'a>(group: &'a AtomGroup, out: &mut Vec<&'a AtomGroup>) {
    if group.atoms().next().is_some() && group.groups().next().is_none() {
        out.push(group);
    } else {
        let mut keys = group.get_group_list();
        sort_nicely(&mut keys);
        for k in keys {
            if let Some(child) = group.get_group(&k) {
                collect_leaf_residues(child, out);
            }
        }
    }
}

fn lookup_template<'a>(db: &'a CcdTemplateDb, comp_id: &str) -> Option<&'a CcdBondTemplate> {
    db.lookup(comp_id).or_else(|| {
        if crate::his_state::is_histidine_comp_id(comp_id) {
            db.lookup("HIS")
        } else {
            None
        }
    })
}

fn ideal_pos(template: &CcdBondTemplate, name: &str) -> Option<Position> {
    let (x, y, z) = template.get_atom(name)?.ideal_xyz?;
    Some(Position::new(x, y, z))
}

/// Returns true when the atom name denotes a hydrogen or deuterium (starts with `H` or `D`).
fn is_hydrogen_name(name: &str) -> bool {
    name.starts_with('H') || name.starts_with('D')
}

fn angle_deg(a: &Position, vertex: &Position, c: &Position) -> Option<f64> {
    let v1 = Position::new(a.x - vertex.x, a.y - vertex.y, a.z - vertex.z);
    let v2 = Position::new(c.x - vertex.x, c.y - vertex.y, c.z - vertex.z);
    let denom = v1.length() * v2.length();
    if denom < 1e-12 {
        return None;
    }
    Some((v1.dot(&v2) / denom).clamp(-1.0, 1.0).acos().to_degrees())
}

/// Checks bond lengths against the template ideal values.
///
/// A bond is reported when `|observed - ideal| / ideal > max_ratio` (e.g. `0.10` for 10%).
/// Bonds whose ideal length is zero or less are skipped, as are bonds whose
/// atoms are absent from the residue, or whose template atoms have no ideal coordinates; residues without a template are skipped. When `include_hydrogens` is
/// false, bonds involving an atom whose name starts with `H` or `D` are not checked.
///
/// # Errors
/// Returns an error if `max_ratio` is negative or NaN.
pub fn check_bond_lengths(
    group: &AtomGroup,
    db: &CcdTemplateDb,
    max_ratio: f64,
    include_hydrogens: bool,
) -> Result<GeometryReport> {
    validate_threshold("max_ratio", max_ratio)?;
    let mut residues = Vec::new();
    collect_leaf_residues(group, &mut residues);

    let mut report = GeometryReport::default();
    for res in residues {
        let Some(template) = lookup_template(db, &res.name) else {
            continue;
        };
        for (n1, n2, _) in &template.bonds {
            if !include_hydrogens && (is_hydrogen_name(n1) || is_hydrogen_name(n2)) {
                continue;
            }
            let (Some(i1), Some(i2)) = (ideal_pos(template, n1), ideal_pos(template, n2)) else {
                continue;
            };
            let (Some(a1), Some(a2)) = (res.get_atom(n1), res.get_atom(n2)) else {
                continue;
            };
            let ideal = i1.distance_from(&i2);
            let observed = a1.xyz.distance_from(&a2.xyz);
            if ideal <= 0.0 {
                continue;
            }
            let deviation = (observed - ideal).abs();
            if deviation / ideal > max_ratio {
                report.violations.push(GeometryViolation {
                    kind: GeometryViolationKind::BondLength,
                    path: res.path().to_string(),
                    comp_id: res.name.clone(),
                    atoms: vec![n1.clone(), n2.clone()],
                    ideal,
                    observed,
                    deviation,
                });
            }
        }
    }
    Ok(report)
}

/// Checks bond angles against the template ideal values.
///
/// For every atom with two or more template neighbors, each neighbor pair forms an angle.
/// An angle is reported when `|observed - ideal| / ideal > max_ratio` (e.g. `0.10` for 10%).
/// Angles whose ideal value is zero or less, with missing atoms, or without ideal
/// coordinates are skipped. When
/// `include_hydrogens` is false, angles involving an atom whose name starts with `H` or
/// `D` are not checked.
///
/// # Errors
/// Returns an error if `max_ratio` is negative or NaN.
pub fn check_bond_angles(
    group: &AtomGroup,
    db: &CcdTemplateDb,
    max_ratio: f64,
    include_hydrogens: bool,
) -> Result<GeometryReport> {
    validate_threshold("max_ratio", max_ratio)?;
    let mut residues = Vec::new();
    collect_leaf_residues(group, &mut residues);

    let mut report = GeometryReport::default();
    for res in residues {
        let Some(template) = lookup_template(db, &res.name) else {
            continue;
        };
        for center in &template.atoms {
            let neighbors: Vec<&String> = template
                .bonds
                .iter()
                .filter_map(|(a, b, _)| {
                    if a == &center.name {
                        Some(b)
                    } else if b == &center.name {
                        Some(a)
                    } else {
                        None
                    }
                })
                .collect();
            for i in 0..neighbors.len() {
                for j in (i + 1)..neighbors.len() {
                    let (n1, n2) = (neighbors[i], neighbors[j]);
                    if !include_hydrogens
                        && [n1.as_str(), center.name.as_str(), n2.as_str()]
                            .iter()
                            .any(|n| is_hydrogen_name(n))
                    {
                        continue;
                    }
                    let (Some(p1), Some(pc), Some(p2)) = (
                        ideal_pos(template, n1),
                        ideal_pos(template, &center.name),
                        ideal_pos(template, n2),
                    ) else {
                        continue;
                    };
                    let (Some(a1), Some(ac), Some(a2)) = (
                        res.get_atom(n1),
                        res.get_atom(&center.name),
                        res.get_atom(n2),
                    ) else {
                        continue;
                    };
                    let (Some(ideal), Some(observed)) = (
                        angle_deg(&p1, &pc, &p2),
                        angle_deg(&a1.xyz, &ac.xyz, &a2.xyz),
                    ) else {
                        continue;
                    };
                    if ideal <= 0.0 {
                        continue;
                    }
                    let deviation = (observed - ideal).abs();
                    if deviation / ideal > max_ratio {
                        report.violations.push(GeometryViolation {
                            kind: GeometryViolationKind::BondAngle,
                            path: res.path().to_string(),
                            comp_id: res.name.clone(),
                            atoms: vec![n1.clone(), center.name.clone(), n2.clone()],
                            ideal,
                            observed,
                            deviation,
                        });
                    }
                }
            }
        }
    }
    Ok(report)
}

/// Checks planarity of aromatic rings (HIS, PHE, TYR, TRP ring atoms).
///
/// The best-fit plane of the ring atoms is computed, and a ring is reported when the
/// maximum out-of-plane distance of its atoms exceeds `max_deviation` (Angstrom).
/// Rings with any missing atom are skipped.
///
/// # Errors
/// Returns an error if `max_deviation` is negative or NaN.
pub fn check_planarity(group: &AtomGroup, max_deviation: f64) -> Result<GeometryReport> {
    validate_threshold("max_deviation", max_deviation)?;
    let mut residues = Vec::new();
    collect_leaf_residues(group, &mut residues);

    let mut report = GeometryReport::default();
    for res in residues {
        let base = if crate::his_state::is_histidine_comp_id(&res.name) {
            "HIS"
        } else {
            res.name.as_str()
        };
        for def in AROMATIC_RINGS.iter().filter(|d| d.residue_name == base) {
            let coords: Option<Vec<Position>> = def
                .atom_names
                .iter()
                .map(|n| res.get_atom(n).map(|a| a.xyz))
                .collect();
            let Some(coords) = coords else {
                continue;
            };
            let (center, normal) = calc_ring_center_and_normal(&coords)?;
            let observed = coords
                .iter()
                .map(|p| {
                    let d = Position::new(p.x - center.x, p.y - center.y, p.z - center.z);
                    d.dot(&normal).abs()
                })
                .fold(0.0_f64, f64::max);
            if observed > max_deviation {
                report.violations.push(GeometryViolation {
                    kind: GeometryViolationKind::Planarity,
                    path: res.path().to_string(),
                    comp_id: res.name.clone(),
                    atoms: def.atom_names.iter().map(|s| s.to_string()).collect(),
                    ideal: 0.0,
                    observed,
                    deviation: observed,
                });
            }
        }
    }
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;

    fn build_residue(comp_id: &str) -> AtomGroup {
        let db = CcdTemplateDb::global();
        let template = db.lookup(comp_id).expect("template should exist");
        let mut res = AtomGroup::with_name(comp_id);
        for (idx, a) in template.atoms.iter().enumerate() {
            let (x, y, z) = a.ideal_xyz.expect("ideal coordinates");
            let mut atom = Atom::new_with_pos(&a.element, Position::new(x, y, z))
                .expect("valid element symbol");
            atom.name = a.name.clone();
            res.set_atom(&format!("{}_{}", idx + 1, a.name), atom);
        }
        res
    }

    fn shift_atom(res: &mut AtomGroup, name: &str, dx: f64, dy: f64, dz: f64) {
        let atom = res
            .atoms_mut()
            .map(|(_, a)| a)
            .find(|a| a.name == name)
            .expect("atom exists");
        atom.xyz = Position::new(atom.xyz.x + dx, atom.xyz.y + dy, atom.xyz.z + dz);
    }

    #[test]
    fn ideal_residue_has_no_violations() {
        let db = CcdTemplateDb::global();
        let res = build_residue("ALA");
        assert!(check_bond_lengths(&res, db, 0.001, true)
            .unwrap()
            .is_empty());
        assert!(check_bond_angles(&res, db, 0.001, true).unwrap().is_empty());
        let phe = build_residue("PHE");
        assert!(check_planarity(&phe, 0.01).unwrap().is_empty());
    }

    #[test]
    fn stretched_bond_is_reported_only_beyond_threshold() {
        let db = CcdTemplateDb::global();
        let mut res = build_residue("ALA");
        // Move the terminal CB along its bond direction by exactly 0.2 A (CA-CB lengthens).
        let ca = res.get_atom("CA").unwrap().xyz;
        let cb = res.get_atom("CB").unwrap().xyz;
        let len = ca.distance_from(&cb);
        let s = 0.2 / len;
        shift_atom(
            &mut res,
            "CB",
            (cb.x - ca.x) * s,
            (cb.y - ca.y) * s,
            (cb.z - ca.z) * s,
        );

        // CA-CB lengthens by 0.2 A (about 13% of the ideal length).
        // ratio 0.05: CA-CB reported; ratio 0.2: nothing.
        let tight = check_bond_lengths(&res, db, 0.05, true).unwrap();
        assert_eq!(tight.len(), 1);
        let v = &tight.violations[0];
        assert_eq!(v.kind, GeometryViolationKind::BondLength);
        assert_eq!(v.atoms, vec!["CA".to_string(), "CB".to_string()]);
        assert!((v.deviation - 0.2).abs() < 1e-6);
        assert!(check_bond_lengths(&res, db, 0.2, true).unwrap().is_empty());
    }

    #[test]
    fn distorted_angle_is_reported() {
        let db = CcdTemplateDb::global();
        let mut res = build_residue("ALA");
        // Displace CB perpendicular-ish to change N-CA-CB and CB-CA-C angles.
        shift_atom(&mut res, "CB", 0.0, 0.0, 0.6);
        let loose = check_bond_angles(&res, db, 0.01, true).unwrap();
        assert!(!loose.is_empty());
        assert!(loose
            .violations
            .iter()
            .all(|v| v.kind == GeometryViolationKind::BondAngle && v.atoms.len() == 3));
        assert!(loose
            .violations
            .iter()
            .any(|v| v.atoms[2] == "CB" || v.atoms[0] == "CB"));
        // A huge threshold reports nothing.
        assert!(check_bond_angles(&res, db, 100.0, true).unwrap().is_empty());
        // The number of violations does not increase when the threshold is relaxed.
        let mid = check_bond_angles(&res, db, 0.05, true).unwrap();
        assert!(mid.len() <= loose.len());
    }

    #[test]
    fn hydrogen_angles_can_be_excluded() {
        let db = CcdTemplateDb::global();
        let mut res = build_residue("ALA");
        // Move the amide hydrogen so that H-N-CA changes, and CB so that a heavy angle changes.
        shift_atom(&mut res, "H", 0.0, 0.0, 0.8);
        shift_atom(&mut res, "CB", 0.0, 0.0, 0.6);
        let has_h = |v: &GeometryViolation| v.atoms.iter().any(|n| is_hydrogen_name(n));

        let with_h = check_bond_angles(&res, db, 0.01, true).unwrap();
        assert!(with_h.violations.iter().any(has_h));

        let without_h = check_bond_angles(&res, db, 0.01, false).unwrap();
        assert!(!without_h.is_empty());
        assert!(!without_h.violations.iter().any(has_h));

        // Heavy-atom angles are identical in both modes.
        let heavy_of_with_h: Vec<_> = with_h.violations.iter().filter(|v| !has_h(v)).collect();
        assert_eq!(heavy_of_with_h.len(), without_h.len());
        for (a, b) in heavy_of_with_h.iter().zip(without_h.violations.iter()) {
            assert_eq!(a.atoms, b.atoms);
            assert!((a.observed - b.observed).abs() < 1e-12);
        }
    }

    #[test]
    fn hydrogen_bonds_can_be_excluded() {
        let db = CcdTemplateDb::global();
        let mut res = build_residue("ALA");
        shift_atom(&mut res, "H", 0.0, 0.0, 0.8);
        let has_h = |v: &GeometryViolation| v.atoms.iter().any(|n| is_hydrogen_name(n));
        let with_h = check_bond_lengths(&res, db, 0.1, true).unwrap();
        assert!(with_h.violations.iter().any(has_h));
        let without_h = check_bond_lengths(&res, db, 0.1, false).unwrap();
        assert!(!without_h.violations.iter().any(has_h));
    }

    #[test]
    fn out_of_plane_ring_atom_is_reported() {
        let mut res = build_residue("PHE");
        // Move CZ along the ring normal.
        let coords: Vec<Position> = AROMATIC_RINGS
            .iter()
            .find(|d| d.residue_name == "PHE")
            .unwrap()
            .atom_names
            .iter()
            .map(|n| res.get_atom(n).unwrap().xyz)
            .collect();
        let (_, normal) = calc_ring_center_and_normal(&coords).unwrap();
        shift_atom(
            &mut res,
            "CZ",
            normal.x * 0.5,
            normal.y * 0.5,
            normal.z * 0.5,
        );

        let tight = check_planarity(&res, 0.1).unwrap();
        assert_eq!(tight.len(), 1);
        assert_eq!(tight.violations[0].kind, GeometryViolationKind::Planarity);
        assert!(tight.violations[0].observed > 0.1);
        assert!(check_planarity(&res, 1.0).unwrap().is_empty());
    }

    #[test]
    fn invalid_threshold_is_an_error() {
        let db = CcdTemplateDb::global();
        let res = build_residue("ALA");
        assert!(check_bond_lengths(&res, db, -0.1, true).is_err());
        assert!(check_bond_lengths(&res, db, f64::NAN, true).is_err());
        assert!(check_bond_angles(&res, db, -0.1, true).is_err());
        assert!(check_bond_angles(&res, db, f64::NAN, true).is_err());
        assert!(check_planarity(&res, -0.1).is_err());
        assert!(check_planarity(&res, f64::NAN).is_err());
    }
}
