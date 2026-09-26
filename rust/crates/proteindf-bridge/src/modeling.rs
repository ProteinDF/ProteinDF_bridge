// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

#![allow(non_snake_case)]

use std::collections::HashMap;
use std::path::Path;

use regex::Regex;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::brd::{load_atomgroup, load_atomgroup_from_bytes};
use crate::error::{BridgeError, Result};
use crate::matrix::Matrix;
use crate::position::Position;
use crate::superposer::Superposer;

const TRANS1_BYTES: &[u8] = include_bytes!("../tests/data/ACE_ALA_NME_trans1.brd");
const TRANS2_BYTES: &[u8] = include_bytes!("../tests/data/ACE_ALA_NME_trans2.brd");
const CIS1_BYTES: &[u8] = include_bytes!("../tests/data/ACE_ALA_NME_cis1.brd");
const CIS2_BYTES: &[u8] = include_bytes!("../tests/data/ACE_ALA_NME_cis2.brd");

/// Structural modeling and capping/neutralization utility corresponding to `proteindf_bridge.modeling.Modeling`.
#[derive(Debug, Clone)]
pub struct Modeling {
    ace_ala_nme: HashMap<String, AtomGroup>,
}

impl Default for Modeling {
    fn default() -> Self {
        Self::new().expect("failed to initialize embedded ACE_ALA_NME reference conformers")
    }
}

impl Modeling {
    /// List of standard ACE-ALA-NME conformer names.
    pub const CONFORMERS: [&'static str; 4] = ["trans1", "trans2", "cis1", "cis2"];

    /// Creates a new `Modeling` instance with embedded reference conformer structures.
    pub fn new() -> Result<Self> {
        let mut map = HashMap::new();
        map.insert(
            "trans1".to_string(),
            load_atomgroup_from_bytes(TRANS1_BYTES)?,
        );
        map.insert(
            "trans2".to_string(),
            load_atomgroup_from_bytes(TRANS2_BYTES)?,
        );
        map.insert("cis1".to_string(), load_atomgroup_from_bytes(CIS1_BYTES)?);
        map.insert("cis2".to_string(), load_atomgroup_from_bytes(CIS2_BYTES)?);
        Ok(Self { ace_ala_nme: map })
    }

    /// Creates a new `Modeling` instance loading reference conformers from a directory.
    pub fn from_data_dir<P: AsRef<Path>>(data_dir: P) -> Result<Self> {
        let dir = data_dir.as_ref();
        let mut map = HashMap::new();
        for conformer in Self::CONFORMERS {
            let path = dir.join(format!("ACE_ALA_NME_{conformer}.brd"));
            let ag = load_atomgroup(&path)?;
            map.insert(conformer.to_string(), ag);
        }
        Ok(Self { ace_ala_nme: map })
    }

    /// Returns a reference to the stored conformer `AtomGroup`.
    pub fn get_ace_ala_nme(&self, conformer: &str) -> Option<&AtomGroup> {
        self.ace_ala_nme.get(conformer)
    }

    // -----------------------------------------------------------------
    // Capping: ACE and NME
    // -----------------------------------------------------------------

    /// Evaluates all reference conformers using `match_fn`, tracking the best RMSD match
    /// and collecting errors so full context is preserved if all conformers fail.
    fn find_best_conformer<F>(
        &self,
        capping_name: &str,
        res: &AtomGroup,
        next_aa: Option<&AtomGroup>,
        match_fn: F,
    ) -> Result<AtomGroup>
    where
        F: Fn(&Self, &AtomGroup, &AtomGroup, Option<&AtomGroup>) -> Result<(AtomGroup, f64)>,
    {
        let mut aan_best: Option<AtomGroup> = None;
        let mut rmsd_min = f64::MAX;
        let mut errors: Vec<(String, String)> = Vec::new();

        for conformer in Self::CONFORMERS {
            let ref_aan = match self.ace_ala_nme.get(conformer) {
                Some(aan) => aan,
                None => {
                    errors.push((
                        conformer.to_string(),
                        format!("reference conformer {conformer} not found"),
                    ));
                    continue;
                }
            };
            match match_fn(self, ref_aan, res, next_aa) {
                Ok((matched, rmsd)) => {
                    if rmsd < rmsd_min {
                        rmsd_min = rmsd;
                        aan_best = Some(matched);
                    }
                }
                Err(err) => {
                    log::warn!("{capping_name} conformer {conformer} match failed: {err}");
                    errors.push((conformer.to_string(), err.to_string()));
                }
            }
        }

        match aan_best {
            Some(best) => {
                if rmsd_min > 1.0 {
                    log::warn!("{capping_name} RMSD value is too large: {rmsd_min}");
                }
                Ok(best)
            }
            None => {
                let err_details = errors
                    .into_iter()
                    .map(|(conf, err)| format!("{conf}: {err}"))
                    .collect::<Vec<_>>()
                    .join("; ");
                Err(BridgeError::general(format!(
                    "no matching {capping_name} conformer found ({err_details})"
                )))
            }
        }
    }

    /// Builds an ACE capping group by fitting against the optimal reference conformer.
    pub fn get_ACE(&self, res: &AtomGroup, next_aa: Option<&AtomGroup>) -> Result<AtomGroup> {
        let best = self.find_best_conformer("ACE", res, next_aa, Self::match_ace)?;
        let ace_group = best.get_group("1").ok_or_else(|| {
            BridgeError::general("ACE group '1' not found in reference structure")
        })?;
        let mut answer = ace_group.clone();
        answer.set_path_with_depth("/ACE".to_string(), 1);
        Ok(answer)
    }

    fn match_ace(
        &self,
        aan: &AtomGroup,
        res: &AtomGroup,
        next_aa: Option<&AtomGroup>,
    ) -> Result<(AtomGroup, f64)> {
        let aan_res2 = aan
            .get_group("2")
            .ok_or_else(|| BridgeError::general("reference group '2' not found in ACE-ALA-NME"))?;
        let (mut aan_part, mut res_part) = self.match_residues(aan_res2, res, -1);

        if let Some(next) = next_aa {
            let aan_res3 = aan.get_group("3").ok_or_else(|| {
                BridgeError::general("reference group '3' not found in ACE-ALA-NME")
            })?;
            if next.has_atom("N") {
                if let Some(n) = aan_res3.get_atom("N") {
                    aan_part.set_atom("N2", n.clone());
                    res_part.set_atom("N2", next.get_atom("N").unwrap().clone());
                }
            }
            if next.has_atom("H") {
                if let Some(h) = aan_res3.get_atom("H") {
                    aan_part.set_atom("NH2", h.clone());
                    res_part.set_atom("NH2", next.get_atom("H").unwrap().clone());
                }
            }
            if next.has_atom("CA") {
                if let Some(ca) = aan_res3.get_atom("CH3") {
                    aan_part.set_atom("CH3", ca.clone());
                    res_part.set_atom("CH3", next.get_atom("CA").unwrap().clone());
                }
            }
        }

        let sp = Superposer::new(&aan_part, &res_part)?;
        let rmsd = sp.rmsd();
        let matched_aan = sp.superimpose(aan)?;
        Ok((matched_aan, rmsd))
    }

    /// Builds an NME capping group by fitting against the optimal reference conformer.
    pub fn get_NME(&self, res: &AtomGroup, next_aa: Option<&AtomGroup>) -> Result<AtomGroup> {
        let best = self.find_best_conformer("NME", res, next_aa, Self::match_nme)?;
        let nme_group = best.get_group("3").ok_or_else(|| {
            BridgeError::general("NME group '3' not found in reference structure")
        })?;
        let mut answer = nme_group.clone();
        answer.set_path_with_depth("/NME".to_string(), 1);
        Ok(answer)
    }

    fn match_nme(
        &self,
        aan: &AtomGroup,
        res: &AtomGroup,
        next_aa: Option<&AtomGroup>,
    ) -> Result<(AtomGroup, f64)> {
        let aan_res2 = aan
            .get_group("2")
            .ok_or_else(|| BridgeError::general("reference group '2' not found in ACE-ALA-NME"))?;
        let (mut aan_part, mut res_part) = self.match_residues(aan_res2, res, -1);

        if let Some(next) = next_aa {
            let aan_res1 = aan.get_group("1").ok_or_else(|| {
                BridgeError::general("reference group '1' not found in ACE-ALA-NME")
            })?;
            if next.has_atom("C") {
                if let Some(c) = aan_res1.get_atom("C") {
                    aan_part.set_atom("C2", c.clone());
                    res_part.set_atom("C2", next.get_atom("C").unwrap().clone());
                }
            }
            if next.has_atom("O") {
                if let Some(o) = aan_res1.get_atom("O") {
                    aan_part.set_atom("O2", o.clone());
                    res_part.set_atom("O2", next.get_atom("O").unwrap().clone());
                }
            }
            if next.has_atom("CA") {
                if let Some(ca) = aan_res1.get_atom("CH3") {
                    aan_part.set_atom("CH3", ca.clone());
                    res_part.set_atom("CH3", next.get_atom("CA").unwrap().clone());
                }
            }
        }

        let sp = Superposer::new(&aan_part, &res_part)?;
        let rmsd = sp.rmsd();
        let matched_aan = sp.superimpose(aan)?;
        Ok((matched_aan, rmsd))
    }

    fn match_residues(
        &self,
        res1: &AtomGroup,
        res2: &AtomGroup,
        mut max_number_of_atoms: i32,
    ) -> (AtomGroup, AtomGroup) {
        let atom_names = ["CA", "O", "C", "N", "CB", "HA"];
        if max_number_of_atoms == -1 {
            max_number_of_atoms = atom_names.len() as i32;
        }
        let mut ans_res1 = AtomGroup::new();
        let mut ans_res2 = AtomGroup::new();

        for atom_name in atom_names {
            let pickup1 = res1.pickup_atoms(atom_name);
            if !pickup1.is_empty() {
                let pickup2 = res2.pickup_atoms(atom_name);
                if !pickup2.is_empty() {
                    ans_res1.set_atom(atom_name, pickup1[0].clone());
                    ans_res2.set_atom(atom_name, pickup2[0].clone());
                }
            }
            if ans_res1.get_number_of_atoms() >= max_number_of_atoms as usize {
                break;
            }
        }

        if (ans_res1.get_number_of_atoms() as i32) < max_number_of_atoms {
            let res1_h = if let Some(h) = res1.get_atom("H") {
                Some(h.clone())
            } else {
                res1.get_atom("CD").cloned()
            };
            let res2_h = if let Some(h) = res2.get_atom("H") {
                Some(h.clone())
            } else {
                res2.get_atom("CD").cloned()
            };
            if let (Some(h1), Some(h2)) = (res1_h, res2_h) {
                ans_res1.set_atom("H", h1);
                ans_res2.set_atom("H", h2);
            }
        }

        (ans_res1, ans_res2)
    }

    /// Turns the neighboring C-alpha position into an ACE methyl group.
    pub fn get_ACE_simple(&self, next_aa: &AtomGroup) -> Result<AtomGroup> {
        let mut answer = AtomGroup::new();

        let cas = next_aa.pickup_atoms("CA");
        if !cas.is_empty() {
            answer.set_atom("CA", cas[0].clone());
        } else {
            return Err(BridgeError::input_error(
                "next_aa",
                "cannot found \"CA\" atom on building ACE.",
            ));
        }

        let cs = next_aa.pickup_atoms("C");
        if !cs.is_empty() {
            answer.set_atom("C", cs[0].clone());
        } else {
            return Err(BridgeError::input_error(
                "next_aa",
                "cannot found \"C\" atom on building ACE.",
            ));
        }

        let os = next_aa.pickup_atoms("O");
        if !os.is_empty() {
            answer.set_atom("O", os[0].clone());
        } else {
            return Err(BridgeError::input_error(
                "next_aa",
                "cannot found \"O\" atom on building ACE.",
            ));
        }

        let ca = answer.get_atom("CA").unwrap().clone();
        let c = answer.get_atom("C").unwrap().clone();
        let methyl = self.add_methyl(&ca, &c)?;
        answer |= methyl;
        answer.set_path("/ACE".to_string());
        Ok(answer)
    }

    /// Turns the neighboring C-alpha position into an NME methyl group.
    pub fn get_NME_simple(&self, next_aa: &AtomGroup) -> Result<AtomGroup> {
        let mut answer = AtomGroup::new();

        let cas = next_aa.pickup_atoms("CA");
        if !cas.is_empty() {
            answer.set_atom("CA", cas[0].clone());
        } else {
            return Err(BridgeError::input_error(
                "next_aa",
                "cannot found \"CA\" atom on building NME.",
            ));
        }

        let ns = next_aa.pickup_atoms("N");
        if !ns.is_empty() {
            answer.set_atom("N", ns[0].clone());
        } else {
            return Err(BridgeError::input_error(
                "next_aa",
                "cannot found \"N\" atom on building NME.",
            ));
        }

        let hs = next_aa.pickup_atoms("H");
        if !hs.is_empty() {
            answer.set_atom("H", hs[0].clone());
        } else {
            let cds = next_aa.pickup_atoms("CD");
            if !cds.is_empty() {
                let mut dummy_h = cds[0].clone();
                dummy_h.set_symbol("H")?;
                dummy_h.name = "H".to_string();
                answer.set_atom("H", dummy_h);
            } else {
                return Err(BridgeError::input_error(
                    "next_aa",
                    "cannot found \"H\" or \"CD\" atom(for proline) on building NME.",
                ));
            }
        }

        let ca = answer.get_atom("CA").unwrap().clone();
        let n = answer.get_atom("N").unwrap().clone();
        let methyl = self.add_methyl(&ca, &n)?;
        answer |= methyl;
        answer.set_path("/NME".to_string());
        Ok(answer)
    }

    // -----------------------------------------------------------------
    // Geometry helpers
    // -----------------------------------------------------------------

    /// Adds hydrogens of -CH3 attached to C1, oriented away from C2.
    pub fn add_methyl(&self, c1: &Atom, c2: &Atom) -> Result<AtomGroup> {
        let mut ethane = AtomGroup::new();
        ethane.set_atom(
            "C1",
            Atom::new_with_pos("C", Position::new(0.00000, 0.00000, 0.00000))?,
        );
        ethane.set_atom(
            "H11",
            Atom::new_with_pos("H", Position::new(-0.85617, -0.58901, -0.35051))?,
        );
        ethane.set_atom(
            "H12",
            Atom::new_with_pos("H", Position::new(-0.08202, 1.03597, -0.35051))?,
        );
        ethane.set_atom(
            "H13",
            Atom::new_with_pos("H", Position::new(0.93818, -0.44696, -0.35051))?,
        );
        ethane.set_atom(
            "C2",
            Atom::new_with_pos("C", Position::new(0.00000, 0.00000, 1.47685))?,
        );
        ethane.set_atom(
            "H21",
            Atom::new_with_pos("H", Position::new(-0.93818, 0.44696, 1.82736))?,
        );
        ethane.set_atom(
            "H22",
            Atom::new_with_pos("H", Position::new(0.85617, 0.58901, 1.82736))?,
        );
        ethane.set_atom(
            "H23",
            Atom::new_with_pos("H", Position::new(0.08202, -1.03597, 1.82736))?,
        );

        let inc21 = c2.xyz - c1.xyz;
        let ref_c2 = ethane.get_atom("C2").unwrap().xyz;
        let ref_c1 = ethane.get_atom("C1").unwrap().xyz;
        let refc21 = ref_c2 - ref_c1;

        let shift = c1.xyz - ref_c1;
        let rot = self.arbitary_rotate_matrix(inc21, refc21)?;

        ethane.rotate(&rot)?;
        ethane.shift_by(shift);

        let mut answer = AtomGroup::new();
        answer.set_atom("H11", ethane.get_atom("H11").unwrap().clone());
        answer.set_atom("H12", ethane.get_atom("H12").unwrap().clone());
        answer.set_atom("H13", ethane.get_atom("H13").unwrap().clone());

        Ok(answer)
    }

    /// Computes the 3x3 rotation matrix that aligns vector `in_a` with `in_b`.
    pub fn arbitary_rotate_matrix(&self, in_a: Position, in_b: Position) -> Result<Matrix> {
        let mut a = in_a;
        let mut b = in_b;
        a.norm()?;
        b.norm()?;

        let cos_theta = a.dot(&b);
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();

        let mut n = a.cross(&b);
        n.norm()?;

        let nx = n.x;
        let ny = n.y;
        let nz = n.z;

        let mut rot = Matrix::new(3, 3);
        rot.set(0, 0, nx * nx * (1.0 - cos_theta) + cos_theta);
        rot.set(0, 1, nx * ny * (1.0 - cos_theta) + nz * sin_theta);
        rot.set(0, 2, nx * nz * (1.0 - cos_theta) - ny * sin_theta);
        rot.set(1, 0, nx * ny * (1.0 - cos_theta) - nz * sin_theta);
        rot.set(1, 1, ny * ny * (1.0 - cos_theta) + cos_theta);
        rot.set(1, 2, nx * nz * (1.0 - cos_theta) + nx * sin_theta);
        rot.set(2, 0, nx * nz * (1.0 - cos_theta) + ny * sin_theta);
        rot.set(2, 1, ny * nz * (1.0 - cos_theta) - nx * sin_theta);
        rot.set(2, 2, nz * nz * (1.0 - cos_theta) + cos_theta);

        Ok(rot)
    }

    /// Generates an NH3 atom group with given angle and bond length.
    pub fn get_NH3(&self, angle: f64, length: f64) -> Result<AtomGroup> {
        let pi23 = std::f64::consts::PI * 2.0 / 3.0;
        let sin23 = pi23.sin();
        let cos23 = pi23.cos();
        let sin_input = angle.sin();
        let cos_input = angle.cos();

        let mut xz_rot = Matrix::new(3, 3);
        xz_rot.set(0, 0, cos_input);
        xz_rot.set(0, 2, -sin_input);
        xz_rot.set(2, 0, sin_input);
        xz_rot.set(2, 2, cos_input);
        xz_rot.set(1, 1, 1.0);

        let mut xy_rot = Matrix::new(3, 3);
        xy_rot.set(0, 0, cos23);
        xy_rot.set(0, 1, -sin23);
        xy_rot.set(1, 0, sin23);
        xy_rot.set(1, 1, cos23);
        xy_rot.set(2, 2, 1.0);

        let mut pos_h1 = Position::new(0.0, 0.0, 1.0);
        pos_h1.rotate(&xz_rot)?;

        let mut pos_h2 = Position::new(0.0, 0.0, 1.0);
        pos_h2.rotate(&xz_rot)?;
        pos_h2.rotate(&xy_rot)?;

        let mut pos_h3 = Position::new(0.0, 0.0, 1.0);
        pos_h3.rotate(&xz_rot)?;
        pos_h3.rotate(&xy_rot)?;
        pos_h3.rotate(&xy_rot)?;

        pos_h1 *= length;
        pos_h2 *= length;
        pos_h3 *= length;

        let mut nh3 = AtomGroup::new();
        nh3.set_atom("N", Atom::new_with_pos("N", Position::new(0.0, 0.0, 0.0))?);
        nh3.set_atom("H1", Atom::new_with_pos("H", pos_h1)?);
        nh3.set_atom("H2", Atom::new_with_pos("H", pos_h2)?);
        nh3.set_atom("H3", Atom::new_with_pos("H", pos_h3)?);

        Ok(nh3)
    }

    /// Selects consecutive amino acid residues from a chain group.
    pub fn select_residues(&self, chain: &AtomGroup, from_resid: i64, to_resid: i64) -> AtomGroup {
        let mut answer = AtomGroup::new();
        for (resid_key, res) in chain.groups() {
            if let Ok(resid_int) = resid_key.parse::<i64>() {
                if resid_int >= from_resid && resid_int <= to_resid {
                    answer |= res.clone();
                }
            }
        }
        answer
    }

    /// Returns the maximum numerical index found in atom keys.
    pub fn get_last_index(&self, res: &AtomGroup) -> usize {
        let mut answer = 0;
        let re = Regex::new("([0-9]+)").unwrap();
        for (key, _) in res.atoms() {
            if let Some(m) = re.find(key) {
                if let Ok(num) = m.as_str().parse::<usize>() {
                    answer = answer.max(num);
                }
            }
        }
        answer
    }

    // -----------------------------------------------------------------
    // Neutralization helpers
    // -----------------------------------------------------------------

    /// Returns a Cl- group to neutralize the N-terminal residue.
    pub fn neutralize_Nterm(&self, res: &AtomGroup) -> Result<AtomGroup> {
        if res.name == "PRO" {
            self.neutralize_nterm_pro(res)
        } else {
            self.neutralize_nterm_general(res)
        }
    }

    fn neutralize_nterm_general(&self, res: &AtomGroup) -> Result<AtomGroup> {
        let mut ag = AtomGroup::new();
        let n = res
            .get_atom("N")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find N"))?;
        let h1 = res
            .get_atom("H1")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find H1"))?;
        let h2 = res
            .get_atom("H2")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find H2"))?;
        ag.set_atom("N", n.clone());
        ag.set_atom("H1", h1.clone());
        ag.set_atom("H2", h2.clone());

        if let Some(hxt) = res.get_atom("HXT") {
            ag.set_atom("H3", hxt.clone());
        } else if let Some(h3) = res.get_atom("H3") {
            ag.set_atom("H3", h3.clone());
        } else {
            return Err(BridgeError::input_error("res", "cannot find HXT or H3"));
        }

        let pos = self.get_neutralize_pos_nh3_type(&ag)?;
        let mut answer = AtomGroup::new();
        let mut cl = Atom::new_with_pos("Cl", pos)?;
        cl.name = "Cl".to_string();
        answer.set_atom("Cl", cl);
        Ok(answer)
    }

    fn neutralize_nterm_pro(&self, res: &AtomGroup) -> Result<AtomGroup> {
        let mut ag = AtomGroup::new();
        let n = res
            .get_atom("N")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find N"))?;
        let h2 = res
            .get_atom("H2")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find H2"))?;
        ag.set_atom("N", n.clone());
        ag.set_atom("H2", h2.clone());

        if let Some(hxt) = res.get_atom("HXT") {
            ag.set_atom("H1", hxt.clone());
        } else if let Some(h3) = res.get_atom("H3") {
            ag.set_atom("H1", h3.clone());
        } else {
            return Err(BridgeError::input_error("res", "cannot find HXT or H3"));
        }

        let pos = self.get_neutralize_pos_nh2_type(&ag)?;
        let mut answer = AtomGroup::new();
        let mut cl = Atom::new_with_pos("Cl", pos)?;
        cl.name = "Cl".to_string();
        answer.set_atom("Cl", cl);
        Ok(answer)
    }

    /// Returns an Na+ group to neutralize the C-terminal residue.
    pub fn neutralize_Cterm(&self, res: &AtomGroup) -> Result<AtomGroup> {
        let mut ag = AtomGroup::new();
        let c = res
            .get_atom("C")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find C"))?;
        let o = res
            .get_atom("O")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find O"))?;
        let oxt = res
            .get_atom("OXT")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find OXT"))?;
        ag.set_atom("C", c.clone());
        ag.set_atom("O1", o.clone());
        ag.set_atom("O2", oxt.clone());

        let pos = self.get_neutralize_pos_coo_type(&ag)?;
        let mut answer = AtomGroup::new();
        let mut na = Atom::new_with_pos("Na", pos)?;
        na.name = "Na".to_string();
        answer.set_atom("Na", na);
        Ok(answer)
    }

    /// Returns an Na+ group to neutralize a GLU side chain.
    pub fn neutralize_GLU(&self, res: &AtomGroup) -> Result<AtomGroup> {
        let mut ag = AtomGroup::new();
        let c = res
            .get_atom("CD")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find CD"))?;
        let o1 = res
            .get_atom("OE1")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find OE1"))?;
        let o2 = res
            .get_atom("OE2")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find OE2"))?;
        ag.set_atom("C", c.clone());
        ag.set_atom("O1", o1.clone());
        ag.set_atom("O2", o2.clone());

        let pos = self.get_neutralize_pos_coo_type(&ag)?;
        let mut answer = AtomGroup::new();
        let key = self.get_last_index(res);
        let mut na = Atom::new_with_pos("Na", pos)?;
        na.name = "Na".to_string();
        answer.set_atom(&format!("{}_Na", key + 1), na);
        Ok(answer)
    }

    /// Returns an Na+ group to neutralize an ASP side chain.
    pub fn neutralize_ASP(&self, res: &AtomGroup) -> Result<AtomGroup> {
        let mut ag = AtomGroup::new();
        let c = res
            .get_atom("CG")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find CG"))?;
        let o1 = res
            .get_atom("OD1")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find OD1"))?;
        let o2 = res
            .get_atom("OD2")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find OD2"))?;
        ag.set_atom("C", c.clone());
        ag.set_atom("O1", o1.clone());
        ag.set_atom("O2", o2.clone());

        let pos = self.get_neutralize_pos_coo_type(&ag)?;
        let mut answer = AtomGroup::new();
        let key = self.get_last_index(res);
        let mut na = Atom::new_with_pos("Na", pos)?;
        na.name = "Na".to_string();
        answer.set_atom(&format!("{}_Na", key + 1), na);
        Ok(answer)
    }

    /// Returns a Cl- group to neutralize a LYS side chain.
    pub fn neutralize_LYS(&self, res: &AtomGroup) -> Result<AtomGroup> {
        let mut ag = AtomGroup::new();
        let n = res
            .get_atom("NZ")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find NZ"))?;
        let h1 = res
            .get_atom("HZ1")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find HZ1"))?;
        let h2 = res
            .get_atom("HZ2")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find HZ2"))?;
        let h3 = res
            .get_atom("HZ3")
            .ok_or_else(|| BridgeError::input_error("res", "cannot find HZ3"))?;
        ag.set_atom("N", n.clone());
        ag.set_atom("H1", h1.clone());
        ag.set_atom("H2", h2.clone());
        ag.set_atom("H3", h3.clone());

        let pos = self.get_neutralize_pos_nh3_type(&ag)?;
        let mut answer = AtomGroup::new();
        let key = self.get_last_index(res);
        let mut cl = Atom::new_with_pos("Cl", pos)?;
        cl.name = "Cl".to_string();
        answer.set_atom(&format!("{}_Cl", key + 1), cl);
        Ok(answer)
    }

    /// Returns a Cl- group to neutralize an ARG side chain.
    /// case: 0 (center), 1 (NH1 side), 2 (NH2 side).
    pub fn neutralize_ARG(&self, res: &AtomGroup, case: usize) -> Result<AtomGroup> {
        let pos = match case {
            0 => {
                let length = 3.0;
                let nh1 = res
                    .get_atom("NH1")
                    .ok_or_else(|| BridgeError::input_error("res", "cannot find NH1"))?;
                let nh2 = res
                    .get_atom("NH2")
                    .ok_or_else(|| BridgeError::input_error("res", "cannot find NH2"))?;
                let cz = res
                    .get_atom("CZ")
                    .ok_or_else(|| BridgeError::input_error("res", "cannot find CZ"))?;
                let m = Position::new(
                    0.5 * (nh1.xyz.x + nh2.xyz.x),
                    0.5 * (nh1.xyz.y + nh2.xyz.y),
                    0.5 * (nh1.xyz.z + nh2.xyz.z),
                );
                let mut v_cm = m - cz.xyz;
                v_cm.norm()?;
                cz.xyz + v_cm * length
            }
            1 => {
                let length = 2.0;
                let hh11 = res
                    .get_atom("HH11")
                    .ok_or_else(|| BridgeError::input_error("res", "cannot find HH11"))?;
                let hh12 = res
                    .get_atom("HH12")
                    .ok_or_else(|| BridgeError::input_error("res", "cannot find HH12"))?;
                let n = res
                    .get_atom("NH1")
                    .ok_or_else(|| BridgeError::input_error("res", "cannot find NH1"))?;
                let m = Position::new(
                    0.5 * (hh11.xyz.x + hh12.xyz.x),
                    0.5 * (hh11.xyz.y + hh12.xyz.y),
                    0.5 * (hh11.xyz.z + hh12.xyz.z),
                );
                let mut v_nm = m - n.xyz;
                v_nm.norm()?;
                n.xyz + v_nm * length
            }
            2 => {
                let length = 2.0;
                let hh21 = res
                    .get_atom("HH21")
                    .ok_or_else(|| BridgeError::input_error("res", "cannot find HH21"))?;
                let hh22 = res
                    .get_atom("HH22")
                    .ok_or_else(|| BridgeError::input_error("res", "cannot find HH22"))?;
                let n = res
                    .get_atom("NH2")
                    .ok_or_else(|| BridgeError::input_error("res", "cannot find NH2"))?;
                let m = Position::new(
                    0.5 * (hh21.xyz.x + hh22.xyz.x),
                    0.5 * (hh21.xyz.y + hh22.xyz.y),
                    0.5 * (hh21.xyz.z + hh22.xyz.z),
                );
                let mut v_nm = m - n.xyz;
                v_nm.norm()?;
                n.xyz + v_nm * length
            }
            other => {
                return Err(BridgeError::value_error(
                    "case",
                    format!("invalid ARG case: {other} (expected 0, 1, or 2)"),
                ));
            }
        };

        let mut answer = AtomGroup::new();
        let key = self.get_last_index(res);
        let mut cl = Atom::new_with_pos("Cl", pos)?;
        cl.name = "Cl".to_string();
        answer.set_atom(&format!("{}_Cl", key + 1), cl);
        Ok(answer)
    }

    /// Neutralizes FAD phosphate groups by adding two Na+ ions.
    /// Fails if OP1/O1P or OP2/O2P are absent (reproducing Python raise behavior).
    pub fn neutralize_FAD(&self, ag: &AtomGroup) -> Result<AtomGroup> {
        let mut poo1 = AtomGroup::new();
        let p = ag
            .get_atom("P")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find P in FAD"))?;
        poo1.set_atom("P", p.clone());

        // amber format: OP1, pdb: O1P
        if let Some(o1p) = ag.get_atom("O1P") {
            poo1.set_atom("O1", o1p.clone());
        } else if let Some(op1) = ag.get_atom("OP1") {
            poo1.set_atom("O1", op1.clone());
        } else {
            return Err(BridgeError::input_error(
                "ag",
                "cannot find O1P or OP1 in FAD",
            ));
        }

        // amber format: OP2, pdb: O2P
        if let Some(o2p) = ag.get_atom("O2P") {
            poo1.set_atom("O2", o2p.clone());
        } else if let Some(op2) = ag.get_atom("OP2") {
            poo1.set_atom("O2", op2.clone());
        } else {
            return Err(BridgeError::input_error(
                "ag",
                "cannot find O2P or OP2 in FAD",
            ));
        }

        let na1_pos = self.get_neutralize_pos_poo_type(&poo1)?;
        let mut na1 = Atom::new_with_pos("Na", na1_pos)?;
        na1.name = "Na".to_string();

        let mut poo2 = AtomGroup::new();
        let pa = ag
            .get_atom("PA")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find PA in FAD"))?;
        let o1a = ag
            .get_atom("O1A")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find O1A in FAD"))?;
        let o2a = ag
            .get_atom("O2A")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find O2A in FAD"))?;
        poo2.set_atom("P", pa.clone());
        poo2.set_atom("O1", o1a.clone());
        poo2.set_atom("O2", o2a.clone());

        let na2_pos = self.get_neutralize_pos_poo_type(&poo2)?;
        let mut na2 = Atom::new_with_pos("Na", na2_pos)?;
        na2.name = "Na".to_string();

        let key = self.get_last_index(ag);
        let mut answer = AtomGroup::new();
        answer.set_atom(&format!("{}_Na1", key + 1), na1);
        answer.set_atom(&format!("{}_Na2", key + 1), na2);
        Ok(answer)
    }

    /// Computes the neutralization ion position for an NH3-type group.
    pub fn get_neutralize_pos_nh3_type(&self, ag: &AtomGroup) -> Result<Position> {
        let length = 3.187;
        let h1 = ag
            .get_atom("H1")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find H1"))?;
        let h2 = ag
            .get_atom("H2")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find H2"))?;
        let h3 = ag
            .get_atom("H3")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find H3"))?;
        let n = ag
            .get_atom("N")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find N"))?;

        let m = Position::new(
            (h1.xyz.x + h2.xyz.x + h3.xyz.x) / 3.0,
            (h1.xyz.y + h2.xyz.y + h3.xyz.y) / 3.0,
            (h1.xyz.z + h2.xyz.z + h3.xyz.z) / 3.0,
        );
        let mut v_nm = m - n.xyz;
        v_nm.norm()?;
        Ok(n.xyz + v_nm * length)
    }

    /// Computes the neutralization ion position for an NH2-type group (e.g. Proline N-term).
    pub fn get_neutralize_pos_nh2_type(&self, ag: &AtomGroup) -> Result<Position> {
        let length = 3.187;
        let h1 = ag
            .get_atom("H1")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find H1"))?;
        let h2 = ag
            .get_atom("H2")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find H2"))?;
        let n = ag
            .get_atom("N")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find N"))?;

        let v_nh1 = h1.xyz - n.xyz;
        let v_nh2 = h2.xyz - n.xyz;
        let mut v_m = (v_nh1 + v_nh2) * 0.5;
        v_m.norm()?;
        Ok(n.xyz + v_m * length)
    }

    /// Computes the neutralization ion position for a carboxylate (COO-) type group.
    pub fn get_neutralize_pos_coo_type(&self, ag: &AtomGroup) -> Result<Position> {
        let length = 2.521;
        let o1 = ag
            .get_atom("O1")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find O1"))?;
        let o2 = ag
            .get_atom("O2")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find O2"))?;
        let c = ag
            .get_atom("C")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find C"))?;

        let m = Position::new(
            0.5 * (o1.xyz.x + o2.xyz.x),
            0.5 * (o1.xyz.y + o2.xyz.y),
            0.5 * (o1.xyz.z + o2.xyz.z),
        );
        let mut v_cm = m - c.xyz;
        v_cm.norm()?;
        Ok(c.xyz + v_cm * length)
    }

    /// Computes the neutralization ion position for a phosphate (POO-) type group.
    pub fn get_neutralize_pos_poo_type(&self, ag: &AtomGroup) -> Result<Position> {
        let length = 2.748;
        let o1 = ag
            .get_atom("O1")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find O1"))?;
        let o2 = ag
            .get_atom("O2")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find O2"))?;
        let p = ag
            .get_atom("P")
            .ok_or_else(|| BridgeError::input_error("ag", "cannot find P"))?;

        let m = Position::new(
            0.5 * (o1.xyz.x + o2.xyz.x),
            0.5 * (o1.xyz.y + o2.xyz.y),
            0.5 * (o1.xyz.z + o2.xyz.z),
        );
        let mut v_pm = m - p.xyz;
        v_pm.norm()?;
        Ok(p.xyz + v_pm * length)
    }
}
