// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

//! Histidine protonation state determination (`HID`, `HIE`, `HIP`, or `MetalBridge`).
//!
//! Evaluates each histidine residue (`HIS`, `HID`, `HIE`, `HIP`) in an [`AtomGroup`]
//! without modifying the structure, following a strict priority order:
//! 0. Metal coordination (`ZN`, `FE`, `CU`, `MG`, `MN`, `CO`, `NI`; `Ca` excluded)
//! 1. Explicit residue name (`HID`, `HIE`, `HIP`)
//! 2. Existing ring hydrogens (`HD1`, `HE2`)
//! 3. Salt bridge to nearby `ASP` (`OD1`/`OD2`) or `GLU` (`OE1`/`OE2`) within 3.5 Å
//! 4. Heavy-atom burial count within `burial_radius` around `ND1` vs. `NE2`
//! 5. Default fallback (`Hid`)

use std::collections::BTreeMap;

use crate::atom::Atom;
use crate::atom_group::AtomGroup;
use crate::brd::sort_nicely;
use crate::position::Position;

/// Default maximum distance in Angstroms from a ring nitrogen (`ND1`/`NE2`) to a
/// coordinating metal atom.
pub const DEFAULT_METAL_DISTANCE: f64 = 2.5;

/// Default radius in Angstroms around a ring nitrogen (`ND1`/`NE2`) within which
/// surrounding heavy atoms of other residues are counted for burial comparison.
pub const DEFAULT_BURIAL_RADIUS: f64 = 8.0;

/// Maximum distance in Angstroms from a ring nitrogen (`ND1`/`NE2`) to an `ASP`
/// (`OD1`/`OD2`) or `GLU` (`OE1`/`OE2`) carboxylate oxygen to form a salt bridge.
pub const SALT_BRIDGE_DISTANCE: f64 = 3.5;

/// Element symbols of transition/divalent metals considered capable of coordinating
/// histidine ring nitrogens (`Ca` is intentionally excluded).
pub const COORDINATING_METALS: &[&str] = &["ZN", "FE", "CU", "MG", "MN", "CO", "NI"];

/// Options for histidine protonation state determination.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HisStateOptions {
    /// Maximum distance in Angstroms from `ND1` or `NE2` to a metal atom to consider
    /// the nitrogen coordinated (default: `2.5`).
    pub metal_distance: f64,
    /// Radius in Angstroms around `ND1` or `NE2` within which heavy atoms of other
    /// residues are counted for burial comparison (default: `8.0`).
    pub burial_radius: f64,
}

impl Default for HisStateOptions {
    fn default() -> Self {
        Self {
            metal_distance: DEFAULT_METAL_DISTANCE,
            burial_radius: DEFAULT_BURIAL_RADIUS,
        }
    }
}

/// Protonation or coordination state assigned to a histidine residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub enum HisState {
    /// Neutral histidine protonated at `ND1` (`HD1` present, `HE2` absent). Default state.
    #[default]
    Hid,
    /// Neutral histidine protonated at `NE2` (`HE2` present, `HD1` absent).
    Hie,
    /// Doubly protonated, positively charged histidine (`HD1` and `HE2` both present).
    Hip,
    /// Both `ND1` and `NE2` are coordinated to metal atoms within `metal_distance`;
    /// protonation state is left undecided.
    MetalBridge,
}

/// Rule that decided the [`HisState`] of a histidine site.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub enum HisDecidedBy {
    /// Decided by metal coordination within `metal_distance` of `ND1` and/or `NE2`.
    Metal,
    /// Decided by an explicit residue name (`HID`, `HIE`, or `HIP`).
    ResidueName,
    /// Decided by existing ring hydrogens (`HD1` and/or `HE2`) in the input residue.
    Hydrogens,
    /// Decided by a salt bridge (`<= 3.5` Å) to `ASP` (`OD1`/`OD2`) or `GLU` (`OE1`/`OE2`).
    SaltBridge,
    /// Decided by comparing heavy-atom burial counts within `burial_radius` around `ND1` and `NE2`.
    Burial,
    /// Assigned the default state ([`HisState::Hid`]) when no earlier rule decided.
    #[default]
    Default,
}

/// Per-residue report of histidine state determination.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct HisSiteReport {
    /// Hierarchical path of the residue group (e.g. `"/model_1/A/1/"`).
    pub path: String,
    /// Chemical component identifier (residue group's `name`, e.g. `"HIS"`).
    pub comp_id: String,
    /// Assigned histidine state.
    pub state: HisState,
    /// Rule that decided `state`.
    pub decided_by: HisDecidedBy,
}

/// Structure-wide report of histidine state determination across residues.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct HisReport {
    /// Per-site results for all histidine residues in traversal order.
    pub sites: Vec<HisSiteReport>,
}

impl HisReport {
    /// Returns a map from [`HisState`] to the number of histidine sites in that state.
    pub fn count_by_state(&self) -> BTreeMap<HisState, usize> {
        let mut counts = BTreeMap::new();
        for site in &self.sites {
            *counts.entry(site.state).or_insert(0) += 1;
        }
        counts
    }
}

/// Returns `true` if `atom` is hydrogen (`H` or `D`).
fn is_hydrogen(atom: &Atom) -> bool {
    atom.atomic_number() == 1
        || atom
            .symbol()
            .is_ok_and(|s| s.trim().eq_ignore_ascii_case("H") || s.trim().eq_ignore_ascii_case("D"))
}

/// Returns `true` if `atom` is one of the coordinating metal elements (`ZN`, `FE`, `CU`, `MG`, `MN`, `CO`, `NI`).
fn is_coordinating_metal(atom: &Atom) -> bool {
    let Ok(sym) = atom.symbol() else {
        return false;
    };
    let trimmed = sym.trim();
    COORDINATING_METALS
        .iter()
        .any(|&m| m.eq_ignore_ascii_case(trimmed))
}

/// Finds an atom in `residue` by trimmed `Atom.name`.
fn find_atom_by_name<'a>(residue: &'a AtomGroup, name: &str) -> Option<&'a Atom> {
    residue
        .atoms()
        .map(|(_, a)| a)
        .find(|a| a.name.trim() == name)
}

/// Checks if `residue` contains an atom with the given trimmed `Atom.name`.
fn has_atom_by_name(residue: &AtomGroup, name: &str) -> bool {
    find_atom_by_name(residue, name).is_some()
}

/// Returns `true` if `comp_id` is a histidine residue name (`HIS`, `HID`, `HIE`, `HIP`).
fn is_histidine_comp_id(comp_id: &str) -> bool {
    let trimmed = comp_id.trim();
    trimmed.eq_ignore_ascii_case("HIS")
        || trimmed.eq_ignore_ascii_case("HID")
        || trimmed.eq_ignore_ascii_case("HIE")
        || trimmed.eq_ignore_ascii_case("HIP")
}

/// Returns `true` if an atom with `atom_name` in a residue with `comp_id` is an
/// acidic carboxylate oxygen capable of forming a salt bridge (`ASP` `OD1`/`OD2` or `GLU` `OE1`/`OE2`).
fn is_acidic_carboxylate_oxygen(comp_id: &str, atom_name: &str) -> bool {
    let comp = comp_id.trim();
    let name = atom_name.trim();
    if comp.eq_ignore_ascii_case("ASP") {
        name == "OD1" || name == "OD2"
    } else if comp.eq_ignore_ascii_case("GLU") {
        name == "OE1" || name == "OE2"
    } else {
        false
    }
}

struct ScopedResidue<'a> {
    residue: &'a AtomGroup,
    scope_id: usize,
}

#[derive(Default)]
struct ScopeTraversal<'a> {
    residues: Vec<ScopedResidue<'a>>,
    num_scopes: usize,
}

impl<'a> ScopeTraversal<'a> {
    fn new_scope(&mut self) -> usize {
        let id = self.num_scopes;
        self.num_scopes += 1;
        id
    }

    fn push_residue(&mut self, residue: &'a AtomGroup, scope_id: usize) {
        self.residues.push(ScopedResidue { residue, scope_id });
    }
}

fn collect_from_chain<'a>(chain: &'a AtomGroup, scope_id: usize, state: &mut ScopeTraversal<'a>) {
    let mut res_keys = chain.get_group_list();
    sort_nicely(&mut res_keys);

    for key in &res_keys {
        if let Some(residue) = chain.get_group(key) {
            if residue.get_number_of_atoms() > 0 {
                state.push_residue(residue, scope_id);
            }
        }
    }
}

/// Traverses an [`AtomGroup`] hierarchy using the exact same residue discovery order
/// and per-model scope isolation rules as `missing_atoms::collect_from_group`.
fn collect_from_group<'a>(
    group: &'a AtomGroup,
    inherited_scope: Option<usize>,
    state: &mut ScopeTraversal<'a>,
) {
    if group.get_number_of_groups() == 0 {
        if group.get_number_of_atoms() > 0 {
            let scope = inherited_scope.unwrap_or_else(|| state.new_scope());
            state.push_residue(group, scope);
        }
        return;
    }

    let all_children_are_leaves = group
        .groups()
        .all(|(_, child)| child.get_number_of_groups() == 0);

    if all_children_are_leaves {
        let scope = inherited_scope.unwrap_or_else(|| state.new_scope());
        collect_from_chain(group, scope, state);
    } else {
        let mut child_keys = group.get_group_list();
        sort_nicely(&mut child_keys);

        let mut shared_chain_scope = inherited_scope;

        for key in &child_keys {
            if let Some(child) = group.get_group(key) {
                if child.get_number_of_groups() == 0 {
                    if child.get_number_of_atoms() > 0 {
                        let scope = *shared_chain_scope.get_or_insert_with(|| state.new_scope());
                        state.push_residue(child, scope);
                    }
                } else if child.groups().all(|(_, gc)| gc.get_number_of_groups() == 0) {
                    let scope = *shared_chain_scope.get_or_insert_with(|| state.new_scope());
                    collect_from_chain(child, scope, state);
                } else {
                    collect_from_group(child, inherited_scope, state);
                }
            }
        }
    }
}

#[derive(Default)]
struct ScopeEnvironment {
    metals: Vec<Position>,
    carboxyl_oxygens: Vec<(usize, Position)>,
    heavy_atoms: Vec<(usize, Position)>,
}

fn decide_single_his(
    residue: &AtomGroup,
    res_idx: usize,
    env: &ScopeEnvironment,
    options: &HisStateOptions,
) -> HisSiteReport {
    let path = residue.path().to_string();
    let comp_id = residue.name.clone();
    let comp_trimmed = comp_id.trim();

    let nd1_pos = find_atom_by_name(residue, "ND1").map(|a| a.xyz);
    let ne2_pos = find_atom_by_name(residue, "NE2").map(|a| a.xyz);

    // 0. Metal coordination
    let nd1_metal = nd1_pos.is_some_and(|p| {
        env.metals
            .iter()
            .any(|m| p.distance_from(m) <= options.metal_distance)
    });
    let ne2_metal = ne2_pos.is_some_and(|p| {
        env.metals
            .iter()
            .any(|m| p.distance_from(m) <= options.metal_distance)
    });

    if nd1_metal && ne2_metal {
        return HisSiteReport {
            path,
            comp_id,
            state: HisState::MetalBridge,
            decided_by: HisDecidedBy::Metal,
        };
    }
    if nd1_metal {
        return HisSiteReport {
            path,
            comp_id,
            state: HisState::Hie,
            decided_by: HisDecidedBy::Metal,
        };
    }
    if ne2_metal {
        return HisSiteReport {
            path,
            comp_id,
            state: HisState::Hid,
            decided_by: HisDecidedBy::Metal,
        };
    }

    // 1. Explicit residue name (HID, HIE, HIP)
    if comp_trimmed.eq_ignore_ascii_case("HID") {
        return HisSiteReport {
            path,
            comp_id,
            state: HisState::Hid,
            decided_by: HisDecidedBy::ResidueName,
        };
    }
    if comp_trimmed.eq_ignore_ascii_case("HIE") {
        return HisSiteReport {
            path,
            comp_id,
            state: HisState::Hie,
            decided_by: HisDecidedBy::ResidueName,
        };
    }
    if comp_trimmed.eq_ignore_ascii_case("HIP") {
        return HisSiteReport {
            path,
            comp_id,
            state: HisState::Hip,
            decided_by: HisDecidedBy::ResidueName,
        };
    }

    // 2. Input hydrogens (HD1, HE2)
    let has_hd1 = has_atom_by_name(residue, "HD1");
    let has_he2 = has_atom_by_name(residue, "HE2");
    if has_hd1 && has_he2 {
        return HisSiteReport {
            path,
            comp_id,
            state: HisState::Hip,
            decided_by: HisDecidedBy::Hydrogens,
        };
    }
    if has_hd1 {
        return HisSiteReport {
            path,
            comp_id,
            state: HisState::Hid,
            decided_by: HisDecidedBy::Hydrogens,
        };
    }
    if has_he2 {
        return HisSiteReport {
            path,
            comp_id,
            state: HisState::Hie,
            decided_by: HisDecidedBy::Hydrogens,
        };
    }

    // 3. Salt bridge to other residues' ASP (OD1/OD2) or GLU (OE1/OE2) within 3.5 Å
    let has_salt_bridge = [nd1_pos, ne2_pos].into_iter().flatten().any(|n_pos| {
        env.carboxyl_oxygens.iter().any(|(other_idx, o_pos)| {
            *other_idx != res_idx && n_pos.distance_from(o_pos) <= SALT_BRIDGE_DISTANCE
        })
    });
    if has_salt_bridge {
        return HisSiteReport {
            path,
            comp_id,
            state: HisState::Hip,
            decided_by: HisDecidedBy::SaltBridge,
        };
    }

    // 4. Heavy-atom burial within burial_radius (other residues only)
    if let (Some(nd1), Some(ne2)) = (nd1_pos, ne2_pos) {
        let mut nd1_burial = 0usize;
        let mut ne2_burial = 0usize;
        for &(other_idx, pos) in &env.heavy_atoms {
            if other_idx == res_idx {
                continue;
            }
            if nd1.distance_from(&pos) <= options.burial_radius {
                nd1_burial += 1;
            }
            if ne2.distance_from(&pos) <= options.burial_radius {
                ne2_burial += 1;
            }
        }
        if nd1_burial > ne2_burial {
            return HisSiteReport {
                path,
                comp_id,
                state: HisState::Hid,
                decided_by: HisDecidedBy::Burial,
            };
        }
        if ne2_burial > nd1_burial {
            return HisSiteReport {
                path,
                comp_id,
                state: HisState::Hie,
                decided_by: HisDecidedBy::Burial,
            };
        }
    }

    // 5. Default fallback (HID)
    HisSiteReport {
        path,
        comp_id,
        state: HisState::Hid,
        decided_by: HisDecidedBy::Default,
    }
}

/// Determines the protonation or coordination state of all histidine residues
/// (`HIS`, `HID`, `HIE`, `HIP`) in `group` without modifying the tree.
pub fn his_states(group: &AtomGroup, options: &HisStateOptions) -> HisReport {
    let mut traversal = ScopeTraversal::default();
    collect_from_group(group, None, &mut traversal);

    let mut envs: Vec<ScopeEnvironment> = (0..traversal.num_scopes)
        .map(|_| ScopeEnvironment::default())
        .collect();

    for (res_idx, item) in traversal.residues.iter().enumerate() {
        let env = &mut envs[item.scope_id];
        let comp_id = item.residue.name.as_str();
        for (_, atom) in item.residue.atoms() {
            if is_hydrogen(atom) {
                continue;
            }
            env.heavy_atoms.push((res_idx, atom.xyz));
            if is_coordinating_metal(atom) {
                env.metals.push(atom.xyz);
            }
            if is_acidic_carboxylate_oxygen(comp_id, &atom.name) {
                env.carboxyl_oxygens.push((res_idx, atom.xyz));
            }
        }
    }

    let mut sites = Vec::new();
    for (res_idx, item) in traversal.residues.iter().enumerate() {
        if is_histidine_comp_id(&item.residue.name) {
            let site = decide_single_his(item.residue, res_idx, &envs[item.scope_id], options);
            sites.push(site);
        }
    }

    HisReport { sites }
}

impl AtomGroup {
    /// Determines the protonation or coordination state of each histidine residue
    /// (`HIS`, `HID`, `HIE`, `HIP`) in this hierarchy without modifying the tree.
    pub fn his_states(&self, options: &HisStateOptions) -> HisReport {
        his_states(self, options)
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;
    use crate::format::pdb::Pdb;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../../proteindf_bridge/data")
            .canonicalize()
            .expect("proteindf_bridge/data directory must exist")
    }

    fn make_atom(symbol: &str, name: &str, xyz: Position) -> Atom {
        let mut atom = Atom::new_with_pos(symbol, xyz).expect("valid element symbol");
        atom.name = name.to_string();
        atom
    }

    /// Builds a minimal histidine residue with `ND1` at `(0, 0, 0)` and `NE2` at `(20, 0, 0)`.
    fn make_his_residue(comp_id: &str) -> AtomGroup {
        let mut res = AtomGroup::with_name(comp_id);
        res.set_atom("ND1", make_atom("N", "ND1", Position::new(0.0, 0.0, 0.0)));
        res.set_atom("NE2", make_atom("N", "NE2", Position::new(20.0, 0.0, 0.0)));
        res
    }

    fn make_chain_with_two_residues(res1: AtomGroup, res2: AtomGroup) -> AtomGroup {
        let mut chain = AtomGroup::with_name("A");
        chain.set_group("1", res1);
        chain.set_group("2", res2);
        chain
    }

    // ① ND1 の近く(2.2 Å)に ZN がある → HIE(Metal)
    #[test]
    fn test_01_metal_near_nd1_gives_hie() {
        let opts = HisStateOptions::default();

        for &metal_sym in COORDINATING_METALS {
            let his = make_his_residue("HIS");
            let mut metal_res = AtomGroup::with_name(metal_sym);
            metal_res.set_atom(
                metal_sym,
                make_atom(metal_sym, metal_sym, Position::new(2.2, 0.0, 0.0)),
            );

            let chain = make_chain_with_two_residues(his, metal_res);
            let report = chain.his_states(&opts);
            assert_eq!(report.sites.len(), 1);
            assert_eq!(report.sites[0].state, HisState::Hie);
            assert_eq!(report.sites[0].decided_by, HisDecidedBy::Metal);
        }
    }

    // ② NE2 の近く → HID(Metal)
    #[test]
    fn test_02_metal_near_ne2_gives_hid() {
        let opts = HisStateOptions::default();
        let his = make_his_residue("HIS");
        let mut zn_res = AtomGroup::with_name("ZN");
        // 2.2 Å from NE2 at (20, 0, 0)
        zn_res.set_atom("ZN", make_atom("Zn", "ZN", Position::new(22.2, 0.0, 0.0)));

        let chain = make_chain_with_two_residues(his, zn_res);
        let report = chain.his_states(&opts);
        assert_eq!(report.sites.len(), 1);
        assert_eq!(report.sites[0].state, HisState::Hid);
        assert_eq!(report.sites[0].decided_by, HisDecidedBy::Metal);
    }

    // ③ 両方の窒素が 2.5 Å 以内に金属 → MetalBridge
    #[test]
    fn test_03_metal_near_both_nitrogens_gives_metal_bridge() {
        let opts = HisStateOptions::default();
        let his = make_his_residue("HIS");
        let mut zn1 = AtomGroup::with_name("ZN");
        zn1.set_atom("ZN", make_atom("Zn", "ZN", Position::new(2.5, 0.0, 0.0)));
        let mut zn2 = AtomGroup::with_name("CU");
        zn2.set_atom("CU", make_atom("Cu", "CU", Position::new(22.4, 0.0, 0.0)));

        let mut chain = AtomGroup::with_name("A");
        chain.set_group("1", his);
        chain.set_group("2", zn1);
        chain.set_group("3", zn2);

        let report = chain.his_states(&opts);
        assert_eq!(report.sites.len(), 1);
        assert_eq!(report.sites[0].state, HisState::MetalBridge);
        assert_eq!(report.sites[0].decided_by, HisDecidedBy::Metal);
        assert_eq!(
            report.count_by_state(),
            BTreeMap::from([(HisState::MetalBridge, 1)])
        );
    }

    // ④ 金属が 3.0 Å なら配位とみなさない
    #[test]
    fn test_04_metal_at_3_0_angstrom_not_coordinated() {
        let opts = HisStateOptions::default();
        let his = make_his_residue("HIS");
        let mut zn_res = AtomGroup::with_name("ZN");
        // 3.0 Å from ND1 (> 2.5 Å metal_distance, but <= 8.0 Å burial_radius)
        zn_res.set_atom("ZN", make_atom("Zn", "ZN", Position::new(3.0, 0.0, 0.0)));

        let chain = make_chain_with_two_residues(his, zn_res);
        let report = chain.his_states(&opts);
        assert_eq!(report.sites.len(), 1);
        assert_ne!(report.sites[0].decided_by, HisDecidedBy::Metal);
        // Since ZN is a heavy atom within 8.0 Å of ND1, burial decides HID
        assert_eq!(report.sites[0].state, HisState::Hid);
        assert_eq!(report.sites[0].decided_by, HisDecidedBy::Burial);
    }

    // ⑤ 残基名 HID → HID(ResidueName) (HIE, HIP も確認)
    #[test]
    fn test_05_residue_name_hid_hie_hip() {
        let opts = HisStateOptions::default();

        let mut chain = AtomGroup::with_name("A");
        chain.set_group("1", make_his_residue("HID"));
        chain.set_group("2", make_his_residue("HIE"));
        chain.set_group("3", make_his_residue("HIP"));

        let report = chain.his_states(&opts);
        assert_eq!(report.sites.len(), 3);
        assert_eq!(report.sites[0].comp_id, "HID");
        assert_eq!(report.sites[0].state, HisState::Hid);
        assert_eq!(report.sites[0].decided_by, HisDecidedBy::ResidueName);

        assert_eq!(report.sites[1].comp_id, "HIE");
        assert_eq!(report.sites[1].state, HisState::Hie);
        assert_eq!(report.sites[1].decided_by, HisDecidedBy::ResidueName);

        assert_eq!(report.sites[2].comp_id, "HIP");
        assert_eq!(report.sites[2].state, HisState::Hip);
        assert_eq!(report.sites[2].decided_by, HisDecidedBy::ResidueName);

        assert_eq!(
            report.count_by_state(),
            BTreeMap::from([(HisState::Hid, 1), (HisState::Hie, 1), (HisState::Hip, 1),])
        );
    }

    // ⑥ 水素 HD1 のみ → HID(Hydrogens) (HE2 のみ・両方ありも確認)
    #[test]
    fn test_06_input_hydrogens_hd1_he2() {
        let opts = HisStateOptions::default();

        let mut his_hd1 = make_his_residue("HIS");
        his_hd1.set_atom("HD1", make_atom("H", "HD1", Position::new(-1.0, 0.0, 0.0)));

        let mut his_he2 = make_his_residue("HIS");
        his_he2.set_atom("HE2", make_atom("H", "HE2", Position::new(21.0, 0.0, 0.0)));

        let mut his_both = make_his_residue("HIS");
        his_both.set_atom("HD1", make_atom("H", "HD1", Position::new(-1.0, 0.0, 0.0)));
        his_both.set_atom("HE2", make_atom("H", "HE2", Position::new(21.0, 0.0, 0.0)));

        let rep_hd1 = his_hd1.his_states(&opts);
        assert_eq!(rep_hd1.sites.len(), 1);
        assert_eq!(rep_hd1.sites[0].state, HisState::Hid);
        assert_eq!(rep_hd1.sites[0].decided_by, HisDecidedBy::Hydrogens);

        let rep_he2 = his_he2.his_states(&opts);
        assert_eq!(rep_he2.sites.len(), 1);
        assert_eq!(rep_he2.sites[0].state, HisState::Hie);
        assert_eq!(rep_he2.sites[0].decided_by, HisDecidedBy::Hydrogens);

        let rep_both = his_both.his_states(&opts);
        assert_eq!(rep_both.sites.len(), 1);
        assert_eq!(rep_both.sites[0].state, HisState::Hip);
        assert_eq!(rep_both.sites[0].decided_by, HisDecidedBy::Hydrogens);
    }

    // ⑦ 水素がなく、密度が等しい → HID(Default)
    #[test]
    fn test_07_no_hydrogens_equal_density_defaults_to_hid() {
        let opts = HisStateOptions::default();

        // Isolated HIS (0 heavy atoms around both ND1 and NE2)
        let his = make_his_residue("HIS");
        let rep_isolated = his.his_states(&opts);
        assert_eq!(rep_isolated.sites.len(), 1);
        assert_eq!(rep_isolated.sites[0].state, HisState::Hid);
        assert_eq!(rep_isolated.sites[0].decided_by, HisDecidedBy::Default);

        // Symmetric heavy atoms (1 near ND1 at 5.0 Å, 1 near NE2 at 5.0 Å)
        let mut ala = AtomGroup::with_name("ALA");
        ala.set_atom("CA", make_atom("C", "CA", Position::new(5.0, 0.0, 0.0)));
        ala.set_atom("CB", make_atom("C", "CB", Position::new(25.0, 0.0, 0.0)));
        let chain = make_chain_with_two_residues(make_his_residue("HIS"), ala);
        let rep_sym = chain.his_states(&opts);
        assert_eq!(rep_sym.sites.len(), 1);
        assert_eq!(rep_sym.sites[0].state, HisState::Hid);
        assert_eq!(rep_sym.sites[0].decided_by, HisDecidedBy::Default);
    }

    // ⑧ 密度が高い方 ND1 → HID(Burial) (NE2 が高いとき HIE(Burial) も確認)
    #[test]
    fn test_08_burial_higher_density_nd1_and_ne2() {
        let opts = HisStateOptions::default();

        // Higher density near ND1 (at origin): two heavy atoms within 8.0 Å of ND1, one near NE2
        let mut env_nd1 = AtomGroup::with_name("ALA");
        env_nd1.set_atom("N", make_atom("N", "N", Position::new(4.0, 0.0, 0.0)));
        env_nd1.set_atom("CA", make_atom("C", "CA", Position::new(5.0, 0.0, 0.0)));
        env_nd1.set_atom("CB", make_atom("C", "CB", Position::new(25.0, 0.0, 0.0)));
        // Hydrogen near NE2 must NOT be counted as heavy atom
        env_nd1.set_atom("HB1", make_atom("H", "HB1", Position::new(24.0, 0.0, 0.0)));
        env_nd1.set_atom("HB2", make_atom("H", "HB2", Position::new(23.0, 0.0, 0.0)));

        let chain_nd1 = make_chain_with_two_residues(make_his_residue("HIS"), env_nd1);
        let rep_nd1 = chain_nd1.his_states(&opts);
        assert_eq!(rep_nd1.sites.len(), 1);
        assert_eq!(rep_nd1.sites[0].state, HisState::Hid);
        assert_eq!(rep_nd1.sites[0].decided_by, HisDecidedBy::Burial);

        // Higher density near NE2 (at (20, 0, 0))
        let mut env_ne2 = AtomGroup::with_name("ALA");
        env_ne2.set_atom("N", make_atom("N", "N", Position::new(24.0, 0.0, 0.0)));
        env_ne2.set_atom("CA", make_atom("C", "CA", Position::new(25.0, 0.0, 0.0)));
        let chain_ne2 = make_chain_with_two_residues(make_his_residue("HIS"), env_ne2);
        let rep_ne2 = chain_ne2.his_states(&opts);
        assert_eq!(rep_ne2.sites.len(), 1);
        assert_eq!(rep_ne2.sites[0].state, HisState::Hie);
        assert_eq!(rep_ne2.sites[0].decided_by, HisDecidedBy::Burial);
    }

    // ⑨a 窒素から 3.5 Å に Asp の OD1 がある → HIP(SaltBridge)
    // ⑨b 3.5 Å を超える → 塩橋とみなさない
    // ⑨c Arg の近くは塩橋とみなさない
    #[test]
    fn test_09_salt_bridge_asp_glu_distance_and_arg_excluded() {
        let opts = HisStateOptions::default();

        // ⑨a: ASP OD1 at exactly 3.5 Å from ND1 -> HIP (SaltBridge)
        let mut asp_within = AtomGroup::with_name("ASP");
        asp_within.set_atom("OD1", make_atom("O", "OD1", Position::new(3.5, 0.0, 0.0)));
        let chain_9a = make_chain_with_two_residues(make_his_residue("HIS"), asp_within);
        let rep_9a = chain_9a.his_states(&opts);
        assert_eq!(rep_9a.sites.len(), 1);
        assert_eq!(rep_9a.sites[0].state, HisState::Hip);
        assert_eq!(rep_9a.sites[0].decided_by, HisDecidedBy::SaltBridge);

        // Also check GLU OE2 within 3.5 Å from NE2 -> HIP (SaltBridge)
        let mut glu_within = AtomGroup::with_name("GLU");
        glu_within.set_atom("OE2", make_atom("O", "OE2", Position::new(23.2, 0.0, 0.0)));
        let chain_glu = make_chain_with_two_residues(make_his_residue("HIS"), glu_within);
        let rep_glu = chain_glu.his_states(&opts);
        assert_eq!(rep_glu.sites.len(), 1);
        assert_eq!(rep_glu.sites[0].state, HisState::Hip);
        assert_eq!(rep_glu.sites[0].decided_by, HisDecidedBy::SaltBridge);

        // ⑨b: ASP OD1 beyond 3.5 Å (3.6 Å) -> NOT SaltBridge
        let mut asp_beyond = AtomGroup::with_name("ASP");
        asp_beyond.set_atom("OD1", make_atom("O", "OD1", Position::new(3.6, 0.0, 0.0)));
        let chain_9b = make_chain_with_two_residues(make_his_residue("HIS"), asp_beyond);
        let rep_9b = chain_9b.his_states(&opts);
        assert_eq!(rep_9b.sites.len(), 1);
        assert_ne!(rep_9b.sites[0].decided_by, HisDecidedBy::SaltBridge);
        assert_eq!(rep_9b.sites[0].state, HisState::Hid);
        assert_eq!(rep_9b.sites[0].decided_by, HisDecidedBy::Burial);

        // ⑨c: ARG atoms within 3.0 Å -> NOT SaltBridge
        let mut arg_near = AtomGroup::with_name("ARG");
        arg_near.set_atom("NH1", make_atom("N", "NH1", Position::new(3.0, 0.0, 0.0)));
        arg_near.set_atom("NH2", make_atom("N", "NH2", Position::new(3.2, 0.0, 0.0)));
        let chain_9c = make_chain_with_two_residues(make_his_residue("HIS"), arg_near);
        let rep_9c = chain_9c.his_states(&opts);
        assert_eq!(rep_9c.sites.len(), 1);
        assert_ne!(rep_9c.sites[0].decided_by, HisDecidedBy::SaltBridge);
        assert_eq!(rep_9c.sites[0].state, HisState::Hid);
        assert_eq!(rep_9c.sites[0].decided_by, HisDecidedBy::Burial);
    }

    // ⑩ Ca は配位とみなさない
    #[test]
    fn test_10_calcium_not_coordinating_metal() {
        let opts = HisStateOptions::default();
        let his = make_his_residue("HIS");
        let mut ca_res = AtomGroup::with_name("CA");
        // Calcium ion 2.2 Å from ND1 (would give HIE if treated as coordinating metal)
        ca_res.set_atom("CA", make_atom("Ca", "CA", Position::new(2.2, 0.0, 0.0)));

        let chain = make_chain_with_two_residues(his, ca_res);
        let report = chain.his_states(&opts);
        assert_eq!(report.sites.len(), 1);
        assert_ne!(report.sites[0].decided_by, HisDecidedBy::Metal);
        // Ca is counted only as a surrounding heavy atom near ND1 -> HID (Burial)
        assert_eq!(report.sites[0].state, HisState::Hid);
        assert_eq!(report.sites[0].decided_by, HisDecidedBy::Burial);
    }

    // ⑪ 空の木と、残基・鎖・モデルから呼んで同じ結果（モデルをまたがないことも確認）
    #[test]
    fn test_11_empty_tree_hierarchy_levels_and_model_isolation() {
        let opts = HisStateOptions::default();

        // Empty tree
        let empty = AtomGroup::new();
        let rep_empty = empty.his_states(&opts);
        assert!(rep_empty.sites.is_empty());
        assert!(rep_empty.count_by_state().is_empty());

        // Single residue in root -> model_1 -> A -> 1
        let mut his = make_his_residue("HIS");
        his.set_atom("HE2", make_atom("H", "HE2", Position::new(21.0, 0.0, 0.0)));

        let mut chain = AtomGroup::with_name("A");
        chain.set_group("1", his);
        let mut model = AtomGroup::with_name("model_1");
        model.set_group("A", chain);
        let mut root = AtomGroup::new();
        root.set_group("model_1", model);

        let rep_root = root.his_states(&opts);
        let rep_model = root.get_group("model_1").unwrap().his_states(&opts);
        let rep_chain = root
            .get_group("model_1")
            .unwrap()
            .get_group("A")
            .unwrap()
            .his_states(&opts);
        let rep_res = root
            .get_group("model_1")
            .unwrap()
            .get_group("A")
            .unwrap()
            .get_group("1")
            .unwrap()
            .his_states(&opts);

        assert_eq!(rep_root.sites.len(), 1);
        assert_eq!(rep_root.sites[0].path, "/model_1/A/1/");
        assert_eq!(rep_root.sites[0].state, HisState::Hie);
        assert_eq!(rep_root.sites[0].decided_by, HisDecidedBy::Hydrogens);
        assert_eq!(rep_root, rep_model);
        assert_eq!(rep_root, rep_chain);
        assert_eq!(rep_root, rep_res);

        // Cross-model isolation: ZN in model_2 at 2.0 Å from model_1's ND1 must NOT affect model_1
        let mut m1_chain = AtomGroup::with_name("A");
        m1_chain.set_group("1", make_his_residue("HIS"));
        let mut m1 = AtomGroup::with_name("model_1");
        m1.set_group("A", m1_chain);

        let mut zn_res = AtomGroup::with_name("ZN");
        zn_res.set_atom("ZN", make_atom("Zn", "ZN", Position::new(2.0, 0.0, 0.0)));
        let mut m2_chain = AtomGroup::with_name("A");
        m2_chain.set_group("1", zn_res);
        let mut m2 = AtomGroup::with_name("model_2");
        m2.set_group("A", m2_chain);

        let mut multi_root = AtomGroup::new();
        multi_root.set_group("model_1", m1);
        multi_root.set_group("model_2", m2);

        let rep_multi = multi_root.his_states(&opts);
        assert_eq!(rep_multi.sites.len(), 1);
        assert_eq!(rep_multi.sites[0].path, "/model_1/A/1/");
        assert_eq!(rep_multi.sites[0].state, HisState::Hid);
        assert_eq!(rep_multi.sites[0].decided_by, HisDecidedBy::Default);
    }

    // ⑫ 木が変わらない、および 1hls の 3 残基の確認
    #[test]
    fn test_12_tree_unchanged_and_1hls_three_residues() {
        let opts = HisStateOptions::default();

        let pdb_path = test_data_dir().join("1hls.pdb");
        let pdb = Pdb::from_file(&pdb_path, None).expect("load 1hls.pdb");
        let ag_1hls = pdb.get_atomgroup(None, None).expect("1hls atomgroup");
        let snapshot = ag_1hls.clone();

        let report = ag_1hls.his_states(&opts);
        assert_eq!(
            ag_1hls, snapshot,
            "tree must not be modified by his_states()"
        );
        assert_eq!(report.sites.len(), 3);

        for site in &report.sites {
            println!(
                "1hls site: path={}, comp_id={}, state={:?}, decided_by={:?}",
                site.path, site.comp_id, site.state, site.decided_by
            );
            assert_eq!(site.comp_id, "HIS");
            assert_eq!(site.state, HisState::Hid);
            assert_eq!(site.decided_by, HisDecidedBy::Hydrogens);
        }
        assert_eq!(
            report.count_by_state(),
            BTreeMap::from([(HisState::Hid, 3)])
        );

        // Also check 1hls when hydrogens are stripped (heavy atoms only)
        fn strip_hydrogens(group: &mut AtomGroup) {
            let h_keys: Vec<String> = group
                .atoms()
                .filter(|(_, a)| is_hydrogen(a))
                .map(|(k, _)| k.to_string())
                .collect();
            for k in h_keys {
                group.remove_atom(&k);
            }
            let g_keys = group.get_group_list();
            for gk in g_keys {
                if let Some(child) = group.get_group_mut(&gk) {
                    strip_hydrogens(child);
                }
            }
        }

        let mut ag_1hls_no_h = ag_1hls.clone();
        strip_hydrogens(&mut ag_1hls_no_h);
        let report_no_h = ag_1hls_no_h.his_states(&opts);
        assert_eq!(report_no_h.sites.len(), 3);
        for site in &report_no_h.sites {
            println!(
                "1hls (no H) site: path={}, comp_id={}, state={:?}, decided_by={:?}",
                site.path, site.comp_id, site.state, site.decided_by
            );
            assert_eq!(site.comp_id, "HIS");
            assert_eq!(site.state, HisState::Hid);
            assert_eq!(site.decided_by, HisDecidedBy::Burial);
        }
    }
}
