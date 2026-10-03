# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import tempfile
import unittest

import proteindf_bridge_rs as rs_br


class TestRsAnalysis(unittest.TestCase):
    """
    Test suite for Python bindings of Analysis features (Phase 8 & 9):
    1. Backbone & Sidechain Hydrogen Bonds
    2. Secondary Structure (DSSP) calculation and in-place application
    3. CH-pi Interactions with tunable thresholds
    4. InteractionSet aggregation, filtering, and MessagePack/YAML roundtrip
    """

    def setUp(self):
        self.data_dir = os.path.join(
            os.path.dirname(__file__), "..", "rust", "crates", "proteindf-bridge", "tests", "data"
        )
        self.pdb_1hls = os.path.join(self.data_dir, "1hls.pdb")

    def _load_1hls(self):
        pdb = rs_br.Pdb()
        pdb.load(self.pdb_1hls)
        return pdb.get_atomgroup()

    # -------------------------------------------------------------------------
    # 1. Backbone Hydrogen Bonds
    # -------------------------------------------------------------------------
    def test_backbone_hbonds_1hls(self):
        """
        Verifies backbone hydrogen bonds on 1hls chain A.
        Baseline: rust/crates/proteindf-bridge/tests/test_hydrogen_bond.rs line 151 (len=9)
        and lines 154-164 for expected residue pairs and Kabsch-Sander energies.
        """
        ag = self._load_1hls()
        chain_a = ag["model_1"]["A"]

        hbonds = rs_br.calc_backbone_hbonds(chain_a)

        # Baseline: test_hydrogen_bond.rs line 151
        self.assertEqual(len(hbonds), 9)

        # Expected 9 hbonds from test_hydrogen_bond.rs lines 154-164:
        # (donor, acceptor, expected_energy)
        expected = [
            ("6", "2", -1.6882),
            ("7", "3", -0.6522),
            ("8", "3", -1.6268),
            ("15", "12", -0.5502),
            ("16", "12", -0.5110),
            ("17", "14", -0.7758),
            ("18", "14", -0.5546),
            ("19", "16", -0.8396),
            ("20", "17", -1.4639),
        ]

        for exp_d, exp_a, exp_e in expected:
            found = [
                hb
                for hb in hbonds
                if hb.donor_residue_key == exp_d and hb.acceptor_residue_key == exp_a
            ]
            self.assertEqual(
                len(found),
                1,
                f"Expected backbone hbond ({exp_d} -> {exp_a}) not found in results",
            )
            hb = found[0]
            self.assertAlmostEqual(
                hb.energy,
                exp_e,
                places=3,
                msg=f"Energy mismatch for ({exp_d} -> {exp_a}): expected {exp_e}, got {hb.energy}",
            )
            self.assertIn("HydrogenBond(", repr(hb))

        # Check residue 1 is not donor (test_hydrogen_bond.rs line 93-96)
        self.assertFalse(any(hb.donor_residue_key == "1" for hb in hbonds))

        # Check |d - a| > 2 exclusion (test_hydrogen_bond.rs lines 124-138)
        for hb in hbonds:
            d_idx = int(hb.donor_residue_key)
            a_idx = int(hb.acceptor_residue_key)
            self.assertGreater(abs(d_idx - a_idx), 2)
            self.assertLess(hb.energy, -0.5)

    def test_backbone_hbonds_determinism(self):
        """
        Verifies that repeated calls to calc_backbone_hbonds return deterministically ordered results.
        """
        ag = self._load_1hls()
        chain_a = ag["model_1"]["A"]

        run1 = rs_br.calc_backbone_hbonds(chain_a)
        run2 = rs_br.calc_backbone_hbonds(chain_a)

        self.assertEqual(len(run1), len(run2))
        for hb1, hb2 in zip(run1, run2):
            self.assertEqual(hb1.donor_residue_key, hb2.donor_residue_key)
            self.assertEqual(hb1.acceptor_residue_key, hb2.acceptor_residue_key)
            self.assertEqual(hb1.energy, hb2.energy)

    # -------------------------------------------------------------------------
    # 2. Sidechain Hydrogen Bonds
    # -------------------------------------------------------------------------
    def test_sidechain_hbonds_1hls_benchmark(self):
        """
        Verifies sidechain hydrogen bonds on 1hls model_1 in heavy-atom-only mode.
        Baseline: rust/crates/proteindf-bridge/tests/test_sidechain_hydrogen_bond.rs lines 33-39
        and lines 71-74 (angle is None in heavy-atom-only mode).
        """
        ag = self._load_1hls()
        model_1 = ag["model_1"]

        hbonds = rs_br.calc_sidechain_hbonds(model_1)

        # Benchmark cases from test_sidechain_hydrogen_bond.rs lines 33-39:
        # (donor_subpath, donor_atom, acceptor_subpath, acceptor_atom, expected_distance)
        benchmark_cases = [
            ("A/8", "OG1", "A/4", "O", 2.5853),
            ("B/16", "NE2", "B/13", "O", 2.6420),
            ("A/2", "N", "A/19", "OH", 2.7983),
            ("B/10", "ND1", "B/7", "O", 2.9566),
            ("B/5", "NE2", "A/8", "O", 3.1916),
        ]

        for d_sub, d_atom, a_sub, a_atom, exp_dist in benchmark_cases:
            found = [
                hb
                for hb in hbonds
                if hb.donor_path.rstrip("/").endswith(d_sub)
                and hb.donor_atom == d_atom
                and hb.acceptor_path.rstrip("/").endswith(a_sub)
                and hb.acceptor_atom == a_atom
            ]
            self.assertEqual(
                len(found),
                1,
                f"Expected sidechain hbond {d_sub}.{d_atom} -> {a_sub}.{a_atom} not found",
            )
            hb = found[0]
            self.assertAlmostEqual(
                hb.distance,
                exp_dist,
                places=3,
                msg=f"Distance mismatch for {d_sub}.{d_atom} -> {a_sub}.{a_atom}",
            )
            # test_sidechain_hydrogen_bond.rs line 71: angle must be None in heavy-atom-only mode
            self.assertIsNone(hb.angle)
            self.assertIn("SidechainHydrogenBond(", repr(hb))

    def test_sidechain_hbonds_explicit_hydrogen_angle(self):
        """
        Verifies sidechain hydrogen bond angle filtering with explicit hydrogens.
        Baseline: rust/crates/proteindf-bridge/tests/test_sidechain_hydrogen_bond.rs lines 79-145:
        - Case 1 (lines 84-110): Angle 180° > 120° is detected with distance 2.5 and angle ~180.0°.
        - Case 2 (lines 112-144): Angle ~68.2° <= 120° is excluded in explicit-H mode (0 detected),
          but detected in heavy-atom-only mode (1 detected).
        """
        # Case 1: Ideal hydrogen angle (180 deg)
        model1 = rs_br.AtomGroup(name="model")
        chain1 = rs_br.AtomGroup(name="A")

        res1 = rs_br.AtomGroup(name="GLU")
        atom_o = rs_br.Atom(name="O", symbol="O", xyz=rs_br.Position(0.0, 0.0, 0.0))
        res1.set_atom("O", atom_o)
        chain1.set_group("1", res1)

        res2 = rs_br.AtomGroup(name="SER")
        atom_og = rs_br.Atom(name="OG", symbol="O", xyz=rs_br.Position(2.5, 0.0, 0.0))
        atom_hg = rs_br.Atom(name="HG", symbol="H", xyz=rs_br.Position(1.5, 0.0, 0.0))
        res2.set_atom("OG", atom_og)
        res2.set_atom("HG", atom_hg)
        chain1.set_group("2", res2)

        model1.set_group("A", chain1)

        # test_sidechain_hydrogen_bond.rs lines 103-109
        hbonds1 = rs_br.calc_sidechain_hbonds_with_options(model1, True)
        self.assertEqual(len(hbonds1), 1)
        self.assertEqual(hbonds1[0].donor_atom, "OG")
        self.assertEqual(hbonds1[0].acceptor_atom, "O")
        self.assertAlmostEqual(hbonds1[0].distance, 2.5, places=5)
        self.assertIsNotNone(hbonds1[0].angle)
        self.assertAlmostEqual(hbonds1[0].angle, 180.0, places=2)

        # Case 2: Unfavorable hydrogen angle (~68 deg <= 120 deg)
        model2 = rs_br.AtomGroup(name="model")
        chain2 = rs_br.AtomGroup(name="A")

        res2_1 = rs_br.AtomGroup(name="GLU")
        atom_o2 = rs_br.Atom(name="O", symbol="O", xyz=rs_br.Position(0.0, 0.0, 0.0))
        res2_1.set_atom("O", atom_o2)
        chain2.set_group("1", res2_1)

        res2_2 = rs_br.AtomGroup(name="SER")
        atom_og2 = rs_br.Atom(name="OG", symbol="O", xyz=rs_br.Position(2.5, 0.0, 0.0))
        atom_hg2 = rs_br.Atom(name="HG", symbol="H", xyz=rs_br.Position(2.5, 1.0, 0.0))
        res2_2.set_atom("OG", atom_og2)
        res2_2.set_atom("HG", atom_hg2)
        chain2.set_group("2", res2_2)

        model2.set_group("A", chain2)

        # test_sidechain_hydrogen_bond.rs lines 130-135: explicit-H mode -> 0 hbonds
        hbonds2_explicit = rs_br.calc_sidechain_hbonds_with_options(model2, True)
        self.assertEqual(len(hbonds2_explicit), 0)

        # test_sidechain_hydrogen_bond.rs lines 138-143: heavy-atom-only mode -> 1 hbond
        hbonds2_heavy = rs_br.calc_sidechain_hbonds(model2)
        self.assertEqual(len(hbonds2_heavy), 1)

    # -------------------------------------------------------------------------
    # 3. Secondary Structure (DSSP)
    # -------------------------------------------------------------------------
    def test_secondary_structure_1hls_chains(self):
        """
        Verifies 3-state secondary structure assignment on 1hls chains A and B.
        Baseline: rust/crates/proteindf-bridge/tests/test_secondary_structure.rs
        - Chain A (lines 34-59): len=21, residues 3..=6 Helix, 17..=19 Helix, rest Loop.
        - Chain B (lines 82-116): len=30, residues 9..=19 Helix, rest Loop.
        """
        ag = self._load_1hls()

        # Chain A
        chain_a = ag["model_1"]["A"]
        ss_a = rs_br.calc_secondary_structure(chain_a)
        # Baseline: test_secondary_structure.rs line 34
        self.assertEqual(len(ss_a), 21)

        # Baseline: test_secondary_structure.rs lines 37-59
        expected_a = [
            ("1", "GLY", "-"),
            ("2", "ILE", "-"),
            ("3", "VAL", "H"),
            ("4", "GLU", "H"),
            ("5", "GLN", "H"),
            ("6", "CYS", "H"),
            ("7", "CYS", "-"),
            ("8", "THR", "-"),
            ("9", "SER", "-"),
            ("10", "ILE", "-"),
            ("11", "CYS", "-"),
            ("12", "SER", "-"),
            ("13", "LEU", "-"),
            ("14", "TYR", "-"),
            ("15", "GLN", "-"),
            ("16", "LEU", "-"),
            ("17", "GLU", "H"),
            ("18", "ASN", "H"),
            ("19", "TYR", "H"),
            ("20", "CYS", "-"),
            ("21", "ASN", "-"),
        ]

        for i, (exp_k, exp_n, exp_c) in enumerate(expected_a):
            self.assertEqual(ss_a[i].residue_key, exp_k)
            self.assertEqual(ss_a[i].residue_name, exp_n)
            self.assertEqual(ss_a[i].code, exp_c)
            self.assertIn("SecondaryStructure(", repr(ss_a[i]))

        # Chain B
        chain_b = ag["model_1"]["B"]
        ss_b = rs_br.calc_secondary_structure(chain_b)
        # Baseline: test_secondary_structure.rs line 82
        self.assertEqual(len(ss_b), 30)

        # Baseline: test_secondary_structure.rs lines 85-116
        for i in range(30):
            res_num = i + 1
            exp_code = "H" if 9 <= res_num <= 19 else "-"
            self.assertEqual(ss_b[i].residue_key, str(res_num))
            self.assertEqual(ss_b[i].code, exp_code)

    def test_apply_secondary_structure_in_place(self):
        """
        Verifies apply_secondary_structure writes secondary_structure field in-place.
        Baseline: rust/crates/proteindf-bridge/tests/test_secondary_structure.rs lines 226-266.
        """
        ag = self._load_1hls()
        chain_a = ag["model_1"]["A"]

        # Before application, secondary_structure is None
        self.assertIsNone(chain_a["1"].secondary_structure)
        self.assertIsNone(chain_a["3"].secondary_structure)

        # Apply via function
        rs_br.apply_secondary_structure(chain_a)

        # Baseline: test_secondary_structure.rs lines 38-41
        self.assertEqual(chain_a["1"].secondary_structure, "-")
        self.assertEqual(chain_a["2"].secondary_structure, "-")
        self.assertEqual(chain_a["3"].secondary_structure, "H")
        self.assertEqual(chain_a["6"].secondary_structure, "H")
        self.assertEqual(chain_a["7"].secondary_structure, "-")

        # Apply via method on Chain B
        chain_b = ag["model_1"]["B"]
        chain_b.apply_secondary_structure()

        # Baseline: test_secondary_structure.rs lines 86-98
        self.assertEqual(chain_b["1"].secondary_structure, "-")
        self.assertEqual(chain_b["9"].secondary_structure, "H")
        self.assertEqual(chain_b["19"].secondary_structure, "H")
        self.assertEqual(chain_b["20"].secondary_structure, "-")

    def test_apply_secondary_structure_subtree_copy_independence(self):
        """
        Verifies that modifying a subtree copy does not affect the original tree.
        Consistent with design principle §4.3 item 2.
        """
        ag = self._load_1hls()
        chain_copy = ag["model_1"]["A"]  # Returns an owned subtree copy

        # Calling apply_secondary_structure on the copy
        chain_copy.apply_secondary_structure()
        self.assertEqual(chain_copy["3"].secondary_structure, "H")

        # The original ag must remain unchanged (None)
        self.assertIsNone(ag["model_1"]["A"]["3"].secondary_structure)

    # -------------------------------------------------------------------------
    # 4. CH-pi Interactions
    # -------------------------------------------------------------------------
    def test_ch_pi_interactions_1hls_benchmark(self):
        """
        Verifies CH-pi interactions detection on 1hls model_1.
        Baseline: rust/crates/proteindf-bridge/tests/test_ch_pi.rs lines 92-145:
        - A/10.CG1 -> B/5: dist 3.400, angle 23.79 (detected)
        - B/15.CB -> B/24: dist 3.793, angle 23.69 (detected)
        - B/17.CA -> B/16: dist 4.059, angle 15.62 (detected)
        - A/10.CD1 -> B/5: dist 4.074, angle 41.11 (excluded at 40° max angle)
        - B/27.C -> B/26: dist 4.578, angle 11.75 (excluded at 4.5 Å max distance)
        """
        ag = self._load_1hls()
        model_1 = ag["model_1"]

        interactions = rs_br.calc_ch_pi_interactions(model_1)

        candidates = [
            ("A/10", "CG1", "B/5", 3.400, 23.79, True),
            ("B/15", "CB", "B/24", 3.793, 23.69, True),
            ("B/17", "CA", "B/16", 4.059, 15.62, True),
            ("A/10", "CD1", "B/5", 4.074, 41.11, False),
            ("B/27", "C", "B/26", 4.578, 11.75, False),
        ]

        for c_sub, c_atom, r_sub, exp_dist, exp_angle, should_detect in candidates:
            found = [
                i
                for i in interactions
                if i.carbon_path.rstrip("/").endswith(c_sub)
                and i.carbon_atom == c_atom
                and i.ring_path.rstrip("/").endswith(r_sub)
            ]
            if should_detect:
                self.assertEqual(
                    len(found),
                    1,
                    f"Expected CH-pi interaction {c_sub}.{c_atom} -> {r_sub} not found",
                )
                item = found[0]
                self.assertAlmostEqual(item.distance, exp_dist, places=3)
                self.assertAlmostEqual(item.angle, exp_angle, places=2)
                self.assertIn("ChPiInteraction(", repr(item))
            else:
                self.assertEqual(
                    len(found),
                    0,
                    f"Candidate {c_sub}.{c_atom} -> {r_sub} should not be detected at default thresholds",
                )

    def test_ch_pi_tunable_thresholds(self):
        """
        Verifies tunable thresholds for CH-pi interaction detection.
        Baseline: rust/crates/proteindf-bridge/tests/test_ch_pi.rs lines 149-179:
        - Relaxing angle threshold to 42.0° detects A/10.CD1 -> B/5 (angle 41.11°).
        - Relaxing distance threshold to 4.6 Å detects B/27.C -> B/26 (dist 4.578 Å).
        """
        ag = self._load_1hls()
        model_1 = ag["model_1"]

        # 1. Relax angle threshold to 42.0 deg (test_ch_pi.rs line 157)
        inter_angle_42 = rs_br.calc_ch_pi_interactions(model_1, max_angle_deg=42.0)
        found_angle = [
            i
            for i in inter_angle_42
            if i.carbon_path.rstrip("/").endswith("A/10")
            and i.carbon_atom == "CD1"
            and i.ring_path.rstrip("/").endswith("B/5")
        ]
        self.assertEqual(len(found_angle), 1)

        # 2. Relax distance threshold to 4.6 A via calc_ch_pi_interactions_with_thresholds (test_ch_pi.rs line 169)
        inter_dist_46 = rs_br.calc_ch_pi_interactions_with_thresholds(model_1, 4.6, 40.0)
        found_dist = [
            i
            for i in inter_dist_46
            if i.carbon_path.rstrip("/").endswith("B/27")
            and i.carbon_atom == "C"
            and i.ring_path.rstrip("/").endswith("B/26")
        ]
        self.assertEqual(len(found_dist), 1)

    # -------------------------------------------------------------------------
    # 5. InteractionSet
    # -------------------------------------------------------------------------
    def test_interaction_set_detect_all_1hls(self):
        """
        Verifies InteractionSet.detect_all on 1hls model_1.
        Baseline: rust/crates/proteindf-bridge/tests/test_interaction_set.rs lines 14-109:
        - Disulfide: 3 (lines 23-25)
        - Salt bridge: 0 (line 35)
        - Hydrogen bonds: 28 total (22 backbone + 6 sidechain) (lines 58-61)
        - CH-pi: 7 (lines 85-86)
        - Total: 38 (line 108)
        """
        ag = self._load_1hls()
        model_1 = ag["model_1"]

        iset = rs_br.InteractionSet.detect_all(model_1)

        # Total interactions: test_interaction_set.rs line 108
        self.assertEqual(len(iset), 38)
        self.assertFalse(iset.is_empty())

        # Disulfide: test_interaction_set.rs line 24-25
        self.assertEqual(iset.count_by_kind("disulfide"), 3)
        disulfides = iset.filter_by_kind("disulfide")
        self.assertEqual(len(disulfides), 3)
        for ds in disulfides:
            self.assertEqual(len(ds.atoms), 2)
            self.assertIsNotNone(ds.distance)
            self.assertTrue(1.9 < ds.distance < 2.3)

        # Salt bridge: test_interaction_set.rs line 35
        self.assertEqual(iset.count_by_kind("salt_bridge"), 0)
        self.assertEqual(len(iset.filter_by_kind("salt_bridge")), 0)

        # Hydrogen bonds: test_interaction_set.rs lines 58-61
        self.assertEqual(iset.count_by_kind("hydrogen_bond"), 28)
        hbonds = iset.filter_by_kind("hydrogen_bond")
        self.assertEqual(len(hbonds), 28)
        bb_count = sum(1 for h in hbonds if h.donor_acceptor_role and "backbone" in h.donor_acceptor_role)
        sc_count = sum(1 for h in hbonds if h.donor_acceptor_role and "sidechain" in h.donor_acceptor_role)
        self.assertEqual(bb_count, 22)
        self.assertEqual(sc_count, 6)

        # CH-pi: test_interaction_set.rs line 85-86
        self.assertEqual(iset.count_by_kind("ch_pi"), 7)
        ch_pi_list = iset.filter_by_kind("ch_pi")
        self.assertEqual(len(ch_pi_list), 7)

        # Repr
        repr_str = repr(iset)
        self.assertIn("InteractionSet(total=38", repr_str)
        self.assertIn("disulfide=3", repr_str)
        self.assertIn("hydrogen_bond=28", repr_str)
        self.assertIn("ch_pi=7", repr_str)

    def test_interaction_set_tuning(self):
        """
        Verifies InteractionSet.detect_all with custom CH-pi thresholds.
        Baseline: rust/crates/proteindf-bridge/tests/test_interaction_set.rs lines 112-129.
        """
        ag = self._load_1hls()
        model_1 = ag["model_1"]

        set_default = rs_br.InteractionSet.detect_all(model_1)
        self.assertEqual(set_default.count_by_kind("ch_pi"), 7)

        set_relaxed = rs_br.InteractionSet.detect_all(model_1, ch_pi_max_distance=5.0, ch_pi_max_angle_deg=50.0)
        self.assertGreater(set_relaxed.count_by_kind("ch_pi"), 7)

        set_strict = rs_br.InteractionSet.detect_all(model_1, ch_pi_max_distance=3.5, ch_pi_max_angle_deg=20.0)
        self.assertLess(set_strict.count_by_kind("ch_pi"), 7)

    def test_interaction_set_msgpack_roundtrip(self):
        """
        Verifies InteractionSet serialization and deserialization via MessagePack (in-memory and file).
        Baseline: rust/crates/proteindf-bridge/tests/test_interaction_set.rs lines 132-155.
        """
        ag = self._load_1hls()
        model_1 = ag["model_1"]
        original = rs_br.InteractionSet.detect_all(model_1)

        # 1. In-memory roundtrip
        bytes_data = original.to_msgpack()
        self.assertIsInstance(bytes_data, bytes)
        restored = rs_br.InteractionSet.from_msgpack(bytes_data)
        self.assertEqual(original, restored)
        self.assertEqual(len(restored), len(original))

        # 2. File roundtrip
        with tempfile.NamedTemporaryFile(suffix=".msgpack", delete=False) as tmp:
            tmp_path = tmp.name

        try:
            original.save_msgpack(tmp_path)
            loaded = rs_br.InteractionSet.load_msgpack(tmp_path)
            self.assertEqual(original, loaded)
            self.assertEqual(len(loaded), len(original))
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_interaction_set_yaml_roundtrip(self):
        """
        Verifies InteractionSet serialization and deserialization via YAML (in-memory and file).
        Baseline: rust/crates/proteindf-bridge/tests/test_interaction_set.rs lines 158-183.
        """
        ag = self._load_1hls()
        model_1 = ag["model_1"]
        original = rs_br.InteractionSet.detect_all(model_1)

        # 1. In-memory roundtrip
        yaml_str = original.to_yaml()
        self.assertIsInstance(yaml_str, str)
        self.assertIn("disulfide", yaml_str)
        self.assertIn("hydrogen_bond", yaml_str)
        self.assertIn("ch_pi", yaml_str)

        restored = rs_br.InteractionSet.from_yaml(yaml_str)
        self.assertEqual(original, restored)
        self.assertEqual(len(restored), len(original))

        # 2. File roundtrip
        with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as tmp:
            tmp_path = tmp.name

        try:
            original.save_yaml(tmp_path)
            loaded = rs_br.InteractionSet.load_yaml(tmp_path)
            self.assertEqual(original, loaded)
            self.assertEqual(len(loaded), len(original))
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_interaction_set_immutability_and_protection(self):
        """
        Verifies that returned interaction lists and atom path lists cannot mutate internal state.
        Follows PR#46 review lessons on collection immutability.
        """
        ag = self._load_1hls()
        model_1 = ag["model_1"]
        iset = rs_br.InteractionSet.detect_all(model_1)

        # 1. Mutating returned interactions list
        interactions = iset.interactions
        init_len = len(interactions)
        interactions.append("fake")
        self.assertEqual(len(interactions), init_len + 1)
        self.assertEqual(len(iset.interactions), init_len)

        # 2. Mutating atoms list of an Interaction
        first_inter = iset.interactions[0]
        init_atoms_len = len(first_inter.atoms)
        first_inter.atoms.append("fake/atom")
        self.assertEqual(len(first_inter.atoms), init_atoms_len)

    def test_interaction_set_error_handling(self):
        """
        Verifies that invalid kinds or corrupt files raise BrValueError / BrError exceptions.
        """
        iset = rs_br.InteractionSet()
        with self.assertRaises(rs_br.BrValueError):
            iset.filter_by_kind("invalid_kind")

        with self.assertRaises(rs_br.BrValueError):
            iset.count_by_kind("invalid_kind")

        with self.assertRaises(rs_br.BrValueError):
            rs_br.InteractionSet.from_yaml("invalid: yaml: [}")

        with self.assertRaises(rs_br.BrValueError):
            rs_br.InteractionSet.from_msgpack(b"corrupt bytes")

        with self.assertRaises(rs_br.BrError):
            rs_br.InteractionSet.load_yaml("/nonexistent/path/interactions.yaml")

        with self.assertRaises(rs_br.BrError):
            rs_br.InteractionSet.load_msgpack("/nonexistent/path/interactions.msgpack")


if __name__ == "__main__":
    unittest.main()
