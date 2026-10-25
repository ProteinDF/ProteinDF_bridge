#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

"""Test suite for missing heavy-atom completion (rs.AtomGroup.complete_missing_atoms).

Validates:
1. ALA: CB completion (added == ["CB"], is_complete, CB-CA distance ~1.52 Å, orig coords unchanged).
2. LYS: multiple heavy-atom completion (CD, CE, NZ).
3. Insufficient anchors: <3 anchors and collinear anchors.
4. PoorFit and max_rmsd: distorted anchor geometry fails with PoorFit; succeeds with relaxed max_rmsd.
5. Unknown component: NoTemplate status with built-in DB; completed with custom CCD DB.
6. Idempotence: complete structure unchanged, second call adds nothing.
7. Empty tree: safe no-op.
8. Hierarchy: callable on leaf residue, chain, model, and root AtomGroup.
9. Serialization: to_dict JSON serializable, __str__, and __repr__.
10. Workflow: inspection -> completion -> hydrogenation reduces hydrogen defects compared to uncompleted.
11. Real fixture 2FB4: heavy-atom completion RMSD < 0.8 Å against ground truth.
"""

import json
import math
import os
import tempfile
import unittest

import proteindf_bridge.rs as rs_br

DATA_DIR = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "rust",
        "crates",
        "proteindf-bridge",
        "tests",
        "data",
    )
)


def build_residue_from_template(
    comp_id: str,
    include_hydrogens: bool,
    exclude: tuple[str, ...] = (),
) -> rs_br.AtomGroup:
    """Builds a residue AtomGroup from the builtin CCD template, excluding given atom names."""
    db = rs_br.CcdTemplateDb.builtin()
    tmpl = db.lookup(comp_id)
    assert tmpl is not None, f"Template for {comp_id} not found in builtin CCD DB"
    group = rs_br.AtomGroup(comp_id)
    exclude_set = set(exclude)
    idx = 1
    for ccd_atom in tmpl.atoms:
        if not include_hydrogens and ccd_atom.is_hydrogen:
            continue
        if ccd_atom.name in exclude_set:
            continue
        xyz = ccd_atom.ideal_xyz or (0.0, 0.0, 0.0)
        atom = rs_br.Atom(symbol=ccd_atom.element, name=ccd_atom.name, xyz=rs_br.Position(*xyz))
        group.set_atom(f"{idx}_{ccd_atom.name}", atom)
        idx += 1
    return group


class TestRsCompleteAtoms(unittest.TestCase):
    """Unit tests for complete_missing_atoms pyo3 binding."""

    def setUp(self):
        self.tmp_files = []

    def tearDown(self):
        for path in self.tmp_files:
            if os.path.exists(path):
                try:
                    os.remove(path)
                except OSError:
                    pass

    def test_ala_cb_completion(self):
        """ALA CB completion, check_missing_atoms is_complete, distance, and coordinate stability."""
        db = rs_br.CcdTemplateDb.builtin()

        # Build ALA missing CB and optional OXT, shifted in space
        ala = build_residue_from_template("ALA", False, exclude=("CB", "OXT"))
        shift = rs_br.Position(3.5, -2.0, 1.25)
        for _, atom in ala.atoms():
            atom.xyz = rs_br.Position(atom.xyz.x + shift.x, atom.xyz.y + shift.y, atom.xyz.z + shift.z)

        orig_coords = {atom.name: rs_br.Position(atom.xyz.x, atom.xyz.y, atom.xyz.z) for _, atom in ala.atoms()}

        rep = ala.complete_missing_atoms(db=db, max_rmsd=1.0)
        self.assertEqual(len(rep.residues), 1)
        self.assertEqual(len(rep.completed), 1)
        self.assertEqual(len(rep.skipped), 0)
        self.assertEqual(rep.total_added, 1)

        r = rep.residues[0]
        self.assertEqual(r.comp_id, "ALA")
        self.assertEqual(r.status, "completed")
        self.assertEqual(r.added, ["CB"])
        self.assertEqual(r.anchors, 4)
        self.assertIsNotNone(r.rmsd)
        self.assertLess(r.rmsd, 1e-5)

        # check_missing_atoms reports is_complete after completion
        check_rep = ala.check_missing_atoms(db)
        self.assertTrue(check_rep.is_complete)
        self.assertEqual(check_rep.total_missing, 0)
        self.assertFalse(ala.has_atom("OXT"))  # Optional atom not added

        # CB atom added into the tree
        self.assertTrue(ala.has_atom("CB"))
        ca = ala.get_atom("CA")
        cb = ala.get_atom("CB")
        dist = cb.xyz.distance_from(ca.xyz)
        self.assertAlmostEqual(dist, 1.524, delta=0.03)

        # Original coordinates remain unchanged
        for name, orig_xyz in orig_coords.items():
            cur_atom = ala.get_atom(name)
            self.assertAlmostEqual(cur_atom.xyz.x, orig_xyz.x, places=5)
            self.assertAlmostEqual(cur_atom.xyz.y, orig_xyz.y, places=5)
            self.assertAlmostEqual(cur_atom.xyz.z, orig_xyz.z, places=5)

    def test_lys_multiple_atoms(self):
        """LYS missing multiple heavy atoms (CD, CE, NZ)."""
        full_lys = build_residue_from_template("LYS", False, exclude=("OXT",))
        lys = build_residue_from_template("LYS", False, exclude=("CD", "CE", "NZ", "OXT"))

        rep = lys.complete_missing_atoms()
        self.assertEqual(len(rep.residues), 1)
        self.assertEqual(len(rep.completed), 1)
        self.assertEqual(len(rep.skipped), 0)
        self.assertEqual(rep.total_added, 3)

        r = rep.residues[0]
        self.assertEqual(r.comp_id, "LYS")
        self.assertEqual(r.status, "completed")
        self.assertEqual(r.added, ["CD", "CE", "NZ"])
        self.assertEqual(r.anchors, 6)
        self.assertIsNotNone(r.rmsd)
        self.assertLess(r.rmsd, 1e-5)

        self.assertTrue(lys.check_missing_atoms().is_complete)

        for name in ["CD", "CE", "NZ"]:
            self.assertTrue(lys.has_atom(name))
            added_xyz = lys.get_atom(name).xyz
            expected_xyz = full_lys.get_atom(name).xyz
            self.assertAlmostEqual(added_xyz.distance_from(expected_xyz), 0.0, places=4)

    def test_insufficient_anchors(self):
        """Insufficient anchors: fewer than 3 anchors and collinear anchors."""
        # 1. Fewer than 3 anchors (only N, CA)
        ala_two = build_residue_from_template("ALA", False, exclude=("C", "O", "CB", "OXT"))
        raw_before = ala_two.get_raw_data()
        rep_two = ala_two.complete_missing_atoms()

        self.assertEqual(len(rep_two.residues), 1)
        self.assertEqual(len(rep_two.completed), 0)
        self.assertEqual(len(rep_two.skipped), 1)
        self.assertEqual(rep_two.total_added, 0)
        self.assertEqual(ala_two.get_raw_data(), raw_before)

        r_two = rep_two.residues[0]
        self.assertEqual(r_two.status, "insufficient_anchors")
        self.assertEqual(r_two.added, [])
        self.assertEqual(r_two.anchors, 2)
        self.assertIsNone(r_two.rmsd)

        # 2. Collinear anchors (N, CA, C on a straight line)
        ala_col = rs_br.AtomGroup("ALA")
        ala_col.set_atom("N", rs_br.Atom(symbol="N", name="N", xyz=rs_br.Position(0.0, 0.0, 0.0)))
        ala_col.set_atom("CA", rs_br.Atom(symbol="C", name="CA", xyz=rs_br.Position(1.5, 0.0, 0.0)))
        ala_col.set_atom("C", rs_br.Atom(symbol="C", name="C", xyz=rs_br.Position(3.0, 0.0, 0.0)))
        raw_col_before = ala_col.get_raw_data()

        rep_col = ala_col.complete_missing_atoms()
        self.assertEqual(len(rep_col.completed), 0)
        self.assertEqual(len(rep_col.skipped), 1)
        self.assertEqual(ala_col.get_raw_data(), raw_col_before)

        r_col = rep_col.residues[0]
        self.assertEqual(r_col.status, "insufficient_anchors")
        self.assertEqual(r_col.added, [])
        self.assertEqual(r_col.anchors, 3)
        self.assertIsNone(r_col.rmsd)

    def test_poor_fit_and_max_rmsd(self):
        """PoorFit when anchor superposition exceeds max_rmsd; succeeds when max_rmsd relaxed."""
        ala = build_residue_from_template("ALA", False, exclude=("CB", "OXT"))
        # Swap N and O coordinates to disrupt geometry
        n_atom = ala.get_atom("N")
        o_atom = ala.get_atom("O")
        n_pos = rs_br.Position(n_atom.xyz.x, n_atom.xyz.y, n_atom.xyz.z)
        o_pos = rs_br.Position(o_atom.xyz.x, o_atom.xyz.y, o_atom.xyz.z)
        n_atom.xyz = o_pos
        o_atom.xyz = n_pos

        raw_before = ala.get_raw_data()

        # Strict default max_rmsd=1.0 -> PoorFit, tree not modified
        rep_strict = ala.complete_missing_atoms()
        self.assertEqual(len(rep_strict.completed), 0)
        self.assertEqual(len(rep_strict.skipped), 1)
        self.assertEqual(rep_strict.total_added, 0)
        self.assertEqual(ala.get_raw_data(), raw_before)

        r_strict = rep_strict.residues[0]
        self.assertEqual(r_strict.status, "poor_fit")
        self.assertEqual(r_strict.added, [])
        self.assertEqual(r_strict.anchors, 4)
        self.assertIsNotNone(r_strict.rmsd)
        self.assertGreater(r_strict.rmsd, 1.0)
        self.assertFalse(ala.has_atom("CB"))

        # Relaxed max_rmsd=5.0 -> Completed
        rep_relaxed = ala.complete_missing_atoms(max_rmsd=5.0)
        self.assertEqual(len(rep_relaxed.completed), 1)
        self.assertEqual(len(rep_relaxed.skipped), 0)
        self.assertEqual(rep_relaxed.total_added, 1)

        r_relaxed = rep_relaxed.residues[0]
        self.assertEqual(r_relaxed.status, "completed")
        self.assertEqual(r_relaxed.added, ["CB"])
        self.assertTrue(ala.has_atom("CB"))

    def test_unknown_component_and_custom_db(self):
        """Unknown component yields no_template in builtin DB; completed with custom CCD DB."""
        lig = rs_br.AtomGroup("LIG")
        # Template will define C1, C2, O1, C3, N1 with ideal coordinates
        # Supply C1, C2, O1, C3 as anchors (4 non-coplanar atoms)
        lig.set_atom("C1", rs_br.Atom(symbol="C", name="C1", xyz=rs_br.Position(0.0, 0.0, 0.0)))
        lig.set_atom("C2", rs_br.Atom(symbol="C", name="C2", xyz=rs_br.Position(1.5, 0.0, 0.0)))
        lig.set_atom("O1", rs_br.Atom(symbol="O", name="O1", xyz=rs_br.Position(1.5, 1.4, 0.0)))
        lig.set_atom("C3", rs_br.Atom(symbol="C", name="C3", xyz=rs_br.Position(0.0, 0.8, 1.2)))

        # Built-in DB: no template
        rep_unk = lig.complete_missing_atoms()
        self.assertEqual(len(rep_unk.completed), 0)
        self.assertEqual(len(rep_unk.skipped), 1)
        self.assertEqual(rep_unk.total_added, 0)
        r_unk = rep_unk.residues[0]
        self.assertEqual(r_unk.status, "no_template")
        self.assertEqual(r_unk.added, [])

        # Create custom CCD DB with LIG template containing missing atom N1
        cif_content = """\
data_comp_LIG
_chem_comp.id LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
LIG C1 C 0.000 0.000 0.000
LIG C2 C 1.500 0.000 0.000
LIG O1 O 1.500 1.400 0.000
LIG C3 C 0.000 0.800 1.200
LIG N1 N 0.000 -1.400 0.500
"""
        with tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False) as f:
            f.write(cif_content)
            tmp_path = f.name
        self.tmp_files.append(tmp_path)

        custom_db = rs_br.CcdTemplateDb.builtin()
        count = custom_db.add_from_file(tmp_path)
        self.assertEqual(count, 1)

        rep_custom = lig.complete_missing_atoms(db=custom_db)
        self.assertEqual(len(rep_custom.completed), 1)
        self.assertEqual(len(rep_custom.skipped), 0)
        self.assertEqual(rep_custom.total_added, 1)

        r_custom = rep_custom.residues[0]
        self.assertEqual(r_custom.status, "completed")
        self.assertEqual(r_custom.added, ["N1"])
        self.assertTrue(lig.has_atom("N1"))
        n1_xyz = lig.get_atom("N1").xyz
        self.assertAlmostEqual(n1_xyz.x, 0.0, places=3)
        self.assertAlmostEqual(n1_xyz.y, -1.4, places=3)
        self.assertAlmostEqual(n1_xyz.z, 0.5, places=3)

    def test_complete_structure_and_idempotence(self):
        """Complete structure is unchanged; repeated call adds nothing."""
        ala = build_residue_from_template("ALA", False, exclude=("OXT",))
        raw_before = ala.get_raw_data()

        rep1 = ala.complete_missing_atoms()
        self.assertEqual(rep1.total_added, 0)
        self.assertEqual(len(rep1.completed), 0)
        self.assertEqual(len(rep1.skipped), 0)
        self.assertEqual(len(rep1.residues), 1)
        self.assertEqual(rep1.residues[0].status, "unchanged")
        self.assertEqual(ala.get_raw_data(), raw_before)

        # Now test idempotence on incomplete residue after completion
        ala_inc = build_residue_from_template("ALA", False, exclude=("CB", "OXT"))
        rep_first = ala_inc.complete_missing_atoms()
        self.assertEqual(rep_first.total_added, 1)
        self.assertEqual(rep_first.residues[0].status, "completed")

        raw_after_first = ala_inc.get_raw_data()
        rep_second = ala_inc.complete_missing_atoms()
        self.assertEqual(rep_second.total_added, 0)
        self.assertEqual(len(rep_second.completed), 0)
        self.assertEqual(len(rep_second.skipped), 0)
        self.assertEqual(rep_second.residues[0].status, "unchanged")
        self.assertEqual(ala_inc.get_raw_data(), raw_after_first)

    def test_empty_tree(self):
        """Empty AtomGroup returns empty report safely."""
        empty = rs_br.AtomGroup("empty")
        rep = empty.complete_missing_atoms()
        self.assertEqual(rep.total_added, 0)
        self.assertEqual(len(rep.residues), 0)
        self.assertEqual(len(rep.completed), 0)
        self.assertEqual(len(rep.skipped), 0)

    def test_callable_across_hierarchy_levels(self):
        """Callable on leaf residue, chain, model, and root levels."""
        # Build hierarchy: root / model_1 / A / 1 (ALA missing CB), 2 (GLY complete)
        root = rs_br.AtomGroup()
        model = rs_br.AtomGroup("model_1")
        chain = rs_br.AtomGroup("A")
        ala = build_residue_from_template("ALA", False, exclude=("CB", "OXT"))
        gly = build_residue_from_template("GLY", False, exclude=("OXT",))

        chain.set_group("1", ala)
        chain.set_group("2", gly)
        model.set_group("A", chain)
        root.set_group("model_1", model)

        # 1. Calling on leaf residue directly
        ala_leaf = build_residue_from_template("ALA", False, exclude=("CB", "OXT"))
        rep_leaf = ala_leaf.complete_missing_atoms()
        self.assertEqual(rep_leaf.total_added, 1)
        self.assertTrue(ala_leaf.has_atom("CB"))

        # 2. Calling on root AtomGroup
        rep_root = root.complete_missing_atoms()
        self.assertEqual(rep_root.total_added, 1)
        self.assertEqual(len(rep_root.completed), 1)
        self.assertEqual(len(rep_root.residues), 2)
        res_map = {r.path: r for r in rep_root.residues}
        self.assertEqual(res_map["/model_1/A/1/"].status, "completed")
        self.assertEqual(res_map["/model_1/A/2/"].status, "unchanged")

        # Tree is mutated in-place
        ala_in_tree = root["model_1"]["A"]["1"]
        self.assertTrue(ala_in_tree.has_atom("CB"))

    def test_serialization_and_formatting(self):
        """to_dict produces valid JSON; __str__ and __repr__ formatting."""
        ala = build_residue_from_template("ALA", False, exclude=("CB", "OXT"))
        lig = rs_br.AtomGroup("LIG")
        lig.set_atom("1_C", rs_br.Atom(symbol="C", name="C", xyz=rs_br.Position(0.0, 0.0, 0.0)))

        root = rs_br.AtomGroup()
        root.set_group("ala", ala)
        root.set_group("lig", lig)

        rep = root.complete_missing_atoms()
        d = rep.to_dict()

        # JSON round-trip
        json_str = json.dumps(d)
        self.assertIsInstance(json_str, str)
        d_loaded = json.loads(json_str)
        self.assertEqual(d_loaded["total_added"], 1)
        self.assertEqual(len(d_loaded["completed"]), 1)
        self.assertEqual(len(d_loaded["skipped"]), 1)
        self.assertEqual(len(d_loaded["residues"]), 2)

        # ResidueCompletion to_dict
        r_comp = rep.completed[0]
        r_dict = r_comp.to_dict()
        self.assertEqual(r_dict["comp_id"], "ALA")
        self.assertEqual(r_dict["status"], "completed")
        self.assertEqual(r_dict["added"], ["CB"])
        self.assertEqual(r_dict["anchors"], 4)
        self.assertIsNotNone(r_dict["rmsd"])

        # __str__ output
        s = str(rep)
        self.assertIn("CompletionReport: 1 heavy atom(s) added", s)
        self.assertIn("Completed residues:", s)
        self.assertIn("Skipped residues:", s)
        self.assertIn("ALA", s)
        self.assertIn("LIG", s)

        # __repr__ output
        r_str = repr(rep)
        self.assertIn("CompletionReport(total_added=1, completed=1, skipped=1, residues=2)", r_str)
        self.assertIn("ResidueCompletion", repr(r_comp))

    def test_inspection_completion_hydrogenation_workflow(self):
        """Inspection -> completion -> hydrogenation reduces hydrogen defects compared to uncompleted."""
        # Build ALA missing CB
        # Case A: add hydrogens without heavy-atom completion
        ala_a = build_residue_from_template("ALA", False, exclude=("CB", "OXT"))
        ala_a.add_missing_hydrogens()
        check_a = ala_a.check_missing_atoms(include_hydrogens=True)
        # Without CB, HB1, HB2, HB3 cannot be placed, plus CB itself is missing
        missing_a = check_a.total_missing
        self.assertGreaterEqual(missing_a, 4)

        # Case B: complete heavy atoms first, then add hydrogens
        ala_b = build_residue_from_template("ALA", False, exclude=("CB", "OXT"))
        rep_comp = ala_b.complete_missing_atoms()
        self.assertEqual(rep_comp.total_added, 1)

        ala_b.add_missing_hydrogens()
        check_b = ala_b.check_missing_atoms(include_hydrogens=True)
        missing_b = check_b.total_missing

        # Defects in case B must be strictly fewer than in case A
        self.assertLess(missing_b, missing_a)
        # In fact, ALA with all standard hydrogens added should be complete (0 missing)
        self.assertEqual(missing_b, 0)
        self.assertTrue(check_b.is_complete)

    def test_real_fixture_2fb4(self):
        """Real fixture 2FB4: remove heavy atoms across 5 residues, complete, and verify RMSD < 0.8 Å."""
        path_2fb4 = os.path.join(DATA_DIR, "2FB4.cif")
        cif_2fb4 = rs_br.SimpleMmcif(path_2fb4)
        ag_2fb4 = cif_2fb4.get_structure_atomgroup()

        targets_2fb4 = [
            ("H", "14", ["CB", "CG"]),  # PRO
            ("H", "33", ["CB"]),        # ALA
            ("H", "40", ["CB"]),        # ALA
            ("H", "41", ["CB", "CG"]),  # PRO
            ("L", "13", ["CB"]),        # THR
        ]

        orig_coords = []
        for chain_id, res_id, atom_names in targets_2fb4:
            res = ag_2fb4["model_1"][chain_id][res_id]
            res_orig = []
            for aname in atom_names:
                # Find atom key in residue
                target_key = None
                for key, atom in res.atoms():
                    if atom.name == aname:
                        target_key = key
                        res_orig.append((aname, rs_br.Position(atom.xyz.x, atom.xyz.y, atom.xyz.z)))
                        break
                self.assertIsNotNone(target_key, f"Target atom {aname} not found in {chain_id}/{res_id}")
                res.remove_atom(target_key)
            orig_coords.append((chain_id, res_id, res_orig))

        total_removed = sum(len(atoms) for _, _, atoms in orig_coords)
        self.assertEqual(total_removed, 7)

        # Inspection confirms 7 missing atoms
        check_before = ag_2fb4.check_missing_atoms()
        self.assertEqual(check_before.total_missing, 7)

        # Complete missing atoms
        rep = ag_2fb4.complete_missing_atoms()
        self.assertEqual(len(rep.residues), 550)
        self.assertEqual(len(rep.completed), 5)
        self.assertEqual(len(rep.skipped), 0)
        self.assertEqual(rep.total_added, 7)

        # Inspection confirms fully complete
        check_after = ag_2fb4.check_missing_atoms()
        self.assertTrue(check_after.is_complete)
        self.assertEqual(check_after.total_missing, 0)

        # Compute reconstructed heavy atom RMSD against ground truth
        total_sq_err = 0.0
        max_atom_dist = 0.0
        max_res_rmsd = 0.0

        for chain_id, res_id, atoms in orig_coords:
            res = ag_2fb4["model_1"][chain_id][res_id]
            res_sq = 0.0
            for aname, orig_pos in atoms:
                self.assertTrue(res.has_atom(aname), f"Reconstructed atom {aname} missing in {chain_id}/{res_id}")
                reconstructed_atom = res.get_atom(aname)
                dist = orig_pos.distance_from(reconstructed_atom.xyz)
                max_atom_dist = max(max_atom_dist, dist)
                res_sq += dist * dist
                total_sq_err += dist * dist
            res_rmsd = math.sqrt(res_sq / len(atoms))
            max_res_rmsd = max(max_res_rmsd, res_rmsd)

        overall_rmsd = math.sqrt(total_sq_err / total_removed)

        print(
            f"\n[2FB4 Completion Summary] removed={total_removed}, added={rep.total_added}, "
            f"completed_residues={len(rep.completed)}, "
            f"overall_rmsd={overall_rmsd:.4f} Å, max_residue_rmsd={max_res_rmsd:.4f} Å, "
            f"max_atom_dist={max_atom_dist:.4f} Å"
        )

        self.assertLess(
            overall_rmsd,
            0.8,
            f"Overall reconstructed RMSD {overall_rmsd:.4f} Å must be < 0.8 Å",
        )
        self.assertLess(max_res_rmsd, 0.8)


if __name__ == "__main__":
    unittest.main()
