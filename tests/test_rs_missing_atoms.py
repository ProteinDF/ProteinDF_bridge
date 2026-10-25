#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

"""Test suite for missing-atoms inspection (rs.AtomGroup.check_missing_atoms).

Validates:
1. ALA: complete without OXT, CB missing, hydrogen mode (HA, HB1-HB3 missing).
2. HOH: complete in heavy-atom mode, incomplete (H1, H2 missing) in hydrogen mode.
3. Unknown components: NoTemplate status with built-in DB; inspection with custom CCD DB.
4. DA: missing N9, 5'-terminal phosphate detection, internal residue lacking phosphate.
5. Unexpected atoms and duplicate atom names handling.
6. CYS disulfide bonds: HG becomes optional within distance threshold across hierarchy.
7. Empty tree, subtrees, and immutability (get_raw_data unchanged).
8. Serialization (to_dict, json.dumps), __str__, and __repr__.
9. Real fixtures 2FB4 and 2MGO with benchmark timing.
"""

import json
import os
import tempfile
import timeit
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


class TestRsMissingAtoms(unittest.TestCase):
    """Unit tests for check_missing_atoms pyo3 binding."""

    def setUp(self):
        self.tmp_files = []

    def tearDown(self):
        for path in self.tmp_files:
            if os.path.exists(path):
                try:
                    os.remove(path)
                except OSError:
                    pass

    def test_ala_heavy_and_hydrogens(self):
        """ALA complete without OXT, CB missing, and hydrogen-mode inspections."""
        db = rs_br.CcdTemplateDb.builtin()

        # 1. ALA without OXT is complete in heavy-atom mode
        ala_no_oxt = build_residue_from_template("ALA", False, exclude=("OXT",))
        rep = ala_no_oxt.check_missing_atoms(db)
        self.assertEqual(len(rep.residues), 1)
        self.assertTrue(rep.is_complete)
        self.assertEqual(rep.total_missing, 0)
        self.assertEqual(len(rep.incomplete), 0)
        self.assertEqual(len(rep.no_template), 0)
        self.assertEqual(rep.no_template_components, {})

        r = rep.residues[0]
        self.assertEqual(r.path, "/")
        self.assertEqual(r.comp_id, "ALA")
        self.assertEqual(r.status, "complete")
        self.assertEqual(r.expected, 5)
        self.assertEqual(r.present, 5)
        self.assertEqual(r.missing, [])
        self.assertEqual(r.unexpected, [])

        # Default db=None uses the global default CCD DB
        rep_default_db = ala_no_oxt.check_missing_atoms()
        self.assertTrue(rep_default_db.is_complete)
        self.assertEqual(rep_default_db.total_missing, 0)

        # 2. ALA missing CB is incomplete
        ala_no_cb = build_residue_from_template("ALA", False, exclude=("CB", "OXT"))
        rep_cb = ala_no_cb.check_missing_atoms(db)
        self.assertFalse(rep_cb.is_complete)
        self.assertEqual(rep_cb.total_missing, 1)
        self.assertEqual(len(rep_cb.incomplete), 1)

        r_cb = rep_cb.residues[0]
        self.assertEqual(r_cb.status, "incomplete")
        self.assertEqual(r_cb.expected, 5)
        self.assertEqual(r_cb.present, 4)
        self.assertEqual(r_cb.missing, ["CB"])
        self.assertEqual(r_cb.unexpected, [])

        # 3. Hydrogen mode: HA, HB1..HB3 are missing; optional H, H2, HXT are not missing
        rep_h = ala_no_oxt.check_missing_atoms(db, include_hydrogens=True)
        self.assertFalse(rep_h.is_complete)
        r_h = rep_h.residues[0]
        self.assertEqual(r_h.status, "incomplete")
        self.assertEqual(r_h.expected, 9)  # N, CA, C, O, CB, HA, HB1, HB2, HB3
        self.assertEqual(r_h.present, 5)
        self.assertEqual(r_h.missing, ["HA", "HB1", "HB2", "HB3"])
        self.assertEqual(r_h.unexpected, [])

        # Optional hydrogens (H1, H3) present are not reported as unexpected
        ala_nterm = build_residue_from_template("ALA", True, exclude=("H", "HXT"))
        ala_nterm.set_atom("90_H1", rs_br.Atom(symbol="H", name="H1", xyz=rs_br.Position(0.0, 0.0, 0.0)))
        ala_nterm.set_atom("91_H3", rs_br.Atom(symbol="H", name="H3", xyz=rs_br.Position(0.0, 0.0, 0.0)))
        rep_nterm = ala_nterm.check_missing_atoms(db, include_hydrogens=True)
        self.assertTrue(rep_nterm.is_complete)
        self.assertEqual(rep_nterm.residues[0].unexpected, [])

    def test_hoh_heavy_and_hydrogens(self):
        """HOH heavy atom mode vs hydrogen mode."""
        hoh = build_residue_from_template("HOH", False)

        # Heavy-atom mode: O is present -> Complete
        rep_heavy = hoh.check_missing_atoms()
        self.assertTrue(rep_heavy.is_complete)
        self.assertEqual(len(rep_heavy.residues), 1)
        r_heavy = rep_heavy.residues[0]
        self.assertEqual(r_heavy.status, "complete")
        self.assertEqual(r_heavy.expected, 1)
        self.assertEqual(r_heavy.present, 1)
        self.assertEqual(r_heavy.missing, [])

        # Hydrogen mode: H1 and H2 are missing -> Incomplete
        rep_h = hoh.check_missing_atoms(include_hydrogens=True)
        self.assertFalse(rep_h.is_complete)
        r_h = rep_h.residues[0]
        self.assertEqual(r_h.status, "incomplete")
        self.assertEqual(r_h.expected, 3)
        self.assertEqual(r_h.present, 1)
        self.assertEqual(r_h.missing, ["H1", "H2"])

    def test_unknown_component_and_custom_ccd_db(self):
        """Unknown component without template and with minimal custom CCD DB."""
        lig = rs_br.AtomGroup("LIG")
        lig.set_atom("1_C1", rs_br.Atom(symbol="C", name="C1", xyz=rs_br.Position(0.0, 0.0, 0.0)))
        lig.set_atom("2_O1", rs_br.Atom(symbol="O", name="O1", xyz=rs_br.Position(1.2, 0.0, 0.0)))

        # Without template in global DB -> NoTemplate
        rep_unk = lig.check_missing_atoms()
        self.assertTrue(rep_unk.is_complete)  # incomplete residues is empty
        self.assertEqual(rep_unk.total_missing, 0)
        self.assertEqual(len(rep_unk.no_template), 1)
        self.assertEqual(rep_unk.no_template_components, {"LIG": 1})

        r_unk = rep_unk.residues[0]
        self.assertEqual(r_unk.status, "no_template")
        self.assertEqual(r_unk.expected, 0)
        self.assertEqual(r_unk.present, 0)
        self.assertEqual(r_unk.missing, [])
        self.assertEqual(r_unk.unexpected, [])

        # Add minimal CCD mmCIF file via add_from_file to a custom DB
        minimal_cif = """\
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
"""
        with tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False) as f:
            f.write(minimal_cif)
            tmp_cif_path = f.name
        self.tmp_files.append(tmp_cif_path)

        custom_db = rs_br.CcdTemplateDb.builtin()
        count = custom_db.add_from_file(tmp_cif_path)
        self.assertEqual(count, 1)

        rep_custom = lig.check_missing_atoms(db=custom_db)
        self.assertFalse(rep_custom.is_complete)
        self.assertEqual(len(rep_custom.no_template), 0)
        self.assertEqual(rep_custom.no_template_components, {})

        r_custom = rep_custom.residues[0]
        self.assertEqual(r_custom.status, "incomplete")
        self.assertEqual(r_custom.expected, 3)  # C1, O1, C2
        self.assertEqual(r_custom.present, 2)
        self.assertEqual(r_custom.missing, ["C2"])

    def test_da_missing_n9_and_five_prime_terminus(self):
        """DA missing N9 and 5'-terminal phosphate detection."""
        # 1. DA missing N9
        da_no_n9 = build_residue_from_template("DA", False, exclude=("OP3", "N9"))
        rep_n9 = da_no_n9.check_missing_atoms()
        self.assertFalse(rep_n9.is_complete)
        r_n9 = rep_n9.residues[0]
        self.assertEqual(r_n9.status, "incomplete")
        self.assertEqual(r_n9.missing, ["N9"])

        # 2. 5'-terminal DA without any of P, OP1, OP2, OP3 -> Complete
        da_5term_no_phos = build_residue_from_template(
            "DA", False, exclude=("P", "OP1", "OP2", "OP3")
        )
        rep_5term = da_5term_no_phos.check_missing_atoms()
        self.assertTrue(rep_5term.is_complete)
        r_5term = rep_5term.residues[0]
        self.assertEqual(r_5term.status, "complete")
        self.assertEqual(r_5term.expected, 18)
        self.assertEqual(r_5term.present, 18)
        self.assertEqual(r_5term.missing, [])

        # 3. 5'-terminal DA missing only P (OP1, OP2 present) -> P is reported missing
        da_5term_no_p = build_residue_from_template("DA", False, exclude=("P", "OP3"))
        rep_no_p = da_5term_no_p.check_missing_atoms()
        self.assertFalse(rep_no_p.is_complete)
        r_no_p = rep_no_p.residues[0]
        self.assertEqual(r_no_p.status, "incomplete")
        self.assertEqual(r_no_p.missing, ["P"])

        # 4. Internal residue lacking P, OP1, OP2, OP3 is incomplete
        chain = rs_br.AtomGroup("A")
        chain.set_group("1", da_5term_no_phos)
        chain.set_group("2", da_5term_no_phos)
        rep_chain = chain.check_missing_atoms()
        self.assertEqual(len(rep_chain.residues), 2)
        self.assertEqual(rep_chain.residues[0].status, "complete")
        self.assertEqual(rep_chain.residues[1].status, "incomplete")
        self.assertEqual(rep_chain.residues[1].missing, ["P", "OP1", "OP2"])

    def test_unexpected_and_duplicate_atoms(self):
        """Unexpected atoms and duplicate atom names handling."""
        ala = build_residue_from_template("ALA", False, exclude=("OXT",))
        # Optional atom OXT present is not unexpected
        ala.set_atom("80_OXT", rs_br.Atom(symbol="O", name="OXT", xyz=rs_br.Position(0.6, 0.4, -1.7)))
        # Non-standard atoms ZN and FOO
        ala.set_atom("91_ZN", rs_br.Atom(symbol="Zn", name="ZN", xyz=rs_br.Position(1.0, 1.0, 1.0)))
        ala.set_atom("92_FOO", rs_br.Atom(symbol="C", name="FOO", xyz=rs_br.Position(2.0, 2.0, 2.0)))

        rep = ala.check_missing_atoms()
        self.assertTrue(rep.is_complete)
        r = rep.residues[0]
        self.assertEqual(r.status, "complete")
        self.assertEqual(r.expected, 5)
        self.assertEqual(r.present, 5)
        self.assertEqual(r.missing, [])
        self.assertEqual(r.unexpected, ["ZN", "FOO"])

        # Duplicate atom names: CA duplicated, unexpected XX duplicated
        ala_dup = build_residue_from_template("ALA", False, exclude=("OXT",))
        ala_dup.set_atom("90_CA", rs_br.Atom(symbol="C", name="CA", xyz=rs_br.Position(0.25, 0.41, 0.69)))
        ala_dup.set_atom("91_XX", rs_br.Atom(symbol="C", name="XX", xyz=rs_br.Position(1.0, 1.0, 1.0)))
        ala_dup.set_atom("92_YY", rs_br.Atom(symbol="O", name="YY", xyz=rs_br.Position(2.0, 2.0, 2.0)))
        ala_dup.set_atom("93_XX", rs_br.Atom(symbol="C", name="XX", xyz=rs_br.Position(3.0, 3.0, 3.0)))

        rep_dup = ala_dup.check_missing_atoms()
        self.assertTrue(rep_dup.is_complete)
        r_dup = rep_dup.residues[0]
        self.assertEqual(r_dup.expected, 5)
        self.assertEqual(r_dup.present, 5)
        self.assertEqual(r_dup.missing, [])
        self.assertEqual(r_dup.unexpected, ["XX", "YY"])

    def test_cys_disulfide_hg_optional(self):
        """CYS disulfide bond detection across models, chains, and subtrees."""
        def make_cys_with_sg(x: float) -> rs_br.AtomGroup:
            cys = build_residue_from_template(
                "CYS", True, exclude=("OXT", "H", "H2", "HXT", "HG")
            )
            # Find and set SG position
            for key, atom in cys.atoms():
                if atom.name == "SG":
                    cys.set_atom(key, rs_br.Atom(symbol="S", name="SG", xyz=rs_br.Position(x, 0.0, 0.0)))
                    break
            return cys

        # 1. Close CYS (SG-SG 2.05 A) in same model -> HG is optional
        model_close = rs_br.AtomGroup("model_1")
        chain_a = rs_br.AtomGroup("A")
        chain_a.set_group("1", make_cys_with_sg(0.0))
        chain_b = rs_br.AtomGroup("B")
        chain_b.set_group("1", make_cys_with_sg(2.05))
        model_close.set_group("A", chain_a)
        model_close.set_group("B", chain_b)

        rep_close = model_close.check_missing_atoms(include_hydrogens=True)
        self.assertEqual(len(rep_close.residues), 2)
        self.assertTrue(rep_close.is_complete)
        self.assertEqual(rep_close.total_missing, 0)
        for r in rep_close.residues:
            self.assertEqual(r.status, "complete")
            self.assertEqual(r.missing, [])

        # 2. Distant CYS (SG-SG 5.0 A) in same chain -> HG is missing on both
        chain_far = rs_br.AtomGroup("A")
        chain_far.set_group("1", make_cys_with_sg(0.0))
        chain_far.set_group("2", make_cys_with_sg(5.0))
        rep_far = chain_far.check_missing_atoms(include_hydrogens=True)
        self.assertFalse(rep_far.is_complete)
        self.assertEqual(rep_far.total_missing, 2)
        for r in rep_far.residues:
            self.assertEqual(r.status, "incomplete")
            self.assertEqual(r.missing, ["HG"])

        # 3. Nearby CYS (2.05 A) across different models -> not disulfide bonded
        root_models = rs_br.AtomGroup()
        m1 = rs_br.AtomGroup("model_1")
        m1_c = rs_br.AtomGroup("A")
        m1_c.set_group("1", make_cys_with_sg(0.0))
        m1.set_group("A", m1_c)

        m2 = rs_br.AtomGroup("model_2")
        m2_c = rs_br.AtomGroup("A")
        m2_c.set_group("1", make_cys_with_sg(2.05))
        m2.set_group("A", m2_c)

        root_models.set_group("model_1", m1)
        root_models.set_group("model_2", m2)

        rep_models = root_models.check_missing_atoms(include_hydrogens=True)
        self.assertEqual(len(rep_models.residues), 2)
        self.assertFalse(rep_models.is_complete)
        self.assertEqual(rep_models.total_missing, 2)
        for r in rep_models.residues:
            self.assertEqual(r.status, "incomplete")
            self.assertEqual(r.missing, ["HG"])

        # 4. Calling on single residue alone cannot see partner -> HG is missing
        single_cys = model_close.get_group("A").get_group("1")
        rep_single = single_cys.check_missing_atoms(include_hydrogens=True)
        self.assertEqual(len(rep_single.residues), 1)
        self.assertFalse(rep_single.is_complete)
        self.assertEqual(rep_single.residues[0].status, "incomplete")
        self.assertEqual(rep_single.residues[0].missing, ["HG"])

    def test_empty_tree_and_hierarchy_and_immutability(self):
        """Empty tree, hierarchy navigation, and tree immutability guarantee."""
        # Empty tree
        empty = rs_br.AtomGroup()
        rep_empty = empty.check_missing_atoms()
        self.assertEqual(len(rep_empty.residues), 0)
        self.assertTrue(rep_empty.is_complete)
        self.assertEqual(rep_empty.total_missing, 0)
        self.assertEqual(len(rep_empty.incomplete), 0)
        self.assertEqual(len(rep_empty.no_template), 0)
        self.assertEqual(rep_empty.no_template_components, {})

        # Hierarchy: root -> model_1 -> A -> "10", "2", "1"
        root = rs_br.AtomGroup()
        model = rs_br.AtomGroup("model_1")
        chain = rs_br.AtomGroup("A")
        chain.set_group("10", build_residue_from_template("ALA", False, exclude=("CB", "OXT")))
        chain.set_group("2", build_residue_from_template("HOH", False))
        chain.set_group("1", build_residue_from_template("ALA", False, exclude=("OXT",)))
        model.set_group("A", chain)
        root.set_group("model_1", model)

        raw_before = root.get_raw_data()
        rep_root = root.check_missing_atoms()
        self.assertEqual(
            root.get_raw_data(),
            raw_before,
            "AtomGroup must remain unmodified after check_missing_atoms",
        )

        # Check sort_nicely path traversal order
        paths = [r.path for r in rep_root.residues]
        self.assertEqual(paths, ["/model_1/A/1/", "/model_1/A/2/", "/model_1/A/10/"])
        self.assertEqual(rep_root.total_missing, 1)
        self.assertEqual(len(rep_root.incomplete), 1)
        self.assertEqual(rep_root.incomplete[0].path, "/model_1/A/10/")

        # Callable from model level
        model_ref = root.get_group("model_1")
        rep_model = model_ref.check_missing_atoms()
        self.assertEqual([r.path for r in rep_model.residues], paths)
        self.assertEqual(rep_model.total_missing, 1)

        # Callable from chain level
        chain_ref = model_ref.get_group("A")
        rep_chain = chain_ref.check_missing_atoms()
        self.assertEqual([r.path for r in rep_chain.residues], paths)
        self.assertEqual(rep_chain.total_missing, 1)

        # Callable from single-residue level
        res1 = chain_ref.get_group("1")
        rep_res1 = res1.check_missing_atoms()
        self.assertEqual(len(rep_res1.residues), 1)
        self.assertEqual(rep_res1.residues[0].path, "/model_1/A/1/")
        self.assertEqual(rep_res1.residues[0].status, "complete")

    def test_to_dict_str_and_repr(self):
        """to_dict JSON serialization, __str__, and __repr__ formatting."""
        chain = rs_br.AtomGroup("A")
        chain.set_group("1", build_residue_from_template("ALA", False, exclude=("CB", "OXT")))
        unk = rs_br.AtomGroup("UNK")
        unk.set_atom("1_X", rs_br.Atom(symbol="C", name="X", xyz=rs_br.Position(0.0, 0.0, 0.0)))
        chain.set_group("2", unk)
        rep = chain.check_missing_atoms()

        # to_dict produces valid JSON
        d = rep.to_dict()
        self.assertIsInstance(d, dict)
        json_str = json.dumps(d)
        loaded = json.loads(json_str)

        self.assertEqual(loaded["total_missing"], 1)
        self.assertFalse(loaded["is_complete"])
        self.assertEqual(loaded["no_template_components"], {"UNK": 1})
        self.assertEqual(len(loaded["residues"]), 2)
        self.assertEqual(len(loaded["incomplete"]), 1)
        self.assertEqual(len(loaded["no_template"]), 1)
        self.assertEqual(loaded["incomplete"][0]["missing"], ["CB"])

        # ResidueAtomCheck.to_dict()
        r0 = rep.residues[0]
        r0_dict = r0.to_dict()
        self.assertEqual(r0_dict["comp_id"], "ALA")
        self.assertEqual(r0_dict["status"], "incomplete")
        self.assertEqual(r0_dict["missing"], ["CB"])

        # __str__ inspection
        s = str(rep)
        self.assertIn("MissingAtomsReport", s)
        self.assertIn("/1/", s)
        self.assertIn("ALA", s)
        self.assertIn("CB", s)
        self.assertIn("UNK", s)

        # Complete __str__
        complete_rep = build_residue_from_template("ALA", False, exclude=("OXT",)).check_missing_atoms()
        s_complete = str(complete_rep)
        self.assertIn("complete (0 missing atoms)", s_complete)

        # __repr__ inspection
        repr_str = repr(rep)
        self.assertIn("MissingAtomsReport(complete=False, total_missing=1", repr_str)
        self.assertIn("ResidueAtomCheck", repr(r0))

    def test_real_fixtures_2fb4_and_2mgo(self):
        """Fixed expectation validation and timing benchmark on real 2FB4 and 2MGO fixtures."""
        path_2fb4 = os.path.join(DATA_DIR, "2FB4.cif")
        path_2mgo = os.path.join(DATA_DIR, "2MGO.cif")

        cif_2fb4 = rs_br.SimpleMmcif(path_2fb4)
        ag_2fb4 = cif_2fb4.get_structure_atomgroup()

        # 2FB4: heavy-atom mode
        raw_2fb4_before = ag_2fb4.get_raw_data()
        rep_2fb4 = ag_2fb4.check_missing_atoms(include_hydrogens=False)
        self.assertEqual(ag_2fb4.get_raw_data(), raw_2fb4_before)
        self.assertEqual(len(rep_2fb4.residues), 550)
        self.assertTrue(rep_2fb4.is_complete)
        self.assertEqual(len(rep_2fb4.incomplete), 0)
        self.assertEqual(rep_2fb4.total_missing, 0)
        self.assertEqual(rep_2fb4.no_template_components, {})

        # 2FB4: hydrogen mode
        rep_2fb4_h = ag_2fb4.check_missing_atoms(include_hydrogens=True)
        self.assertEqual(len(rep_2fb4_h.residues), 550)
        self.assertFalse(rep_2fb4_h.is_complete)
        self.assertEqual(len(rep_2fb4_h.incomplete), 550)
        self.assertEqual(rep_2fb4_h.total_missing, 3031)
        self.assertEqual(rep_2fb4_h.no_template_components, {})

        # 2MGO: heavy-atom mode
        cif_2mgo = rs_br.SimpleMmcif(path_2mgo)
        ag_2mgo = cif_2mgo.get_structure_atomgroup()
        raw_2mgo_before = ag_2mgo.get_raw_data()
        rep_2mgo = ag_2mgo.check_missing_atoms(include_hydrogens=False)
        self.assertEqual(ag_2mgo.get_raw_data(), raw_2mgo_before)
        self.assertEqual(len(rep_2mgo.residues), 180)
        self.assertTrue(rep_2mgo.is_complete)
        self.assertEqual(len(rep_2mgo.incomplete), 0)
        self.assertEqual(rep_2mgo.total_missing, 0)
        self.assertEqual(rep_2mgo.no_template_components, {})

        # 2MGO: hydrogen mode (NMR structure with hydrogens)
        rep_2mgo_h = ag_2mgo.check_missing_atoms(include_hydrogens=True)
        self.assertEqual(len(rep_2mgo_h.residues), 180)
        self.assertTrue(rep_2mgo_h.is_complete)
        self.assertEqual(len(rep_2mgo_h.incomplete), 0)
        self.assertEqual(rep_2mgo_h.total_missing, 0)
        self.assertEqual(rep_2mgo_h.no_template_components, {})

        # Benchmark: 2FB4 check_missing_atoms() average time over 20 runs
        runs = 20
        total_time = timeit.timeit(
            lambda: ag_2fb4.check_missing_atoms(include_hydrogens=False),
            number=runs,
        )
        avg_ms = (total_time / runs) * 1000.0
        print(f"\n[BENCHMARK] 2FB4 check_missing_atoms(heavy) 20-run average: {avg_ms:.2f} ms")


if __name__ == "__main__":
    unittest.main()
