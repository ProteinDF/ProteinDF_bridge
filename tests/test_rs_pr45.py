#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Test suite for PR#45:
Validation of PyO3 bindings for Foundation and Bond Resolution:
1. AtomGroup.setup() and AtomGroup.setup_with_db(db)
2. CcdTemplateDb class (builtin, add_from_file, lookup, merge)
3. Schema validation (validate_schema, is_model_level, is_chain_level, is_residue_level)
4. Secondary structure field reading and writing on AtomGroup
5. mmCIF get_structure_atomgroup_with_report (unresolved _struct_conn reporting)
"""

import os
import tempfile
import unittest

import proteindf_bridge_rs as rs_br

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


class TestRsPr45BondResolutionAndFoundation(unittest.TestCase):
    """Tests for PR#45: Foundation & Bond Resolution bindings."""

    def setUp(self):
        self.tmp_files = []

    def tearDown(self):
        for f in self.tmp_files:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError:
                    pass

    def _get_tmp_path(self, suffix=".cif"):
        fd, path = tempfile.mkstemp(suffix=suffix)
        os.close(fd)
        self.tmp_files.append(path)
        return path

    def test_setup_1hls_real_pdb(self):
        """
        PR#45 Criterion 1: AtomGroup.setup() on 1hls.pdb.
        Referenced Rust test:
          test_ccd_templates.rs::test_atomgroup_setup_1hls_real_pdb (lines 710-802)
        Baseline criteria verified from Rust test:
          1. Before setup(), ag.get_bond_list() is empty (len == 0).
          2. GLU 4 C=O bond has order 2 and appears exactly once.
          3. ARG 22 in Chain B has at least 2 double bonds (C=O and CZ=NH2).
          4. Inter-residue peptide bond GLU 4 - PRO 5 has order 1.
          5. No duplicate bond pairs across all resolved bonds.
        """
        pdb_path = os.path.join(DATA_DIR, "1hls.pdb")
        pdb = rs_br.Pdb()
        pdb.load(pdb_path)
        ag = pdb.get_atomgroup()

        # 1. Before setup, no bonds exist in raw PDB loader output
        self.assertEqual(len(ag.get_bond_list()), 0)

        # Execute setup()
        ag.setup()
        final_bonds = ag.get_bond_list()
        self.assertGreater(len(final_bonds), 0)

        # 2. GLU 4 C=O bond: order == 2, exactly 1 occurrence
        # (Reference: test_ccd_templates.rs lines 728-750)
        glu4_co_bonds = [
            b
            for b in final_bonds
            if (
                ("/A/4/" in b[0] and (b[0].endswith("_C") or b[0].endswith("/C")))
                and ("/A/4/" in b[1] and (b[1].endswith("_O") or b[1].endswith("/O")))
            )
            or (
                ("/A/4/" in b[0] and (b[0].endswith("_O") or b[0].endswith("/O")))
                and ("/A/4/" in b[1] and (b[1].endswith("_C") or b[1].endswith("/C")))
            )
        ]
        self.assertEqual(len(glu4_co_bonds), 1, "GLU 4 C=O bond should appear exactly once")
        self.assertEqual(glu4_co_bonds[0][2], 2, "GLU 4 C=O bond should have order 2 from CCD template")

        # 3. ARG 22 in Chain B has multiple double bonds
        # (Reference: test_ccd_templates.rs lines 753-763)
        arg_double_bonds = [
            b
            for b in final_bonds
            if ("/B/22/" in b[0] or "/B/22/" in b[1]) and b[2] == 2
        ]
        self.assertGreaterEqual(
            len(arg_double_bonds), 2, "ARG 22 should have at least 2 double bonds"
        )

        # 4. Inter-residue peptide bond between GLU 4 and PRO 5
        # (Reference: test_ccd_templates.rs lines 766-785)
        peptide_bonds = [
            b
            for b in final_bonds
            if (
                ("/A/4/" in b[0] and (b[0].endswith("_C") or b[0].endswith("/C")))
                and ("/A/5/" in b[1] and (b[1].endswith("_N") or b[1].endswith("/N")))
            )
            or (
                ("/A/5/" in b[0] and (b[0].endswith("_N") or b[0].endswith("/N")))
                and ("/A/4/" in b[1] and (b[1].endswith("_C") or b[1].endswith("/C")))
            )
        ]
        self.assertEqual(len(peptide_bonds), 1, "Peptide bond GLU 4 - PRO 5 must be found")
        self.assertEqual(peptide_bonds[0][2], 1, "Peptide bond order must be 1")

        # 5. No duplicates across all final bonds
        # (Reference: test_ccd_templates.rs lines 788-800)
        seen = set()
        for b in final_bonds:
            key = (b[0], b[1]) if b[0] <= b[1] else (b[1], b[0])
            self.assertNotIn(key, seen, f"Duplicate bond found: {key}")
            seen.add(key)

    def test_setup_coexistence_with_file_bonds(self):
        """
        PR#45 Criterion 1: Coexistence of file-derived inter-residue bonds with ag.setup().
        Referenced Rust test:
          test_mmcif_writer.rs::test_setup_coexistence_with_file_bonds (lines 1127-1185)
          test_mmcif_writer.rs::test_roundtrip_inter_residue_bonds_1wct (lines 1187-1198)
        Baseline criteria verified from Rust test:
          1. 1WCT.cif has 10 initial inter-residue bonds.
          2. Total bond count strictly increases after setup() (intra-residue bonds added).
          3. All initial file-derived inter-residue bonds are preserved with identical orders.
          4. No duplicate bonds between any atom pair.
        """
        test_cases = ["1WCT.cif", "1HLS.cif", "2FB4.cif"]

        def _is_inter_residue(b):
            parts1 = [p for p in b[0].split("/") if p]
            parts2 = [p for p in b[1].split("/") if p]
            # depth 4: [model, chain, res, atom]
            if len(parts1) >= 3 and len(parts2) >= 3:
                return (parts1[0], parts1[1], parts1[2]) != (parts2[0], parts2[1], parts2[2])
            return False

        for filename in test_cases:
            cif_path = os.path.join(DATA_DIR, filename)
            cif = rs_br.SimpleMmcif(cif_path)
            ag = cif.get_structure_atomgroup(select_model=1)

            init_bonds = ag.get_bond_list()
            init_inter_bonds = [b for b in init_bonds if _is_inter_residue(b)]
            init_total_count = len(init_bonds)

            if filename == "1WCT.cif":
                # Baseline: 10 inter-residue bonds in 1WCT (2 disulf + 8 covale)
                # (Reference: test_mmcif_writer.rs lines 1196-1198)
                self.assertEqual(len(init_inter_bonds), 10)

            # Perform setup()
            ag.setup()

            post_bonds = ag.get_bond_list()
            post_inter_bonds = [b for b in post_bonds if _is_inter_residue(b)]
            post_total_count = len(post_bonds)

            # 1. Total bond count increases
            self.assertGreater(
                post_total_count,
                init_total_count,
                f"{filename}: setup() should add intra-residue bonds",
            )

            # 2. File-derived inter-residue bonds preserved
            norm_post_inter = {
                (min(b[0], b[1]), max(b[0], b[1])): b[2] for b in post_inter_bonds
            }
            for init_b in init_inter_bonds:
                key = (min(init_b[0], init_b[1]), max(init_b[0], init_b[1]))
                self.assertIn(key, norm_post_inter, f"{filename}: missing bond {key}")
                self.assertEqual(
                    norm_post_inter[key],
                    init_b[2],
                    f"{filename}: bond order mismatch for {key}",
                )

            # 3. No duplicate bonds
            seen = set()
            for b in post_bonds:
                key = (min(b[0], b[1]), max(b[0], b[1]))
                self.assertNotIn(key, seen, f"{filename}: duplicate bond found {key}")
                seen.add(key)

    def test_user_supplied_ccd_db_and_setup_with_db(self):
        """
        PR#45 Criterion 2: Load user-supplied CCD file (ALA.cif) and use with setup_with_db().
        Referenced Rust test:
          test_ccd_templates.rs::test_from_mmcif_block_ala_cif_matches_embedded (lines 320-353)
          test_ccd_templates.rs::test_resolve_bonds_preserves_existing_file_bonds (lines 804-839)
        Baseline criteria verified from Rust test:
          1. Builtin DB contains 29 standard templates.
          2. Empty DB starts with 0 templates.
          3. Loading ALA.cif adds 1 template with 13 atoms and 12 bonds.
          4. ALA template contains canonical C=O double bond (order 2).
          5. setup_with_db() resolves intra-residue bonds from loaded DB.
        """
        # 1. Builtin DB baseline
        builtin_db = rs_br.CcdTemplateDb.builtin()
        self.assertEqual(len(builtin_db), 29)
        self.assertFalse(builtin_db.is_empty())
        self.assertIn("ALA", builtin_db)
        self.assertIn("GLY", builtin_db)

        # 2. Empty DB
        empty_db = rs_br.CcdTemplateDb.empty()
        self.assertEqual(len(empty_db), 0)
        self.assertTrue(empty_db.is_empty())
        self.assertNotIn("ALA", empty_db)

        # 3. Load user-supplied CCD file (ALA.cif)
        # (Reference: test_ccd_templates.rs lines 324-329)
        ala_cif_path = os.path.join(DATA_DIR, "ALA.cif")
        count = empty_db.add_from_file(ala_cif_path)
        self.assertEqual(count, 1)
        self.assertEqual(len(empty_db), 1)
        self.assertIn("ALA", empty_db)

        # Lookup ALA template
        ala_tmpl = empty_db.lookup("ALA")
        self.assertIsNotNone(ala_tmpl)
        self.assertEqual(ala_tmpl.comp_id, "ALA")
        # Baseline: ALA has 13 atoms and 12 bonds in CCD template
        # (Reference: test_ccd_templates.rs lines 336-347)
        self.assertEqual(len(ala_tmpl.atoms), 13)
        self.assertEqual(len(ala_tmpl.bonds), 12)

        ca_atom = ala_tmpl.get_atom("CA")
        self.assertIsNotNone(ca_atom)
        self.assertEqual(ca_atom.name, "CA")
        self.assertEqual(ca_atom.element, "C")
        self.assertIsNotNone(ca_atom.ideal_xyz)
        self.assertFalse(ca_atom.is_hydrogen)

        # 4. Use setup_with_db on a synthetic ALA residue
        # (Reference: test_ccd_templates.rs lines 806-856)
        ag = rs_br.AtomGroup(name="ALA")
        ag.path = "/model_1/A/1/"

        n = rs_br.Atom("N")
        n.xyz = [0.0, 0.0, 0.0]
        n.name = "N"
        ca = rs_br.Atom("C")
        ca.xyz = [1.46, 0.0, 0.0]
        ca.name = "CA"
        c = rs_br.Atom("C")
        c.xyz = [2.0, 1.4, 0.0]
        c.name = "C"
        o = rs_br.Atom("O")
        o.xyz = [1.3, 2.4, 0.0]
        o.name = "O"
        cb = rs_br.Atom("C")
        cb.xyz = [2.0, -0.7, 1.2]
        cb.name = "CB"

        ag["N"] = n
        ag["CA"] = ca
        ag["C"] = c
        ag["O"] = o
        ag["CB"] = cb

        self.assertEqual(len(ag.get_bond_list()), 0)
        ag.setup_with_db(empty_db)

        bonds = ag.get_bond_list()
        # Baseline: 4 heavy atom bonds for ALA (N-CA, CA-C, C-O, CA-CB)
        # (Reference: test_ccd_templates.rs line 855)
        self.assertEqual(len(bonds), 4)

        co_bonds = [
            b
            for b in bonds
            if (b[0].endswith("/C") and b[1].endswith("/O"))
            or (b[0].endswith("/O") and b[1].endswith("/C"))
        ]
        self.assertEqual(len(co_bonds), 1)
        self.assertEqual(co_bonds[0][2], 2, "C=O bond should have order 2 from CCD template")

    def test_subtree_copy_setup_does_not_mutate_original(self):
        """
        PR#45 Criterion 3 & RUST_PORT_SPEC.md §4.3「共通の設計方針」2:
        Calling setup() on a subtree copy mutates the copy in place, but does NOT affect the original tree.
        """
        pdb_path = os.path.join(DATA_DIR, "1hls.pdb")
        pdb = rs_br.Pdb()
        pdb.load(pdb_path)
        original_ag = pdb.get_atomgroup()

        # Original has 0 bonds before setup()
        self.assertEqual(len(original_ag.get_bond_list()), 0)

        # Obtain a subtree copy via __getitem__
        subtree_copy = original_ag["model_1"]["A"]

        # Call setup() on the subtree copy
        subtree_copy.setup()

        # Subtree copy now has bonds resolved
        subtree_bonds = subtree_copy.get_bond_list()
        self.assertGreater(len(subtree_bonds), 0)

        # Original tree MUST NOT be mutated (still 0 bonds)
        original_bonds = original_ag.get_bond_list()
        self.assertEqual(
            len(original_bonds),
            0,
            "Original tree must remain unaffected when setup() is invoked on a subtree copy",
        )

    def test_validate_schema_and_hierarchy_levels(self):
        """
        PR#45 Criterion 4: Schema validation and hierarchy level predicates.
        Referenced Rust test:
          test_schema.rs::test_schema_level_checks_normal_hierarchy (lines 18-62)
          test_schema.rs::test_schema_violation_direct_atoms_in_chain (lines 64-96)
          test_schema.rs::test_schema_violation_subgroup_in_residue_and_depth (lines 98-142)
          test_schema.rs::test_schema_real_pdb_fixture_1hls (lines 271-285)
        Baseline criteria verified from Rust test:
          1. Normal hierarchy has path depths 0..3 and correct is_*_level flags.
          2. Normal hierarchy and real 1hls.pdb produce 0 schema violations.
          3. Direct atoms in chain (depth 2) trigger DirectAtomsAtNonResidueLevel.
          4. Subgroups inside residue (depth 3) trigger SubgroupsInResidue,
             ExcessiveDepth (depth 4), and DirectAtomsAtNonResidueLevel (depth 4).
        """
        # 1. Normal hierarchy level checks
        root = rs_br.AtomGroup(name="")
        root.path = "/"
        model = rs_br.AtomGroup(name="model_1")
        model.path = "/model_1/"
        chain = rs_br.AtomGroup(name="A")
        chain.path = "/model_1/A/"
        residue = rs_br.AtomGroup(name="1")
        residue.path = "/model_1/A/1/"
        atom = rs_br.Atom("C")
        atom.xyz = [0.0, 0.0, 0.0]
        atom.name = "CA"

        residue["CA"] = atom
        chain["1"] = residue
        model["A"] = chain
        root["model_1"] = model

        # Root: depth 0
        self.assertEqual(root.path_depth(), 0)
        self.assertFalse(root.is_model_level())
        self.assertFalse(root.is_chain_level())
        self.assertFalse(root.is_residue_level())

        # Model: depth 1
        m = root["model_1"]
        self.assertEqual(m.path_depth(), 1)
        self.assertTrue(m.is_model_level())
        self.assertFalse(m.is_chain_level())
        self.assertFalse(m.is_residue_level())

        # Chain: depth 2
        c = m["A"]
        self.assertEqual(c.path_depth(), 2)
        self.assertFalse(c.is_model_level())
        self.assertTrue(c.is_chain_level())
        self.assertFalse(c.is_residue_level())

        # Residue: depth 3
        r = c["1"]
        self.assertEqual(r.path_depth(), 3)
        self.assertFalse(r.is_model_level())
        self.assertFalse(r.is_chain_level())
        self.assertTrue(r.is_residue_level())

        # 2. Zero violations on valid hierarchy
        violations = root.validate_schema()
        self.assertEqual(len(violations), 0)

        # 1hls.pdb fixture has zero violations
        # (Reference: test_schema.rs lines 271-285)
        pdb = rs_br.Pdb()
        pdb.load(os.path.join(DATA_DIR, "1hls.pdb"))
        ag_1hls = pdb.get_atomgroup()
        self.assertEqual(len(ag_1hls.validate_schema()), 0)

        # 3. Violation: direct atom in chain (depth 2)
        # (Reference: test_schema.rs lines 107-144)
        root_viol1 = rs_br.AtomGroup(name="")
        root_viol1.path = "/"
        model_viol1 = rs_br.AtomGroup(name="model_1")
        model_viol1.path = "/model_1/"
        chain_viol1 = rs_br.AtomGroup(name="A")
        chain_viol1.path = "/model_1/A/"
        res_viol1 = rs_br.AtomGroup(name="1")
        res_viol1.path = "/model_1/A/1/"
        ca1 = rs_br.Atom("C")
        ca1.name = "CA"
        res_viol1["CA"] = ca1
        chain_viol1["1"] = res_viol1

        het = rs_br.Atom("O")
        het.name = "HOH_100"
        chain_viol1["HOH_100"] = het

        model_viol1["A"] = chain_viol1
        root_viol1["model_1"] = model_viol1

        violations_chain = root_viol1.validate_schema()
        self.assertEqual(len(violations_chain), 1)
        v = violations_chain[0]
        self.assertEqual(v.violation_type, "DirectAtomsAtNonResidueLevel")
        self.assertEqual(v.path, "/model_1/A/")
        self.assertEqual(v.depth, 2)
        self.assertEqual(v.atom_keys, ["HOH_100"])
        self.assertIn("DirectAtomsAtNonResidueLevel", repr(v))
        self.assertIn("atoms are only allowed at residue level", str(v))

        # 4. Violation: subgroup inside residue (depth 3)
        # (Reference: test_schema.rs lines 147-190)
        root_viol2 = rs_br.AtomGroup(name="")
        root_viol2.path = "/"
        model_viol2 = rs_br.AtomGroup(name="model_1")
        model_viol2.path = "/model_1/"
        chain_viol2 = rs_br.AtomGroup(name="A")
        chain_viol2.path = "/model_1/A/"
        res_viol2 = rs_br.AtomGroup(name="6")
        res_viol2.path = "/model_1/A/6/"
        sub = rs_br.AtomGroup(name="sub")
        sub.path = "/model_1/A/6/sub/"
        sub_atom = rs_br.Atom("C")
        sub_atom.name = "C1"
        sub["C1"] = sub_atom
        res_viol2["sub"] = sub
        chain_viol2["6"] = res_viol2
        model_viol2["A"] = chain_viol2
        root_viol2["model_1"] = model_viol2

        violations_res = root_viol2.validate_schema()
        types = [viol.violation_type for viol in violations_res]
        self.assertIn("SubgroupsInResidue", types)
        self.assertIn("ExcessiveDepth", types)
        self.assertIn("DirectAtomsAtNonResidueLevel", types)

    def test_secondary_structure_field(self):
        """
        PR#45: Secondary structure field read and write on residue-level AtomGroup.
        Baseline criteria:
          1. Default secondary_structure is None.
          2. Accepts 'H', 'E', '-', and None.
          3. Rejects invalid code strings with BrValueError / ValueError.
        """
        ag = rs_br.AtomGroup(name="1")
        ag.path = "/model_1/A/1/"

        # Default is None
        self.assertIsNone(ag.secondary_structure)

        # Set Helix
        ag.secondary_structure = "H"
        self.assertEqual(ag.secondary_structure, "H")

        # Set Strand
        ag.secondary_structure = "E"
        self.assertEqual(ag.secondary_structure, "E")

        # Set Loop
        ag.secondary_structure = "-"
        self.assertEqual(ag.secondary_structure, "-")

        # Reset to None
        ag.secondary_structure = None
        self.assertIsNone(ag.secondary_structure)

        # Invalid assignment raises ValueError
        with self.assertRaises((ValueError, rs_br.BrValueError)):
            ag.secondary_structure = "X"

    def test_mmcif_with_report(self):
        """
        PR#45 Criterion 4 & 5: SimpleMmcif.get_structure_atomgroup_with_report().
        Referenced Rust tests:
          test_mmcif.rs::test_get_structure_atomgroup_with_report_synthetic_unresolved (lines 850-966)
          test_mmcif.rs::test_get_structure_atomgroup_with_report_real_pdb_fixtures (lines 1090-1148)
        Baseline criteria verified from Rust test:
          1. Real files 1HLS.cif, 2FB4.cif, 1WCT.cif produce 0 unresolved struct_conns.
          2. Synthetic CIF with altloc 'B' partner under altloc 'A' filter produces exactly
             1 unresolved struct_conn (disulf1, AtomNotFound on partner 2 SG).
          3. Under altloc 'B' filter, 0 unresolved struct_conns are reported.
        """
        # 1. Real fixtures: zero unresolved connections
        # (Reference: test_mmcif.rs lines 1095-1125)
        for filename in ["1HLS.cif", "2FB4.cif", "1WCT.cif"]:
            cif_path = os.path.join(DATA_DIR, filename)
            cif = rs_br.SimpleMmcif(cif_path)
            report = cif.get_structure_atomgroup_with_report()
            self.assertFalse(
                report.has_unresolved(),
                f"{filename} should have zero unresolved struct_conns",
            )
            self.assertEqual(len(report.unresolved_struct_conns), 0)
            self.assertGreater(report.atomgroup.get_number_of_all_atoms(), 0)

        # 2. Synthetic fixture with unresolved connection due to altloc
        # (Reference: test_mmcif.rs lines 851-965)
        cif_content = """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_asym_id
_atom_site.auth_comp_id
_atom_site.auth_seq_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . CYS A 1 ? 0.000 0.000 0.000 A CYS 1 N 1
ATOM 2 C CA . CYS A 1 ? 1.000 0.000 0.000 A CYS 1 CA 1
ATOM 3 C C . CYS A 1 ? 2.000 0.000 0.000 A CYS 1 C 1
ATOM 4 O O . CYS A 1 ? 3.000 0.000 0.000 A CYS 1 O 1
ATOM 5 C CB . CYS A 1 ? 1.000 1.000 0.000 A CYS 1 CB 1
ATOM 6 S SG . CYS A 1 ? 1.000 2.000 0.000 A CYS 1 SG 1
ATOM 7 N N . CYS A 2 ? 0.000 0.000 5.000 A CYS 2 N 1
ATOM 8 C CA . CYS A 2 ? 1.000 0.000 5.000 A CYS 2 CA 1
ATOM 9 C C . CYS A 2 ? 2.000 0.000 5.000 A CYS 2 C 1
ATOM 10 O O . CYS A 2 ? 3.000 0.000 5.000 A CYS 2 O 1
ATOM 11 C CB . CYS A 2 ? 1.000 1.000 5.000 A CYS 2 CB 1
ATOM 12 S SG B CYS A 2 ? 1.000 2.000 5.000 A CYS 2 SG 1

loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.pdbx_ptnr1_PDB_ins_code
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_ptnr2_PDB_ins_code
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.pdbx_value_order
disulf1 disulf A CYS 1 SG ? A CYS 1 A CYS 2 SG ? A CYS 2 sing
"""
        tmp_cif = self._get_tmp_path(".cif")
        with open(tmp_cif, "w", encoding="utf-8") as f:
            f.write(cif_content)

        cif_synth = rs_br.SimpleMmcif(tmp_cif)

        # Case A: Filter altloc "A" excludes CYS 2 SG (which is altloc "B")
        # (Reference: test_mmcif.rs lines 910-928)
        report_a = cif_synth.get_structure_atomgroup_with_report(select_altloc="A")
        self.assertTrue(report_a.has_unresolved())
        self.assertEqual(len(report_a.unresolved_struct_conns), 1)

        unres = report_a.unresolved_struct_conns[0]
        self.assertEqual(unres.conn_id, "disulf1")
        self.assertEqual(unres.conn_type_id, "disulf")
        self.assertEqual(unres.model_name, "model_1")
        self.assertIsNone(unres.ptnr1_unresolved)
        self.assertIsNotNone(unres.ptnr2_unresolved)
        self.assertEqual(unres.ptnr2_unresolved.reason, "AtomNotFound")
        self.assertEqual(unres.ptnr2_unresolved.chain_id, "A")
        self.assertEqual(unres.ptnr2_unresolved.res_key, "2")
        self.assertEqual(unres.ptnr2_unresolved.atom_name, "SG")
        self.assertIn("atom 'SG' in residue 'A/2' not found", str(unres.ptnr2_unresolved))
        self.assertIn("disulf1", repr(unres))

        # Case B: Filter altloc "B" retains CYS 2 SG
        # (Reference: test_mmcif.rs lines 955-965)
        report_b = cif_synth.get_structure_atomgroup_with_report(select_altloc="B")
        self.assertFalse(report_b.has_unresolved())
        self.assertEqual(len(report_b.unresolved_struct_conns), 0)


if __name__ == "__main__":
    unittest.main()
