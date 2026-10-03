#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Test suite for PR#46:
Validation of PyO3 bindings for Hydrogenation:
1. AtomGroup.add_missing_hydrogens(db=None) on 1hls.pdb and mmCIF roundtrip
2. OverallHydrogenationReport and per-residue HydrogenationReport inspection
3. Resilience to crystal waters (skipped_residues)
4. Detection of missing sidechain templates on modified residues (step_errors)
5. Non-destructive behavior on subtree copies (in-place modification guarantee)
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


def count_all_hydrogens(group: rs_br.AtomGroup) -> int:
    """Recursively counts all hydrogen atoms in an AtomGroup."""
    direct = sum(
        1
        for _key, atom in group.atoms()
        if atom.atomic_number == 1 or atom.symbol == "H"
    )
    sub = sum(count_all_hydrogens(subgroup) for _key, subgroup in group.groups())
    return direct + sub


class TestRsHydrogenation(unittest.TestCase):
    """Tests for Phase 10 / PR#46: Hydrogenation bindings and reports."""

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

    def _load_1hls(self) -> rs_br.AtomGroup:
        pdb_path = os.path.join(DATA_DIR, "1hls.pdb")
        pdb = rs_br.Pdb()
        pdb.load(pdb_path)
        return pdb.get_atomgroup()

    def test_hydrogenation_pipeline_roundtrip_1hls(self):
        """
        PR#46 Criterion 1:
        1hls.pdb load -> setup() -> add_missing_hydrogens() -> SimpleMmcif.save_structure() -> reload
        roundtrip verification against Rust baseline criteria.

        Referenced Rust tests:
          test_mmcif_writer.rs::test_hydrogenation_pipeline_roundtrip (lines 275-309)
          test_orchestrator.rs::test_hydrogenate_partial_existing_hydrogens (lines 122-143)
          test_orchestrator.rs::test_hydrogenate_1hls_full_structure (lines 45-119)

        Baseline criteria verified from Rust tests:
          1. 1hls.pdb natively has 379 hydrogens (lacking terminal HE2/HD1/etc in 15 residues).
          2. add_missing_hydrogens() adds exactly 15 missing hydrogens (394 - 379 = 15).
          3. Total added = 15, total removed = 0.
          4. 15 residues modified (e.g. GLU A4 receives HE2).
          5. Final total hydrogen count across all 51 residues = 394.
          6. Saving to mmCIF and reloading retains exact hydrogen count 394 and all atoms.
        """
        structure = self._load_1hls()

        # Baseline: 1hls.pdb has 379 native hydrogens
        init_h = count_all_hydrogens(structure)
        self.assertEqual(
            init_h,
            379,
            "1hls.pdb must have exactly 379 native hydrogens (test_orchestrator.rs line 132)",
        )

        # Bond setup
        structure.setup()

        # Hydrogenation with default built-in DB
        report = structure.add_missing_hydrogens()

        # Verify report values against test_orchestrator.rs lines 136-143
        self.assertEqual(report.total_added_hydrogens, 15)
        self.assertEqual(report.total_removed_hydrogens, 0)
        self.assertEqual(
            report.hydrogenated_residues,
            15,
            "15 residues were missing 1 hydrogen each in 1hls.pdb natively",
        )
        self.assertEqual(len(report.skipped_residues), 0)
        self.assertEqual(len(report.step_errors), 0)

        # Final structure hydrogen count = 394
        final_h = count_all_hydrogens(structure)
        self.assertEqual(
            final_h,
            394,
            "1hls must have 394 total hydrogens after hydrogenation (test_orchestrator.rs line 142)",
        )

        # Verify residue addition: GLU A4 received HE2
        self.assertIn("/model_1/A/4/", report.residue_reports)
        res_a4_rep = report.residue_reports["/model_1/A/4/"]
        self.assertEqual(res_a4_rep.added_hydrogens, 1)
        self.assertEqual(res_a4_rep.added_atom_names, ["HE2"])

        # Save to mmCIF and reload (test_mmcif_writer.rs lines 294-309)
        tmp_cif = self._get_tmp_path(".cif")
        rs_br.SimpleMmcif.save_structure(structure, tmp_cif)

        reloaded_cif = rs_br.SimpleMmcif(tmp_cif)
        reloaded_struct = reloaded_cif.get_structure_atomgroup()

        reloaded_h = count_all_hydrogens(reloaded_struct)
        self.assertEqual(
            reloaded_h,
            394,
            "Reloaded mmCIF structure must preserve all 394 hydrogens (test_mmcif_writer.rs line 306)",
        )
        self.assertEqual(
            reloaded_struct.get_number_of_all_atoms(),
            structure.get_number_of_all_atoms(),
            "Reloaded mmCIF structure must have identical total atom count",
        )

    def test_hydrogenate_water_hoh_resilient_skipped_residues(self):
        """
        PR#46 Criterion 2:
        Resilience against components with insufficient heavy atoms (such as HOH water).
        Water must be safely skipped and recorded in skipped_residues, without aborting
        the hydrogenation of neighboring amino acids.

        Referenced Rust test:
          test_orchestrator.rs::test_hydrogenate_water_hoh_resilient (lines 150-208)

        Baseline criteria verified from Rust test:
          1. ALA residue (path /model_1/A/1/) is hydrogenated (hydrogenated_residues == 1).
          2. HOH water (path /model_1/A/2/) is recorded in skipped_residues (len == 1).
          3. Skipped reason explains heavy atom deficiency ("common heavy atoms" or "Sidechain/general hydrogenation").
          4. HOH is NOT in residue_reports (no double-counting).
        """
        chain = rs_br.AtomGroup(name="A")
        chain.path = "/model_1/A/"

        # Residue 1: ALA (N, CA, C, O, CB)
        res_ala = rs_br.AtomGroup(name="ALA")
        res_ala.path = "/model_1/A/1/"

        n = rs_br.Atom("N")
        n.name = "N"
        n.xyz = [0.0, 0.0, 0.0]
        res_ala["N"] = n

        ca = rs_br.Atom("C")
        ca.name = "CA"
        ca.xyz = [1.46, 0.0, 0.0]
        res_ala["CA"] = ca

        c = rs_br.Atom("C")
        c.name = "C"
        c.xyz = [2.0, 1.4, 0.0]
        res_ala["C"] = c

        o = rs_br.Atom("O")
        o.name = "O"
        o.xyz = [3.2, 1.5, 0.0]
        res_ala["O"] = o

        cb = rs_br.Atom("C")
        cb.name = "CB"
        cb.xyz = [2.0, -0.7, 1.2]
        res_ala["CB"] = cb

        chain["1"] = res_ala

        # Residue 2: HOH (single 'O' heavy atom)
        res_hoh = rs_br.AtomGroup(name="HOH")
        res_hoh.path = "/model_1/A/2/"
        hoh_o = rs_br.Atom("O")
        hoh_o.name = "O"
        hoh_o.xyz = [10.0, 10.0, 10.0]
        res_hoh["O"] = hoh_o

        chain["2"] = res_hoh

        report = chain.add_missing_hydrogens()

        # ALA must be hydrogenated
        res1 = chain.get_group("1")
        self.assertIsNotNone(res1)
        self.assertTrue(res1.has_atom("H1"))
        self.assertTrue(res1.has_atom("HA"))
        self.assertEqual(report.hydrogenated_residues, 1)
        self.assertIn("/model_1/A/1/", report.residue_reports)

        # HOH must be recorded in skipped_residues (not in residue_reports)
        self.assertEqual(len(report.skipped_residues), 1)
        skipped_path, reason = report.skipped_residues[0]
        self.assertEqual(skipped_path, "/model_1/A/2/")
        self.assertTrue(
            "common heavy atoms" in reason or "Sidechain/general hydrogenation" in reason,
            f"Reason should explain heavy atom deficiency: {reason}",
        )
        self.assertNotIn("/model_1/A/2/", report.residue_reports)

    def test_hydrogenate_unknown_residue_step_errors(self):
        """
        PR#46 Criterion 2:
        Detection of missing sidechain templates on partially modified residues (step_errors).
        When an unknown residue UNK with backbone N/CA/C has backbone hydrogens added,
        it must be recorded in residue_reports, NOT in skipped_residues, and its missing
        sidechain template must be explicitly recorded in step_errors.

        Referenced Rust test:
          test_orchestrator.rs::test_partially_modified_unknown_residue_recorded_in_reports_only (lines 325-406)

        Baseline criteria verified from Rust test:
          1. UNK residue receives backbone H (is modified).
          2. UNK (/model_1/A/2/) appears in residue_reports.
          3. UNK does NOT appear in skipped_residues.
          4. UNK appears in step_errors with message containing "No CCD template found for residue 'UNK'".
        """
        chain = rs_br.AtomGroup(name="A")
        chain.path = "/model_1/A/"

        # Residue 1: ALA
        res_ala = rs_br.AtomGroup(name="ALA")
        res_ala.path = "/model_1/A/1/"

        n = rs_br.Atom("N")
        n.name = "N"
        n.xyz = [0.0, 0.0, 0.0]
        res_ala["N"] = n

        ca = rs_br.Atom("C")
        ca.name = "CA"
        ca.xyz = [1.46, 0.0, 0.0]
        res_ala["CA"] = ca

        c = rs_br.Atom("C")
        c.name = "C"
        c.xyz = [2.0, 1.4, 0.0]
        res_ala["C"] = c

        o = rs_br.Atom("O")
        o.name = "O"
        o.xyz = [3.2, 1.5, 0.0]
        res_ala["O"] = o

        cb = rs_br.Atom("C")
        cb.name = "CB"
        cb.xyz = [2.0, -0.7, 1.2]
        res_ala["CB"] = cb

        chain["1"] = res_ala

        # Residue 2: UNK (has N, CA, C so backbone H is added, but no CCD template for sidechain)
        res_unk = rs_br.AtomGroup(name="UNK")
        res_unk.path = "/model_1/A/2/"

        unk_n = rs_br.Atom("N")
        unk_n.name = "N"
        unk_n.xyz = [1.3, 2.5, 0.0]
        res_unk["N"] = unk_n

        unk_ca = rs_br.Atom("C")
        unk_ca.name = "CA"
        unk_ca.xyz = [1.8, 3.8, 0.0]
        res_unk["CA"] = unk_ca

        unk_c = rs_br.Atom("C")
        unk_c.name = "C"
        unk_c.xyz = [3.3, 3.9, 0.0]
        res_unk["C"] = unk_c

        chain["2"] = res_unk

        report = chain.add_missing_hydrogens()

        # UNK must have backbone H added
        res_unk_after = chain.get_group("2")
        self.assertIsNotNone(res_unk_after)
        self.assertTrue(res_unk_after.has_atom("H"))

        # UNK must be in residue_reports, NOT in skipped_residues
        self.assertIn("/model_1/A/2/", report.residue_reports)
        self.assertFalse(
            any(p == "/model_1/A/2/" for p, _ in report.skipped_residues),
            "Modified UNK residue must NEVER appear in skipped_residues",
        )

        # UNK missing template must be in step_errors
        unk_errors = [err for p, err in report.step_errors if p == "/model_1/A/2/"]
        self.assertEqual(
            len(unk_errors),
            1,
            f"Missing CCD template on modified UNK must be recorded in step_errors: {report.step_errors}",
        )
        self.assertIn("No CCD template found for residue 'UNK'", unk_errors[0])

    def test_hydrogenation_reports_identity_and_properties(self):
        """
        PR#45 Review Revision Rule:
        Report collections (residue_reports, skipped_residues, step_errors) must return
        the identical Python object on repeated accesses without reallocating copies.
        Also verifies all read-only properties and repr() of reports.
        """
        chain = rs_br.AtomGroup(name="A")
        chain.path = "/model_1/A/"

        res_ala = rs_br.AtomGroup(name="ALA")
        res_ala.path = "/model_1/A/1/"

        n = rs_br.Atom("N")
        n.name = "N"
        n.xyz = [0.0, 0.0, 0.0]
        res_ala["N"] = n

        ca = rs_br.Atom("C")
        ca.name = "CA"
        ca.xyz = [1.46, 0.0, 0.0]
        res_ala["CA"] = ca

        c = rs_br.Atom("C")
        c.name = "C"
        c.xyz = [2.0, 1.4, 0.0]
        res_ala["C"] = c

        o = rs_br.Atom("O")
        o.name = "O"
        o.xyz = [3.2, 1.5, 0.0]
        res_ala["O"] = o

        cb = rs_br.Atom("C")
        cb.name = "CB"
        cb.xyz = [2.0, -0.7, 1.2]
        res_ala["CB"] = cb

        chain["1"] = res_ala

        res_hoh = rs_br.AtomGroup(name="HOH")
        res_hoh.path = "/model_1/A/2/"
        hoh_o = rs_br.Atom("O")
        hoh_o.name = "O"
        hoh_o.xyz = [10.0, 10.0, 10.0]
        res_hoh["O"] = hoh_o

        chain["2"] = res_hoh

        report = chain.add_missing_hydrogens()

        # 1. Identity on repeated access
        self.assertIs(
            report.residue_reports,
            report.residue_reports,
            "report.residue_reports must return identical dict instance across accesses",
        )
        self.assertIs(
            report.skipped_residues,
            report.skipped_residues,
            "report.skipped_residues must return identical list instance across accesses",
        )
        self.assertIs(
            report.step_errors,
            report.step_errors,
            "report.step_errors must return identical list instance across accesses",
        )

        # 2. Overall repr
        rep_str = repr(report)
        self.assertIn("OverallHydrogenationReport(", rep_str)
        self.assertIn("added=", rep_str)
        self.assertIn("skipped=1", rep_str)

        # 3. Per-residue report inspection
        res_rep = report.residue_reports["/model_1/A/1/"]
        self.assertGreater(res_rep.added_hydrogens, 0)
        self.assertIn("H1", res_rep.added_atom_names)
        self.assertEqual(res_rep.removed_hydrogens, 0)
        self.assertEqual(res_rep.removed_atom_names, [])

        res_rep_str = repr(res_rep)
        self.assertIn("HydrogenationReport(added=", res_rep_str)

    def test_subtree_copy_independence(self):
        """
        RUST_PORT_SPEC.md §4.3 Common Design Policy 2:
        Subtree copies returned by get_group() are independent; calling add_missing_hydrogens()
        on a subtree copy must not modify the parent tree.
        """
        structure = self._load_1hls()
        structure.setup()

        init_h = count_all_hydrogens(structure)
        self.assertEqual(init_h, 379)

        # Obtain subtree copy
        model_copy = structure.get_group("model_1")
        if model_copy is None:
            model_copy = structure.get_group("1")
        self.assertIsNotNone(model_copy)

        # Hydrogenate the subtree copy
        sub_rep = model_copy.add_missing_hydrogens()
        self.assertEqual(sub_rep.total_added_hydrogens, 15)

        # Subtree copy now has 394 hydrogens
        self.assertEqual(count_all_hydrogens(model_copy), 394)

        # Original parent tree must remain unchanged at 379 hydrogens
        self.assertEqual(
            count_all_hydrogens(structure),
            379,
            "Modifying a subtree copy must not alter the original parent AtomGroup",
        )

    def test_add_missing_hydrogens_with_explicit_db(self):
        """
        Test calling add_missing_hydrogens with an explicit CcdTemplateDb instance.
        """
        structure = self._load_1hls()
        structure.setup()

        db = rs_br.CcdTemplateDb.builtin()
        report = structure.add_missing_hydrogens(db=db)
        self.assertEqual(report.total_added_hydrogens, 15)
        self.assertEqual(count_all_hydrogens(structure), 394)


if __name__ == "__main__":
    unittest.main()
