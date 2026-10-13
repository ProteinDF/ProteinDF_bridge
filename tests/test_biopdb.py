import os
import unittest
import doctest

from proteindf_bridge.biopdb import Pdb

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "proteindf_bridge", "data"))


class PdbTests(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_init(self):
        pdb = Pdb(os.path.join(DATA_DIR, "2MGO.pdb"))

    def test_get_atomgroup(self):
        pdb = Pdb(os.path.join(DATA_DIR, "2MGO.pdb"))
        ag = pdb.get_atomgroup()

        # model
        self.assertEqual(ag.get_number_of_groups(), 20)
        # chain
        self.assertEqual(ag["model_1"].get_number_of_groups(), 1)
        # residues
        self.assertEqual(ag["model_1"]["A"].get_number_of_groups(), 9)

        self.assertEqual(ag["model_1"]["A"][1].get_number_of_atoms(), 12)  # CYS1
        self.assertEqual(ag["model_1"]["A"][1].name, "CYS")
        self.assertEqual(ag["model_1"]["A"][2].get_number_of_atoms(), 21)  # TYR2
        self.assertEqual(ag["model_1"]["A"][2].name, "TYR")
        self.assertEqual(ag["model_1"]["A"][3].get_number_of_atoms(), 19)  # ILE3
        self.assertEqual(ag["model_1"]["A"][3].name, "ILE")
        self.assertEqual(ag["model_1"]["A"][4].get_number_of_atoms(), 17)  # GLN4
        self.assertEqual(ag["model_1"]["A"][4].name, "GLN")
        self.assertEqual(ag["model_1"]["A"][5].get_number_of_atoms(), 14)  # ASN5
        self.assertEqual(ag["model_1"]["A"][5].name, "ASN")
        self.assertEqual(ag["model_1"]["A"][6].get_number_of_atoms(), 10)  # CYS6
        self.assertEqual(ag["model_1"]["A"][6].name, "CYS")
        self.assertEqual(ag["model_1"]["A"][7].get_number_of_atoms(), 14)  # PRO7
        self.assertEqual(ag["model_1"]["A"][7].name, "PRO")
        self.assertEqual(ag["model_1"]["A"][8].get_number_of_atoms(), 19)  # LEU8
        self.assertEqual(ag["model_1"]["A"][8].name, "LEU")
        self.assertEqual(ag["model_1"]["A"][9].get_number_of_atoms(),  8)  # GLY9
        self.assertEqual(ag["model_1"]["A"][9].name, "GLY")

    def test_set_by_atomgroup(self):
        pass

    def test_split_residue_key(self):
        from proteindf_bridge import split_residue_key

        self.assertEqual(split_residue_key("101"), (101, ""))
        self.assertEqual(split_residue_key("101A"), (101, "A"))
        self.assertEqual(split_residue_key("-1"), (-1, ""))
        self.assertEqual(split_residue_key("-1B"), (-1, "B"))

        for invalid in ["", "ABC", "101AB", " 101", "101 ", "A101"]:
            with self.assertRaises(ValueError):
                split_residue_key(invalid)

    def test_select_residues_with_icode(self):
        from proteindf_bridge import AtomGroup, Atom, Position, Modeling

        chain = AtomGroup()
        chain.name = "H"

        for key, resname in [("100", "ALA"), ("101", "GLY"), ("101A", "HIS"), ("101B", "GLY"), ("102", "VAL")]:
            res = AtomGroup()
            res.name = resname
            res.set_atom("CA", Atom(name="CA", symbol="C", xyz=Position(0.0, 0.0, 0.0)))
            chain.set_group(key, res)

        modeling = Modeling()
        selected = modeling.select_residues(chain, 101, 101)

        # In Python AtomGroup, residues matching range were merged into selected AtomGroup.
        # 101, 101A, 101B all have sequence number 101, so their CA atoms were merged.
        self.assertEqual(len(selected.get_atom_list()), 1)
        self.assertTrue(selected.has_atom("CA"))

        # Range 101 to 102 should select 101, 101A, 101B, and 102
        selected_101_102 = modeling.select_residues(chain, 101, 102)
        self.assertEqual(len(selected_101_102.get_atom_list()), 1)

        # Range 100 to 100 selects only 100
        selected_100 = modeling.select_residues(chain, 100, 100)
        self.assertEqual(len(selected_100.get_atom_list()), 1)

        # Range outside (e.g. 103 to 105) selects none
        selected_none = modeling.select_residues(chain, 103, 105)
        self.assertEqual(len(selected_none.get_atom_list()), 0)

    def test_2fb4_pdb_icode(self):
        import proteindf_bridge.rs as rs_br

        pdb_path = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "rust",
                "crates",
                "proteindf-bridge",
                "tests",
                "data",
                "2FB4.pdb",
            )
        )
        self.assertTrue(os.path.exists(pdb_path), f"2FB4.pdb not found at {pdb_path}")

        # 1. Read with pure Python Pdb
        py_pdb = Pdb(pdb_path)
        py_ag = py_pdb.get_atomgroup()

        # Check H chain residues 101 to 101D
        chain_h = py_ag["model_1"]["H"]
        for key, expected_name in [
            ("101", "GLY"),
            ("101A", "HIS"),
            ("101B", "GLY"),
            ("101C", "PHE"),
            ("101D", "CYS"),
        ]:
            self.assertTrue(chain_h.has_group(key), f"Chain H missing residue {key}")
            self.assertEqual(chain_h[key].name, expected_name)

        # 2. Check that the set of residue keys across all models and chains matches Rust Pdb
        rs_pdb = rs_br.Pdb(pdb_path)
        rs_ag = rs_pdb.get_atomgroup()

        def collect_residue_keys(ag):
            keys = set()
            for model_key, model in ag.groups():
                for chain_key, chain in model.groups():
                    for res_key, _ in chain.groups():
                        keys.add((model_key, chain_key, res_key))
            return keys

        py_keys = collect_residue_keys(py_ag)
        rs_keys = collect_residue_keys(rs_ag)
        self.assertEqual(py_keys, rs_keys)

        # 3. Write out (set_by_atomgroup) and re-read, verifying icode residues are preserved
        out_pdb = Pdb()
        out_pdb.set_by_atomgroup(py_ag)
        pdb_text = str(out_pdb)

        # Save to temp file and reload with Pdb
        import tempfile
        with tempfile.NamedTemporaryFile("w+", suffix=".pdb", delete=False) as f:
            f.write(pdb_text)
            temp_path = f.name

        try:
            reloaded_pdb = Pdb(temp_path)
            reloaded_ag = reloaded_pdb.get_atomgroup()

            reloaded_chain_h = reloaded_ag["model_1"]["H"]
            for key, expected_name in [
                ("101", "GLY"),
                ("101A", "HIS"),
                ("101B", "GLY"),
                ("101C", "PHE"),
                ("101D", "CYS"),
            ]:
                self.assertTrue(reloaded_chain_h.has_group(key), f"Reloaded Chain H missing residue {key}")
                self.assertEqual(reloaded_chain_h[key].name, expected_name)

            reloaded_keys = collect_residue_keys(reloaded_ag)
            self.assertEqual(reloaded_keys, py_keys)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)


def load_tests(loader, tests, ignore):
    from proteindf_bridge import biopdb
    tests.addTests(doctest.DocTestSuite(biopdb))
    return tests


if __name__ == '__main__':
    unittest.main()
