import os
import tempfile
import unittest

import proteindf_bridge.rs as rs


class TestRsAgDirectStrCsv(unittest.TestCase):
    def _read_csv(self, ag):
        with tempfile.NamedTemporaryFile(mode="r", delete=False) as f:
            path = f.name
        try:
            ag.save_csv(path)
            with open(path, "rb") as f:
                return f.read().decode("utf-8")
        finally:
            if os.path.exists(path):
                os.remove(path)

    def test_empty_tree(self):
        """Verify str() and save_csv() on an empty AtomGroup."""
        ag = rs.AtomGroup("empty")
        expected_str = "# group key= name=empty\n"
        expected_csv = ""

        self.assertEqual(str(ag), expected_str)
        self.assertEqual(self._read_csv(ag), expected_csv)

    def test_single_atom_tree(self):
        """Verify str() and save_csv() on an AtomGroup with a single atom."""
        ag = rs.AtomGroup("single")
        atom = rs.Atom("C", name="C1", xyz=[1.0, 2.0, 3.0], charge=-0.5)
        ag["C1"] = atom

        expected_str = (
            "# group key= name=single\n"
            "C (C1  )    1.000    2.000    3.000, -0.50,    0.000    0.000    0.000 /C1\n"
        )
        expected_csv = "C,1.0,2.0,3.0,-0.5\r\n"

        self.assertEqual(str(ag), expected_str)
        self.assertEqual(self._read_csv(ag), expected_csv)

    def test_nested_tree_and_subtrees(self):
        """Verify str() and save_csv() on nested tree with bonds, secondary structure, dummy atom, and charge."""
        root = rs.AtomGroup("complex_root")
        sub1 = rs.AtomGroup("chain_A")
        sub1.secondary_structure = "H"
        sub2 = rs.AtomGroup("res_1")

        a1 = rs.Atom("C", name="CA", xyz=[1.5, 2.5, 3.5], charge=0.1)
        a2 = rs.Atom("C", name="CB", xyz=[2.0, 3.0, 4.0], charge=-0.2)
        dummy = rs.Atom(0, name="DUM", xyz=[0.0, 0.0, 0.0], charge=0.0)

        sub2["CA"] = a1
        sub2["CB"] = a2
        sub2["DUM"] = dummy
        sub1["sub2"] = sub2
        root["sub1"] = sub1

        res1 = root["sub1"]["sub2"]
        ca = res1["CA"]
        cb = res1["CB"]
        res1.add_bond(ca, cb, 1)

        # 1. Root group
        expected_root_str = (
            "# group key= name=complex_root\n"
            "  # group key=sub1 name=chain_A\n"
            "    # group key=sub2 name=res_1\n"
            "    C (CA  )    1.500    2.500    3.500,  0.10,    0.000    0.000    0.000 /sub1/sub2/CA\n"
            "    C (CB  )    2.000    3.000    4.000, -0.20,    0.000    0.000    0.000 /sub1/sub2/CB\n"
            "    X (DUM )    0.000    0.000    0.000,  0.00,    0.000    0.000    0.000 /sub1/sub2/DUM\n"
            "    bond CA <-1-> CB\n"
        )
        expected_root_csv = (
            "sub1,chain_A,sub2,res_1,C,1.5,2.5,3.5,0.1\r\n"
            "sub1,chain_A,sub2,res_1,C,2.0,3.0,4.0,-0.2\r\n"
            "sub1,chain_A,sub2,res_1,X,0.0,0.0,0.0,0.0\r\n"
        )
        self.assertEqual(str(root), expected_root_str)
        self.assertEqual(self._read_csv(root), expected_root_csv)

        # 2. Non-root subtree (sub1)
        sub1_in_tree = root["sub1"]
        expected_sub1_str = (
            "# group key= name=chain_A\n"
            "  # group key=sub2 name=res_1\n"
            "  C (CA  )    1.500    2.500    3.500,  0.10,    0.000    0.000    0.000 /sub1/sub2/CA\n"
            "  C (CB  )    2.000    3.000    4.000, -0.20,    0.000    0.000    0.000 /sub1/sub2/CB\n"
            "  X (DUM )    0.000    0.000    0.000,  0.00,    0.000    0.000    0.000 /sub1/sub2/DUM\n"
            "  bond CA <-1-> CB\n"
        )
        expected_sub1_csv = (
            "sub2,res_1,C,1.5,2.5,3.5,0.1\r\n"
            "sub2,res_1,C,2.0,3.0,4.0,-0.2\r\n"
            "sub2,res_1,X,0.0,0.0,0.0,0.0\r\n"
        )
        self.assertEqual(str(sub1_in_tree), expected_sub1_str)
        self.assertEqual(self._read_csv(sub1_in_tree), expected_sub1_csv)

        # 3. Non-root subtree (sub2 / res1)
        expected_sub2_str = (
            "# group key= name=res_1\n"
            "C (CA  )    1.500    2.500    3.500,  0.10,    0.000    0.000    0.000 /sub1/sub2/CA\n"
            "C (CB  )    2.000    3.000    4.000, -0.20,    0.000    0.000    0.000 /sub1/sub2/CB\n"
            "X (DUM )    0.000    0.000    0.000,  0.00,    0.000    0.000    0.000 /sub1/sub2/DUM\n"
            "bond CA <-1-> CB\n"
        )
        expected_sub2_csv = (
            "C,1.5,2.5,3.5,0.1\r\n"
            "C,2.0,3.0,4.0,-0.2\r\n"
            "X,0.0,0.0,0.0,0.0\r\n"
        )
        self.assertEqual(str(res1), expected_sub2_str)
        self.assertEqual(self._read_csv(res1), expected_sub2_csv)


if __name__ == "__main__":
    unittest.main()
