#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.bond import Bond
from proteindf_bridge.atomgroup import AtomGroup
from proteindf_bridge.atom import Atom
from proteindf_bridge.position import Position


class TestBond(unittest.TestCase):
    def test_bond_setup(self):
        ag = AtomGroup("mol")
        # 共有結合距離（〜1.5Å）にある2つの炭素原子
        c1 = Atom(symbol="C", xyz=Position([0.0, 0.0, 0.0]))
        c2 = Atom(symbol="C", xyz=Position([1.5, 0.0, 0.0]))
        # 遠く離れた炭素原子（10.0Å）
        c3 = Atom(symbol="C", xyz=Position([10.0, 0.0, 0.0]))

        ag.set_atom("1", c1)
        ag.set_atom("2", c2)
        ag.set_atom("3", c3)

        bond = Bond()
        bond.setup(ag)

        self.assertIsNotNone(bond._distmat)
        self.assertIsNotNone(bond._bondmat)
        # c1 と c2 の距離は約 1.5
        self.assertAlmostEqual(bond._distmat.get(1, 0), 1.5)
        # c1 と c2 は結合あり(1)、c1/c2 と c3 は結合なし(0)
        self.assertEqual(bond._bondmat.get(1, 0), 1)
        self.assertEqual(bond._bondmat.get(2, 0), 0)
        self.assertEqual(bond._bondmat.get(2, 1), 0)

    def test_covalent_criteria_and_no_duplicate(self):
        # 1. from proteindf_bridge import Bond ができる
        from proteindf_bridge import Bond as ExportedBond

        self.assertIs(ExportedBond, Bond)

        # 2. C-C 3.0 Å (旧基準では結合 1.7+1.7+0.4=3.8 Å) が結合にならない
        ag_cc = AtomGroup("mol_cc")
        c1 = Atom(symbol="C", xyz=Position([0.0, 0.0, 0.0]))
        c2 = Atom(symbol="C", xyz=Position([3.0, 0.0, 0.0]))
        ag_cc.set_atom("1", c1)
        ag_cc.set_atom("2", c2)

        bond = Bond()
        bond.setup(ag_cc)
        self.assertEqual(bond._bondmat.get(1, 0), 0)
        self.assertEqual(len(ag_cc.get_bond_list()), 0)

        # 3. C-H 1.1 Å (cov_C(0.76) + cov_H(0.31) + 0.45 = 1.52 Å) が結合になる
        ag_ch = AtomGroup("mol_ch")
        c = Atom(symbol="C", xyz=Position([0.0, 0.0, 0.0]))
        h = Atom(symbol="H", xyz=Position([1.1, 0.0, 0.0]))
        ag_ch.set_atom("1", c)
        ag_ch.set_atom("2", h)

        bond = Bond()
        bond.setup(ag_ch)
        self.assertEqual(bond._bondmat.get(1, 0), 1)
        self.assertEqual(len(ag_ch.get_bond_list()), 1)
        self.assertEqual(ag_ch.get_bond_list()[0][2], 1)

        # 4. 既存の結合が重複しない
        # 既存結合(例えば次数2)がある状態でsetupを実行しても二重登録されず、既存の次数が保持される
        ag_dup = AtomGroup("mol_dup")
        c_a = Atom(symbol="C", xyz=Position([0.0, 0.0, 0.0]))
        c_b = Atom(symbol="C", xyz=Position([1.5, 0.0, 0.0]))
        ag_dup.set_atom("1", c_a)
        ag_dup.set_atom("2", c_b)
        ag_dup.add_bond(ag_dup["1"], ag_dup["2"], order=2)
        self.assertEqual(len(ag_dup.get_bond_list()), 1)
        self.assertEqual(ag_dup.get_bond_list()[0][2], 2)

        bond = Bond()
        bond.setup(ag_dup)
        self.assertEqual(len(ag_dup.get_bond_list()), 1)
        self.assertEqual(ag_dup.get_bond_list()[0][2], 2)

    def test_setup_1hls_pdb_bond_count_matches_rust(self):
        import os
        from proteindf_bridge.biopdb import Pdb as PyPdb
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
                "1hls.pdb",
            )
        )
        self.assertTrue(os.path.exists(pdb_path), f"File not found: {pdb_path}")

        # 純Python版 Pdb で読み込み、Bond().setup
        py_pdb = PyPdb(pdb_path)
        py_ag = py_pdb.get_atomgroup()
        py_bond = Bond()
        py_bond.setup(py_ag)
        py_bonds = py_ag.get_bond_list()

        # Rust版 (proteindf_bridge.rs) の Pdb と Bond
        rs_pdb = rs_br.Pdb(pdb_path)
        rs_ag = rs_pdb.get_atomgroup()
        rs_bond = rs_br.Bond()
        rs_bond.setup(rs_ag)
        rs_bonds = rs_ag.get_bond_list()

        self.assertEqual(len(py_bonds), len(rs_bonds))

    def test_bond_cross_model_avoided(self):
        # 1モデル単独の構造: 2つの炭素原子が共有結合距離
        single_m = AtomGroup("model_1")
        c1 = Atom(symbol="C", xyz=Position([0.0, 0.0, 0.0]))
        c2 = Atom(symbol="C", xyz=Position([1.5, 0.0, 0.0]))
        single_m.set_atom("C1", c1)
        single_m.set_atom("C2", c2)
        single_root = AtomGroup("mol_single")
        single_root.set_group("model_1", single_m)
        single_bond = Bond()
        single_bond.setup(single_root)
        single_bonds = single_root.get_bond_list()
        self.assertEqual(len(single_bonds), 1)

        # 2モデルを同じ座標に重ねた構造
        root = AtomGroup("mol")
        m1 = AtomGroup("model_1")
        m2 = AtomGroup("model_2")

        c1_1 = Atom(symbol="C", xyz=Position([0.0, 0.0, 0.0]))
        c2_1 = Atom(symbol="C", xyz=Position([1.5, 0.0, 0.0]))
        m1.set_atom("C1", c1_1)
        m1.set_atom("C2", c2_1)

        c1_2 = Atom(symbol="C", xyz=Position([0.0, 0.0, 0.0]))
        c2_2 = Atom(symbol="C", xyz=Position([1.5, 0.0, 0.0]))
        m2.set_atom("C1", c1_2)
        m2.set_atom("C2", c2_2)

        root.set_group("model_1", m1)
        root.set_group("model_2", m2)

        bond = Bond()
        bond.setup(root)

        bonds = root.get_bond_list()
        self.assertEqual(len(bonds), len(single_bonds) * 2)

        # 各結合が同一モデル内の原子どうしであることを確認
        from proteindf_bridge.bond import _extract_model

        for b in bonds:
            m_a = _extract_model(b[0])
            m_b = _extract_model(b[1])
            self.assertIsNotNone(m_a)
            self.assertEqual(m_a, m_b)

        # 密行列 (bondmat) でも別モデル間が 0 であることを確認
        self.assertIsNotNone(bond._bondmat)
        for p in range(len(bond._atoms)):
            m_p = _extract_model(bond._atoms[p].path)
            for q in range(p):
                m_q = _extract_model(bond._atoms[q].path)
                if m_p is not None and m_q is not None and m_p != m_q:
                    self.assertEqual(bond._bondmat.get(p, q), 0)

    def test_setup_2mgo_pdb_cross_model_avoided(self):
        import os
        from proteindf_bridge.biopdb import Pdb as PyPdb
        import proteindf_bridge.rs as rs_br
        from proteindf_bridge.bond import _extract_model

        pdb_path = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "rust",
                "crates",
                "proteindf-bridge",
                "tests",
                "data",
                "2MGO.pdb",
            )
        )
        self.assertTrue(os.path.exists(pdb_path), f"File not found: {pdb_path}")

        # 純Python版 Pdb で読み込み、Bond().setup
        py_pdb = PyPdb(pdb_path)
        py_ag = py_pdb.get_atomgroup()
        py_bond = Bond()
        py_bond.setup(py_ag)
        py_bonds = py_ag.get_bond_list()

        # 純Python版: 結合が2720本で、別モデル間の結合が0本
        self.assertEqual(len(py_bonds), 2720)
        py_cross_bonds = [
            b for b in py_bonds if _extract_model(b[0]) != _extract_model(b[1])
        ]
        self.assertEqual(len(py_cross_bonds), 0)

        # Rust版 (proteindf_bridge.rs) の Pdb と Bond
        rs_pdb = rs_br.Pdb(pdb_path)
        rs_ag = rs_pdb.get_atomgroup()
        rs_bond = rs_br.Bond()
        rs_bond.setup(rs_ag)
        rs_bonds = rs_ag.get_bond_list()

        # Rust版: 結合が2720本で、別モデル間の結合が0本
        self.assertEqual(len(rs_bonds), 2720)
        rs_cross_bonds = [
            b for b in rs_bonds if _extract_model(b[0]) != _extract_model(b[1])
        ]
        self.assertEqual(len(rs_cross_bonds), 0)


if __name__ == "__main__":
    unittest.main()
