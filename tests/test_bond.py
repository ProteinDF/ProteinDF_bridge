#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from proteindf_bridge.bond import Bond
from proteindf_bridge.atomgroup import AtomGroup
from proteindf_bridge.atom import Atom
from proteindf_bridge.position import Position
from proteindf_bridge.format import Format


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

        self.assertIsNotNone(bond.distmat)
        self.assertIsNotNone(bond.bondmat)
        # c1 と c2 の距離は約 1.5
        self.assertAlmostEqual(bond.distmat.get(1, 0), 1.5)
        # c1 と c2 は結合あり(1)、c1/c2 と c3 は結合なし(0)
        self.assertEqual(bond.bondmat.get(1, 0), 1)
        self.assertEqual(bond.bondmat.get(2, 0), 0)
        self.assertEqual(bond.bondmat.get(2, 1), 0)

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
        self.assertEqual(bond.bondmat.get(1, 0), 0)
        self.assertEqual(len(ag_cc.get_bond_list()), 0)

        # 3. C-H 1.1 Å (cov_C(0.76) + cov_H(0.31) + 0.45 = 1.52 Å) が結合になる
        ag_ch = AtomGroup("mol_ch")
        c = Atom(symbol="C", xyz=Position([0.0, 0.0, 0.0]))
        h = Atom(symbol="H", xyz=Position([1.1, 0.0, 0.0]))
        ag_ch.set_atom("1", c)
        ag_ch.set_atom("2", h)

        bond = Bond()
        bond.setup(ag_ch)
        self.assertEqual(bond.bondmat.get(1, 0), 1)
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

        self.assertEqual(len(py_bonds), 793)
        self.assertEqual(len(rs_bonds), 793)

    def test_bond_cross_model_avoided_with_custom_names(self):
        # (a) 子グループの名前が model_* でなくても、is_models を満たす2モデルを重ねると別モデル間に結合しない
        def make_test_protein(name, pos1, pos2):
            model = AtomGroup(name)
            chain = AtomGroup("A")
            res = AtomGroup("ALA_1")
            c1 = Atom(symbol="C", xyz=Position(pos1))
            c2 = Atom(symbol="C", xyz=Position(pos2))
            res.set_atom("C1", c1)
            res.set_atom("C2", c2)
            chain.set_group("ALA_1", res)
            model.set_group("A", chain)
            return model

        # 1モデル単独の構造
        single_m = make_test_protein("prot_A", [0.0, 0.0, 0.0], [1.5, 0.0, 0.0])
        single_root = AtomGroup("single_root")
        single_root.set_group("prot_A", single_m)
        self.assertTrue(Format.is_models(single_root))
        single_bond = Bond()
        single_bond.setup(single_root)
        single_bonds = single_root.get_bond_list()
        self.assertEqual(len(single_bonds), 1)

        # 2モデルを同じ座標に重ねた構造
        root = AtomGroup("root")
        m1 = make_test_protein("prot_A", [0.0, 0.0, 0.0], [1.5, 0.0, 0.0])
        m2 = make_test_protein("prot_B", [0.0, 0.0, 0.0], [1.5, 0.0, 0.0])
        root.set_group("prot_A", m1)
        root.set_group("prot_B", m2)
        self.assertTrue(Format.is_models(root))

        bond = Bond()
        bond.setup(root)

        bonds = root.get_bond_list()
        self.assertEqual(len(bonds), len(single_bonds) * 2)

        # 各結合が同一モデル内の原子どうしであることを確認
        def get_model_from_path(path):
            parts = [p for p in str(path).split("/") if p and p != "root"]
            return parts[0] if parts else None

        for b in bonds:
            m_a = get_model_from_path(b[0])
            m_b = get_model_from_path(b[1])
            self.assertIsNotNone(m_a)
            self.assertEqual(m_a, m_b)

        # 密行列 (bondmat) でも別モデル間が 0 であることを確認
        self.assertIsNotNone(bond.bondmat)
        atoms = root.get_atom_list()
        for p in range(len(atoms)):
            m_p = get_model_from_path(atoms[p].path)
            for q in range(p):
                m_q = get_model_from_path(atoms[q].path)
                if m_p != m_q:
                    self.assertEqual(bond.bondmat.get(p, q), 0)

    def test_bond_cross_group_allowed_when_not_models(self):
        # (b) is_models を満たさない木では、別の最上位グループの原子どうしでも距離が近ければ結合する
        root = AtomGroup("root")
        g1 = AtomGroup("g1")
        g2 = AtomGroup("g2")
        c1 = Atom(symbol="C", xyz=Position([0.0, 0.0, 0.0]))
        c2 = Atom(symbol="C", xyz=Position([1.5, 0.0, 0.0]))
        g1.set_atom("C1", c1)
        g2.set_atom("C2", c2)
        root.set_group("g1", g1)
        root.set_group("g2", g2)
        # ルート直下に原子を追加して is_models を確実に False にする
        root.set_atom("C_root", Atom(symbol="C", xyz=Position([100.0, 100.0, 100.0])))

        self.assertFalse(Format.is_models(root))

        bond = Bond()
        bond.setup(root)

        # g1 の C1 と g2 の C2 の距離は 1.5Å なので結合が形成される
        bonds = root.get_bond_list()
        cross_bonds = [
            b for b in bonds
            if ("g1" in b[0] and "g2" in b[1]) or ("g2" in b[0] and "g1" in b[1])
        ]
        self.assertGreaterEqual(len(cross_bonds), 1)

    def test_setup_2mgo_pdb_cross_model_avoided(self):
        # (c) 2MGOで両版とも2720本・別モデル間0本
        # (d) ag["model_3"]だけを渡すと136本
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
                "2MGO.pdb",
            )
        )
        self.assertTrue(os.path.exists(pdb_path), f"File not found: {pdb_path}")

        def get_model_name(path):
            parts = [p for p in str(path).split("/") if p]
            for part in parts:
                if part.startswith("model_"):
                    return part
            return None

        # 純Python版 Pdb で読み込み、Bond().setup
        py_pdb = PyPdb(pdb_path)
        py_ag = py_pdb.get_atomgroup()
        self.assertTrue(Format.is_models(py_ag))
        py_bond = Bond()
        py_bond.setup(py_ag)
        py_bonds = py_ag.get_bond_list()

        # 純Python版: 結合が2720本で、別モデル間の結合が0本
        self.assertEqual(len(py_bonds), 2720)
        py_cross_bonds = [
            b for b in py_bonds if get_model_name(b[0]) != get_model_name(b[1])
        ]
        self.assertEqual(len(py_cross_bonds), 0)

        # 純Python版: ag["model_3"] だけを渡すと 136本
        py_model_3 = py_ag["model_3"]
        self.assertFalse(Format.is_models(py_model_3))
        py_bond_3 = Bond()
        py_bond_3.setup(py_model_3)
        self.assertEqual(len(py_model_3.get_bond_list()), 136)

        # Rust版 (proteindf_bridge.rs) の Pdb と Bond
        rs_pdb = rs_br.Pdb(pdb_path)
        rs_ag = rs_pdb.get_atomgroup()
        rs_bond = rs_br.Bond()
        rs_bond.setup(rs_ag)
        rs_bonds = rs_ag.get_bond_list()

        # Rust版: 結合が2720本で、別モデル間の結合が0本
        self.assertEqual(len(rs_bonds), 2720)
        rs_cross_bonds = [
            b for b in rs_bonds if get_model_name(b[0]) != get_model_name(b[1])
        ]
        self.assertEqual(len(rs_cross_bonds), 0)

        # Rust版: ag["model_3"] だけを渡すと 136本
        rs_model_3 = rs_ag["model_3"]
        rs_bond_3 = rs_br.Bond()
        rs_bond_3.setup(rs_model_3)
        self.assertEqual(len(rs_model_3.get_bond_list()), 136)


if __name__ == "__main__":
    unittest.main()
