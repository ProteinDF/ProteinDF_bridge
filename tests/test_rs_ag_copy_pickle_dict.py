#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite validating rs.AtomGroup copy, pickle, and dictionary roundtrip (TASK_rs-ag-copy-pickle-dict).
"""

import copy
import os
import pickle
import tempfile
import unittest

from proteindf_bridge import Pdb as PyPdb, load_msgpack, save_msgpack
import proteindf_bridge.rs as rs


class TestRsAgCopyPickleDict(unittest.TestCase):
    """Test suite for rs.AtomGroup copy, deepcopy, pickle, and dict conversions."""

    @classmethod
    def setUpClass(cls):
        # Locate 1hls.pdb fixture
        candidate_paths = [
            os.path.join(
                os.path.dirname(__file__), "..", "rust", "crates", "proteindf-bridge", "tests", "data", "1hls.pdb"
            ),
            os.path.join(
                os.path.dirname(__file__), "..", "proteindf_bridge", "data", "1hls.pdb"
            ),
        ]
        cls.pdb_path = None
        for p in candidate_paths:
            if os.path.exists(p):
                cls.pdb_path = os.path.abspath(p)
                break
        assert cls.pdb_path is not None, "1hls.pdb fixture not found"

        # Generate 1hls raw dict via pure-Python Pdb and load_msgpack
        py_pdb = PyPdb()
        py_pdb.load(cls.pdb_path)
        cls.py_1hls_ag = py_pdb.get_atomgroup()

        with tempfile.NamedTemporaryFile(suffix=".brd", delete=False) as f:
            cls.brd_temp_path = f.name
        try:
            save_msgpack(cls.py_1hls_ag.get_raw_data(), cls.brd_temp_path)
            cls.raw_1hls = load_msgpack(cls.brd_temp_path)
        finally:
            if os.path.exists(cls.brd_temp_path):
                os.remove(cls.brd_temp_path)

    def test_criterion1_1hls_dict_roundtrip(self):
        """Criterion 1: rs.AtomGroup(raw).get_raw_data() == raw."""
        raw = self.raw_1hls

        # Construct rs.AtomGroup from raw dict
        rs_ag = rs.AtomGroup(raw)
        rs_raw = rs_ag.get_raw_data()

        # rs.AtomGroup(raw).get_raw_data() == raw
        self.assertEqual(rs_raw, raw)

    def test_criterion2_copy_deepcopy_independence(self):
        """Criterion 2: copy.copy and copy.deepcopy preserve get_raw_data() and are deeply independent."""
        raw = self.raw_1hls
        rs_ag = rs.AtomGroup(raw)

        # Shallow copy
        cp = copy.copy(rs_ag)
        self.assertIsNone(cp.parent)
        self.assertEqual(cp.path, rs_ag.path)
        self.assertEqual(cp.get_raw_data(), raw)

        # Deep copy
        dcp = copy.deepcopy(rs_ag)
        self.assertIsNone(dcp.parent)
        self.assertEqual(dcp.path, rs_ag.path)
        self.assertEqual(dcp.get_raw_data(), raw)

        # Mutate copy atom coordinates -> original remains unchanged
        orig_x = rs_ag["model_1"]["A"]["1"]["1_N"].xyz.x
        cp["model_1"]["A"]["1"]["1_N"].xyz.x = orig_x + 100.0
        self.assertEqual(rs_ag["model_1"]["A"]["1"]["1_N"].xyz.x, orig_x)
        self.assertEqual(cp["model_1"]["A"]["1"]["1_N"].xyz.x, orig_x + 100.0)

        # Mutate deepcopy atom coordinates -> original remains unchanged
        dcp["model_1"]["A"]["1"]["1_N"].xyz.x = orig_x + 200.0
        self.assertEqual(rs_ag["model_1"]["A"]["1"]["1_N"].xyz.x, orig_x)
        self.assertEqual(dcp["model_1"]["A"]["1"]["1_N"].xyz.x, orig_x + 200.0)

        # Mutate child groups in copy -> original remains unchanged
        cp.del_group("model_1")
        self.assertTrue(rs_ag.has_group("model_1"))
        self.assertFalse(cp.has_group("model_1"))

        dcp.del_group("model_1")
        self.assertTrue(rs_ag.has_group("model_1"))
        self.assertFalse(dcp.has_group("model_1"))

    def test_criterion2_pickle_roundtrip_and_independence(self):
        """Criterion 2: pickle roundtrip preserves get_raw_data(), has no parent, and is independent."""
        raw = self.raw_1hls
        rs_ag = rs.AtomGroup(raw)

        dumped = pickle.dumps(rs_ag)
        restored = pickle.loads(dumped)

        self.assertIsNone(restored.parent)
        self.assertEqual(restored.path, rs_ag.path)
        self.assertEqual(restored.name, rs_ag.name)
        self.assertEqual(restored.get_raw_data(), raw)

        # Independence: mutate restored -> original is unchanged
        orig_x = rs_ag["model_1"]["A"]["1"]["1_N"].xyz.x
        restored["model_1"]["A"]["1"]["1_N"].xyz.x = orig_x + 300.0
        self.assertEqual(rs_ag["model_1"]["A"]["1"]["1_N"].xyz.x, orig_x)
        self.assertEqual(restored["model_1"]["A"]["1"]["1_N"].xyz.x, orig_x + 300.0)

    def test_criterion2_subtree_copy_and_pickle_path_preservation(self):
        """Criterion 2: Subtree ag['model_1']['A'] copy/pickle retains path='/model_1/A/' and has parent=None."""
        raw = self.raw_1hls
        rs_ag = rs.AtomGroup(raw)

        chain_a = rs_ag["model_1"]["A"]
        self.assertEqual(chain_a.path, "/model_1/A/")
        self.assertIsNotNone(chain_a.parent)

        # 1. copy.copy on subtree
        cp = copy.copy(chain_a)
        self.assertIsNone(cp.parent)
        self.assertEqual(cp.path, "/model_1/A/")
        self.assertEqual(cp["1"].path, "/model_1/A/1/")
        self.assertEqual(cp.get_raw_data(), chain_a.get_raw_data())

        # 2. copy.deepcopy on subtree
        dcp = copy.deepcopy(chain_a)
        self.assertIsNone(dcp.parent)
        self.assertEqual(dcp.path, "/model_1/A/")
        self.assertEqual(dcp["1"].path, "/model_1/A/1/")
        self.assertEqual(dcp.get_raw_data(), chain_a.get_raw_data())

        # 3. pickle on subtree
        pickled = pickle.loads(pickle.dumps(chain_a))
        self.assertIsNone(pickled.parent)
        self.assertEqual(pickled.path, "/model_1/A/")
        self.assertEqual(pickled["1"].path, "/model_1/A/1/")
        self.assertEqual(pickled.get_raw_data(), chain_a.get_raw_data())

        # Subtree mutations do not affect original
        orig_x = chain_a["1"]["1_N"].xyz.x
        cp["1"]["1_N"].xyz.x = orig_x + 10.0
        dcp["1"]["1_N"].xyz.x = orig_x + 20.0
        pickled["1"]["1_N"].xyz.x = orig_x + 30.0
        self.assertEqual(chain_a["1"]["1_N"].xyz.x, orig_x)

    def test_criterion3_pickle_preserves_secondary_structure(self):
        """Criterion 3: pickle roundtrip preserves secondary structure across the hierarchy."""
        pdb = rs.Pdb()
        pdb.load(self.pdb_path)
        ag = pdb.get_atomgroup()
        chain_a = ag["model_1"]["A"]
        chain_a.apply_secondary_structure()

        self.assertEqual(chain_a["1"].secondary_structure, "-")
        self.assertEqual(chain_a["3"].secondary_structure, "H")

        # Pickle full tree
        restored = pickle.loads(pickle.dumps(ag))
        restored_res1 = restored["model_1"]["A"]["1"]
        restored_res3 = restored["model_1"]["A"]["3"]
        self.assertEqual(restored_res1.secondary_structure, "-")
        self.assertEqual(restored_res3.secondary_structure, "H")

        # Pickle subtree with secondary structure
        restored_chain = pickle.loads(pickle.dumps(chain_a))
        self.assertEqual(restored_chain.path, "/model_1/A/")
        self.assertEqual(restored_chain["1"].secondary_structure, "-")
        self.assertEqual(restored_chain["3"].secondary_structure, "H")

    def test_set_by_dict_data(self):
        """Verify set_by_dict_data populates an empty AtomGroup and supports chaining."""
        raw = self.raw_1hls
        ag = rs.AtomGroup("empty")
        ret = ag.set_by_dict_data(raw)
        self.assertIs(ret, ag)
        self.assertEqual(ag.get_raw_data(), raw)
        self.assertTrue(ag.has_group("model_1"))


if __name__ == "__main__":
    unittest.main()
