#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite validating rs.AtomGroup mutators, sorting properties, and charge assignment
(TASK_rs-ag-mutators).
"""

import copy
import inspect
import os
import pickle
import unittest
import warnings

try:
    import pytest
except ImportError:
    pytest = None

import proteindf_bridge.rs as rs
from proteindf_bridge.vector import Vector as PyVector


class TestRsAgMutators(unittest.TestCase):
    """Test suite for erase_*, sort_*, and assign_charges methods of rs.AtomGroup."""

    def setUp(self):
        self.data_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "proteindf_bridge", "data")
        )

    def test_criterion1_erase_atom(self):
        """Criterion 1: erase_atom() raises DeprecationWarning and removes the atom."""
        ag = rs.AtomGroup("root")
        a1 = rs.Atom(name="CA", symbol="C")
        a2 = rs.Atom(name="CB", symbol="C")
        ag.set_atom("1", a1)
        ag.set_atom("2", a2)
        self.assertEqual(ag.get_number_of_atoms(), 2)

        # Check pytest.warns(DeprecationWarning) if pytest is available
        if pytest is not None:
            with pytest.warns(
                DeprecationWarning,
                match="erase_atom\\(\\) is deprecated, use remove_atom\\(\\) instead.",
            ):
                ret = ag.erase_atom("1")
        else:
            with warnings.catch_warnings(record=True) as recorded:
                warnings.simplefilter("always")
                ret = ag.erase_atom("1")
                self.assertEqual(len(recorded), 1)
                self.assertEqual(recorded[0].category, DeprecationWarning)
                self.assertIn("erase_atom() is deprecated", str(recorded[0].message))

        self.assertIsNone(ret)
        self.assertEqual(ag.get_number_of_atoms(), 1)
        self.assertFalse(ag.has_atom("1"))
        self.assertTrue(ag.has_atom("2"))

        # Verify warning origin points to the caller line
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always")
            target_line = inspect.currentframe().f_lineno + 1
            ag.erase_atom("2")
            self.assertEqual(len(recorded), 1)
            self.assertEqual(recorded[0].category, DeprecationWarning)
            self.assertIn("erase_atom() is deprecated", str(recorded[0].message))
            self.assertEqual(recorded[0].filename, __file__)
            self.assertEqual(recorded[0].lineno, target_line)

        self.assertEqual(ag.get_number_of_atoms(), 0)

    def test_criterion1_erase_group(self):
        """Criterion 1: erase_group() raises DeprecationWarning and removes the group."""
        ag = rs.AtomGroup("root")
        g1 = rs.AtomGroup("sub1")
        g2 = rs.AtomGroup("sub2")
        ag.set_group("1", g1)
        ag.set_group("2", g2)
        self.assertEqual(ag.get_number_of_groups(), 2)

        # Check pytest.warns(DeprecationWarning) if pytest is available
        if pytest is not None:
            with pytest.warns(
                DeprecationWarning,
                match="erase_group\\(\\) is deprecated, use remove_group\\(\\) instead.",
            ):
                ret = ag.erase_group("1")
        else:
            with warnings.catch_warnings(record=True) as recorded:
                warnings.simplefilter("always")
                ret = ag.erase_group("1")
                self.assertEqual(len(recorded), 1)
                self.assertEqual(recorded[0].category, DeprecationWarning)
                self.assertIn("erase_group() is deprecated", str(recorded[0].message))

        self.assertIsNone(ret)
        self.assertEqual(ag.get_number_of_groups(), 1)
        self.assertFalse(ag.has_group("1"))
        self.assertTrue(ag.has_group("2"))

        # Verify warning origin points to the caller line
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always")
            target_line = inspect.currentframe().f_lineno + 1
            ag.erase_group("2")
            self.assertEqual(len(recorded), 1)
            self.assertEqual(recorded[0].category, DeprecationWarning)
            self.assertIn("erase_group() is deprecated", str(recorded[0].message))
            self.assertEqual(recorded[0].filename, __file__)
            self.assertEqual(recorded[0].lineno, target_line)

        self.assertEqual(ag.get_number_of_groups(), 0)

    def _build_tree(self):
        """Build Rust tree with insertion order 1, 10, 2."""
        rs_root = rs.AtomGroup("root")

        for key in ["1", "10", "2"]:
            rs_grp = rs.AtomGroup(f"grp_{key}")

            for a_key in ["1", "10", "2"]:
                rs_atm = rs.Atom(name=f"A{a_key}", symbol="C")
                rs_grp.set_atom(a_key, rs_atm)

            rs_root.set_group(key, rs_grp)

            # Also top-level atoms
            rs_root.set_atom(f"top_{key}", rs.Atom(name=f"T{key}", symbol="N"))

        return rs_root

    def test_criterion2_sort_properties_default_and_nice(self):
        """Criterion 2: Tree with keys '1', '10', '2' returns 1, 2, 10 with 'nice',

        returns 1, 10, 2 with None.
        """
        rs_root = self._build_tree()

        # Initial state: None
        self.assertIsNone(rs_root.sort_atoms)
        self.assertIsNone(rs_root.sort_groups)

        # 1. Default (None) ordering: insertion order '1', '10', '2'
        self.assertEqual(rs_root.get_group_list(), ["1", "10", "2"])
        self.assertEqual([k for k, _ in rs_root.groups()], ["1", "10", "2"])

        self.assertEqual(rs_root.get_atom_keys(), ["top_1", "top_10", "top_2"])
        self.assertEqual([k for k, _ in rs_root.atoms()], ["top_1", "top_10", "top_2"])

        # get_atom_list / get_path_list
        expected_default_paths = [
            "/1/1", "/1/10", "/1/2",
            "/10/1", "/10/10", "/10/2",
            "/2/1", "/2/10", "/2/2",
            "/top_1", "/top_10", "/top_2",
        ]
        self.assertEqual(
            [a.path for a in rs_root.get_atom_list()],
            expected_default_paths,
        )
        self.assertEqual(rs_root.get_path_list(), expected_default_paths)

        # pickup_atoms
        self.assertEqual([a.path for a in rs_root.pickup_atoms("C")], [])
        self.assertEqual([a.path for a in rs_root.pickup_atoms("A1")], ["/1/1", "/10/1", "/2/1"])

        def get_coords_lines(xyz_text):
            return [l for l in xyz_text.strip().split("\n") if not l.startswith("#") and not l.strip().isdigit()]

        self.assertEqual(len(get_coords_lines(rs_root.get_xyz())), 12)

        # 2. Switch to 'nice'
        rs_root.sort_groups = "nice"
        rs_root.sort_atoms = "nice"

        for k in ["1", "10", "2"]:
            rs_root[k].sort_atoms = "nice"

        self.assertEqual(rs_root.sort_groups, "nice")
        self.assertEqual(rs_root.sort_atoms, "nice")

        # Natural order: '1', '2', '10'
        self.assertEqual(rs_root.get_group_list(), ["1", "2", "10"])
        self.assertEqual([k for k, _ in rs_root.groups()], ["1", "2", "10"])

        self.assertEqual(rs_root.get_atom_keys(), ["top_1", "top_2", "top_10"])
        self.assertEqual([k for k, _ in rs_root.atoms()], ["top_1", "top_2", "top_10"])

        expected_nice_paths = [
            "/1/1", "/1/2", "/1/10",
            "/2/1", "/2/2", "/2/10",
            "/10/1", "/10/2", "/10/10",
            "/top_1", "/top_2", "/top_10",
        ]
        self.assertEqual(
            [a.path for a in rs_root.get_atom_list()],
            expected_nice_paths,
        )
        self.assertEqual(rs_root.get_path_list(), expected_nice_paths)

        self.assertEqual([a.path for a in rs_root.pickup_atoms("C")], [])
        self.assertEqual([a.path for a in rs_root.pickup_atoms("A1")], ["/1/1", "/2/1", "/10/1"])
        self.assertEqual(len(get_coords_lines(rs_root.get_xyz())), 12)

        # 3. Switch back to None -> restore insertion order
        rs_root.sort_groups = None
        rs_root.sort_atoms = None
        for k in ["1", "10", "2"]:
            rs_root[k].sort_atoms = None

        self.assertEqual(rs_root.get_group_list(), ["1", "10", "2"])
        self.assertEqual(rs_root.get_atom_keys(), ["top_1", "top_10", "top_2"])
        self.assertEqual(rs_root.get_path_list(), expected_default_paths)
        self.assertEqual(len(get_coords_lines(rs_root.get_xyz())), 12)

        # 4. Flat tree get_xyz matches expected
        rs_flat = rs.AtomGroup("flat")
        for k in ["1", "10", "2"]:
            rs_flat.set_atom(k, rs.Atom(name=f"A{k}", symbol="C"))

        self.assertEqual(rs_flat.get_atom_keys(), ["1", "10", "2"])
        rs_flat.sort_atoms = "nice"
        self.assertEqual(rs_flat.get_atom_keys(), ["1", "2", "10"])
        rs_flat.sort_atoms = None
        self.assertEqual(rs_flat.get_atom_keys(), ["1", "10", "2"])

    def test_criterion3_copy_and_pickle_preserve_sort_settings(self):
        """Criterion 3: sort_atoms and sort_groups are preserved in copy, deepcopy, and pickle."""
        rs_root = self._build_tree()
        rs_root.sort_groups = "nice"
        rs_root.sort_atoms = "nice"
        rs_root["1"].sort_atoms = "nice"

        # copy.copy
        c_copy = copy.copy(rs_root)
        self.assertEqual(c_copy.sort_groups, "nice")
        self.assertEqual(c_copy.sort_atoms, "nice")
        self.assertEqual(c_copy["1"].sort_atoms, "nice")
        self.assertEqual(c_copy.get_group_list(), ["1", "2", "10"])

        # copy.deepcopy
        c_deep = copy.deepcopy(rs_root)
        self.assertEqual(c_deep.sort_groups, "nice")
        self.assertEqual(c_deep.sort_atoms, "nice")
        self.assertEqual(c_deep["1"].sort_atoms, "nice")
        self.assertEqual(c_deep.get_group_list(), ["1", "2", "10"])

        # pickle
        dumped = pickle.dumps(rs_root)
        c_pickle = pickle.loads(dumped)
        self.assertEqual(c_pickle.sort_groups, "nice")
        self.assertEqual(c_pickle.sort_atoms, "nice")
        self.assertEqual(c_pickle["1"].sort_atoms, "nice")
        self.assertEqual(c_pickle.get_group_list(), ["1", "2", "10"])

    def test_criterion3_assign_charges_2mgo(self):
        """Criterion 3: assign_charges() on 2MGO.pdb assigns expected charges

        with list, rs.Vector, and PyVector, including sort_atoms='nice',
        and raises ValueError on size mismatch.
        """
        pdb_path = os.path.join(self.data_dir, "2MGO.pdb")
        self.assertTrue(os.path.exists(pdb_path), "2MGO.pdb fixture not found")

        rs_ag = rs.Pdb(pdb_path).get_atomgroup()

        n_atoms = rs_ag.get_number_of_all_atoms()
        self.assertEqual(n_atoms, 2680)

        # Generate unique test charges: 0.001 * (i + 1)
        charge_values = [0.001 * (i + 1) for i in range(n_atoms)]
        py_vec = PyVector(charge_values)
        rs_vec = rs.Vector(charge_values)

        # 1. Assign charges using PyVector
        rs_ag.assign_charges(py_vec)
        rs_atoms = rs_ag.get_atom_list()
        self.assertEqual(len(rs_atoms), n_atoms)
        for i, ra in enumerate(rs_atoms):
            self.assertAlmostEqual(ra.charge, charge_values[i], places=9)

        # 2. Assign charges with rs.Vector
        rs_ag2 = rs.Pdb(pdb_path).get_atomgroup()
        rs_ag2.assign_charges(rs_vec)
        for i, ra in enumerate(rs_ag2.get_atom_list()):
            self.assertAlmostEqual(ra.charge, charge_values[i], places=9)

        # 3. Assign charges with Python list
        rs_ag3 = rs.Pdb(pdb_path).get_atomgroup()
        rs_ag3.assign_charges(charge_values)
        for i, ra in enumerate(rs_ag3.get_atom_list()):
            self.assertAlmostEqual(ra.charge, charge_values[i], places=9)

        # 4. With sort_atoms="nice"
        rs_ag_nice = rs.Pdb(pdb_path).get_atomgroup()

        def set_sort_atoms_nice(ag):
            ag.sort_atoms = "nice"
            for _, grp in ag.groups():
                set_sort_atoms_nice(grp)

        set_sort_atoms_nice(rs_ag_nice)
        rs_ag_nice.assign_charges(rs_vec)

        for i, ra in enumerate(rs_ag_nice.get_atom_list()):
            self.assertAlmostEqual(ra.charge, charge_values[i], places=9)

        # 5. Length mismatch raises ValueError
        wrong_charges_short = charge_values[:-1]
        wrong_charges_long = charge_values + [0.5]

        with self.assertRaises(ValueError):
            rs_ag.assign_charges(wrong_charges_short)

        with self.assertRaises(ValueError):
            rs_ag.assign_charges(wrong_charges_long)


if __name__ == "__main__":
    unittest.main()
