#!/usr/bin/env python3
# Times common operations of the installed proteindf_bridge on one structure file.
# Run by devtool/bench-switch.sh inside a venv that has one release installed;
# prints one JSON object (version, file, atom count, seconds per operation).
#
# Usage: python -I bench_switch.py <file.pdb> [--repeat N] [--ops op1,op2,...]
import argparse
import copy
import json
import os
import pickle
import statistics
import sys
import tempfile
import time

import proteindf_bridge as bridge


def walk(ag):
    # Visit every atom and read its coordinates and name, as user scripts do.
    n = 0
    for _, group in ag.groups():
        n += walk(group)
    for _, atom in ag.atoms():
        p = atom.xyz
        _ = (p.x, p.y, p.z, atom.name)
        n += 1
    return n


def make_ops(path, tmpdir):
    ag0 = bridge.Pdb(path).get_atomgroup()
    brd = os.path.join(tmpdir, "bench.brd")
    bridge.save_atomgroup(ag0, brd)
    raw = ag0.get_raw_data()
    blob = pickle.dumps(ag0)

    def write_pdb(ag):
        pdb = bridge.Pdb()
        pdb.set_by_atomgroup(ag)
        return str(pdb)

    # (name, function taking ag0 or a fresh copy, needs a fresh copy)
    ops = [
        ("pdb_load", lambda _: bridge.Pdb(path).get_atomgroup(), False),
        ("pdb_write", lambda ag: write_pdb(ag), False),
        ("brd_save", lambda ag: bridge.save_atomgroup(ag, brd), False),
        ("brd_load", lambda _: bridge.load_atomgroup(brd), False),
        ("dict_from", lambda _: bridge.AtomGroup(raw), False),
        ("dict_to", lambda ag: ag.get_raw_data(), False),
        ("pickle_dumps", lambda ag: pickle.dumps(ag), False),
        ("pickle_loads", lambda _: pickle.loads(blob), False),
        ("deepcopy", lambda ag: copy.deepcopy(ag), False),
        ("walk_atoms", walk, False),
        ("path_list", lambda ag: ag.get_path_list(), False),
        ("count_atoms", lambda ag: ag.get_number_of_all_atoms(), False),
        ("select_name_CA", lambda ag: ag.select(bridge.Select_Name("CA")), False),
        ("pickup_atoms_CA", lambda ag: ag.pickup_atoms("CA"), False),
        ("center", lambda ag: ag.center(), False),
        ("box", lambda ag: ag.box(), False),
        ("shift_by", lambda ag: ag.shift_by(bridge.Position(1.0, 0.0, 0.0)), True),
        ("bond_setup", lambda ag: bridge.Bond().setup(ag), True),
    ]
    # rs.Bond accepts only the Rust AtomGroup, so it is timed only after the switch.
    if bridge.AtomGroup is getattr(getattr(bridge, "rs", None), "AtomGroup", None):
        ops.append(("rs_bond_setup", lambda ag: bridge.rs.Bond().setup(ag), True))
    return ag0, ops


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--repeat", type=int, default=5)
    ap.add_argument("--ops", default="")
    args = ap.parse_args()
    wanted = set(filter(None, args.ops.split(",")))

    with tempfile.TemporaryDirectory() as tmpdir:
        ag0, ops = make_ops(args.file, tmpdir)
        result = {
            "version": bridge.__version__,
            "file": os.path.basename(args.file),
            "atoms": ag0.get_number_of_all_atoms(),
            "repeat": args.repeat,
            "seconds": {},
        }
        for name, func, fresh in ops:
            if wanted and name not in wanted:
                continue
            times = []
            for _ in range(args.repeat):
                ag = copy.deepcopy(ag0) if fresh else ag0
                t0 = time.perf_counter()
                func(ag)
                times.append(time.perf_counter() - t0)
            result["seconds"][name] = {"median": statistics.median(times), "min": min(times)}
            print(f"  {name}: {statistics.median(times):.4f} s", file=sys.stderr)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
