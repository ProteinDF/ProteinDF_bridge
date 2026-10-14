#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Calculate Ramachandran angles (phi/psi) from bridge file."""

import argparse
import csv
import sys

import proteindf_bridge.rs as rs


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="calculate Ramachandran angles (phi/psi) from bridge file"
    )
    parser.add_argument("FILE", help="bridge file")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="output CSV file (default: stdout)",
    )
    parser.add_argument(
        "--plot",
        type=str,
        default=None,
        metavar="PATH",
        help="output plot file path (e.g. .png, .svg, .pdf)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="verbose output",
    )
    return parser.parse_args(argv)


def _format_angle(val):
    if val is None:
        return ""
    return f"{val:.3f}"


def calculate_angles(models):
    """Calculate phi/psi angles for each model and chain.

    Returns a list of dicts with keys:
    model, chain, residue_key, residue_name, phi, psi
    """
    records = []
    for model_key, model in models.groups():
        for chain_key, chain in model.groups():
            angles = rs.calc_phi_psi(chain)
            for angle in angles:
                records.append(
                    {
                        "model": model_key,
                        "chain": chain_key,
                        "residue_key": angle.residue_key,
                        "residue_name": angle.residue_name,
                        "phi": angle.phi,
                        "psi": angle.psi,
                    }
                )
    return records


def write_csv(records, output_path=None):
    """Write phi/psi records to CSV file or stdout."""
    headers = ["model", "chain", "residue_key", "residue_name", "phi", "psi"]
    out_f = None
    try:
        if output_path:
            out_f = open(output_path, "w", newline="", encoding="utf-8")
            writer = csv.writer(out_f)
        else:
            writer = csv.writer(sys.stdout)

        writer.writerow(headers)
        for rec in records:
            writer.writerow(
                [
                    rec["model"],
                    rec["chain"],
                    rec["residue_key"],
                    rec["residue_name"],
                    _format_angle(rec["phi"]),
                    _format_angle(rec["psi"]),
                ]
            )
    finally:
        if out_f is not None:
            out_f.close()


def _ensure_matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        sys.stderr.write(
            "matplotlib is required for plotting. Please install it with: pip install proteindf_bridge[plot]\n"
        )
        sys.exit(1)


def generate_plot(records, plot_path, verbose=False):
    """Generate Ramachandran scatter plot and save to plot_path."""
    plt = _ensure_matplotlib()

    if verbose:
        print(f"generating plot: {plot_path}", file=sys.stderr)

    general_phi, general_psi = [], []
    gly_phi, gly_psi = [], []
    pro_phi, pro_psi = [], []

    for rec in records:
        phi = rec["phi"]
        psi = rec["psi"]
        if phi is None or psi is None:
            continue
        res_name = (rec["residue_name"] or "").strip().upper()
        if res_name == "GLY":
            gly_phi.append(phi)
            gly_psi.append(psi)
        elif res_name == "PRO":
            pro_phi.append(phi)
            pro_psi.append(psi)
        else:
            general_phi.append(phi)
            general_psi.append(psi)

    fig, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(general_phi, general_psi, c="tab:blue", marker="o", s=20, alpha=0.7, label="General")
    ax.scatter(gly_phi, gly_psi, c="tab:green", marker="^", s=25, alpha=0.8, label="GLY")
    ax.scatter(pro_phi, pro_psi, c="tab:orange", marker="s", s=25, alpha=0.8, label="PRO")

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ticks = [-180, -120, -60, 0, 60, 120, 180]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xlabel(r"$\phi$ (degrees)")
    ax.set_ylabel(r"$\psi$ (degrees)")
    ax.set_aspect("equal")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(plot_path)
    plt.close(fig)


def main(argv=None):
    args = parse_args(argv)
    brd_path = args.FILE
    output_path = args.output
    plot_path = args.plot
    verbose = args.verbose

    if plot_path:
        _ensure_matplotlib()

    if verbose:
        print(f"reading: {brd_path}", file=sys.stderr)

    models = rs.load_atomgroup(brd_path)
    records = calculate_angles(models)

    if verbose and output_path:
        print(f"writing CSV: {output_path}", file=sys.stderr)

    write_csv(records, output_path)

    if plot_path:
        generate_plot(records, plot_path, verbose=verbose)


if __name__ == "__main__":
    main()
