"""CAHRA_metrics:
A set of algorithms and utilities for calculation of metrics
for the 2025 CAHRA-Challenge.
Copyright (C) 2026  Joel Greer(UKRI), Tom Burnley (UKRI),
David Silva Sanchez (Flatiron Institute), Pilar Cossio (Flatiron Institute)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

https://opensource.org/license/mit

"""

import argparse
import os

# import json

# import numpy as np
import matplotlib.pyplot as plt
import gemmi

from cahra_metrics.utilities.utils import (
    # make_dir,
    shorten_labels,
)


def add_arguments(parser: argparse.ArgumentParser):
    """Set up arguments for parsing

    Args:
        parser (argparse.ArgumentParser): Waymarking argument parser

    Returns:
        argparse.ArgumentParser:  Waymarking argument parser with arguments
    """

    parser.add_argument(
        "--pred_pdb_rmsf",
        "-pe",
        help=("Path to directory holding pdb files with RMSF values"),
        type=str,
        nargs="+",  # load multiple ensembles to compare against each other
    )

    parser.add_argument(
        "--pred_ensemble_label",
        "-pl",
        help=(
            "Label for predicted ensemble. Used to set output directories."
            " Use one value per predicted ensemble directory."
        ),
        type=str,
        nargs="+",  # one for each pred ensemble
    )

    parser.add_argument(
        "--verbose",
        "-v",
        help="Increase output verbosity",
        action="store_true",
    )


def make_plot_rmsf(
    ad_rmsf_list: list[list[float]],
    labels: list[str],
    plot_dir: str,
    plot_name: str = "rmsf_residue",
):
    """Plot multiple RMSF curves on the same axes.

    Args:
        ad_rmsf_list: List of RMSF value lists, one per ensemble.
        labels: Label for each RMSF curve (must match length of ad_rmsf_list).
        plot_dir: Directory to save the plot in.
        plot_name: Base name for the output files.
    """
    # plot each RMSF curve
    for ad_rmsf, label in zip(ad_rmsf_list, labels):
        plt.plot(
            range(len(ad_rmsf)),
            ad_rmsf,
            label=label,
        )

    plt.xlabel("Atom index")
    plt.ylabel("RMSF ($\\AA$)")
    plt.axvspan(10, 265, zorder=0, alpha=0.2, color="orange", label="D1")
    plt.axvspan(539, 541, zorder=0, alpha=0.2, color="blue", label="LNK")
    plt.axvspan(832, 1014, zorder=0, alpha=0.2, color="green", label="D4")

    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, plot_name + ".png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, plot_name + ".pdf"))
    plt.clf()


def make_plot_rmsf_switch(
    ad_rmsf_list: list[list[float]],
    labels: list[str],
    plot_dir: str,
    plot_name: str = "rmsf_residue_switch",
):
    """Plot multiple RMSF curves on the same axes.

    Args:
        ad_rmsf_list: List of RMSF value lists, one per ensemble.
        labels: Label for each RMSF curve (must match length of ad_rmsf_list).
        plot_dir: Directory to save the plot in.
        plot_name: Base name for the output files.
    """
    # plot each RMSF curve
    for ad_rmsf, label in zip(ad_rmsf_list, labels):
        plt.plot(
            range(len(ad_rmsf)),
            ad_rmsf,
            label=label,
        )

    plt.xlabel("Atom index")
    plt.ylabel("RMSF ($\\AA$)")
    plt.axvspan(10, 265, zorder=0, alpha=0.2, color="orange", label="D1")
    plt.axvspan(539, 541, zorder=0, alpha=0.2, color="blue", label="LNK")
    plt.axvspan(832, 1014, zorder=0, alpha=0.2, color="green", label="D4")

    plt.axvspan(151, 202, zorder=0, alpha=0.2, color="red", label="Switch A")
    plt.axvspan(
        398, 413, zorder=0, alpha=0.2, color="yellow", label="Switch B"
    )
    plt.axvspan(
        483, 527, zorder=0, alpha=0.2, color="purple", label="Switch C"
    )

    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, plot_name + ".png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, plot_name + ".pdf"))
    plt.clf()


def get_per_res_rmsf(pdbfile: str) -> list[float]:
    st = gemmi.read_structure(pdbfile)
    model = st[0]
    # we have single model, single chain
    res_bfacs = []
    # Iterate through chains, residues, and atoms
    for chain in model:
        for res in chain:
            # Collect filled B-factors for all atoms in the residue
            b_factors = [atom.b_iso for atom in res]
            avg_b = sum(b_factors) / len(b_factors)
            res_bfacs.append(avg_b)
    return res_bfacs


def get_name():
    return "plot_rmsf"


def main(args):
    # STEP 1 - load the json for each pred ensemble
    if args.verbose:
        print(args.pred_pdb_rmsf)
        print(args.pred_ensemble_label)

    # list full of per residue rmsf values for each ensemble
    res_rmsfs: list[list[float]] = []
    # loop through the pdb files and get per residue rmsf
    for idx, (pdb, label) in enumerate(
        zip(args.pred_pdb_rmsf, args.pred_ensemble_label)
    ):
        res_rmsfs.append(get_per_res_rmsf(pdb))
    make_plot_rmsf(
        res_rmsfs,
        shorten_labels(args.pred_ensemble_label),
        os.path.dirname(
            os.path.dirname(os.path.dirname(args.pred_ensemble_label[0]))
        ),
    )
    make_plot_rmsf_switch(
        res_rmsfs,
        shorten_labels(args.pred_ensemble_label),
        os.path.dirname(
            os.path.dirname(os.path.dirname(args.pred_ensemble_label[0]))
        ),
    )
