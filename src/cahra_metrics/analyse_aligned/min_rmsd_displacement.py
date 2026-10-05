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
import json

import numpy as np
import matplotlib.pyplot as plt
import gemmi

from cahra_metrics.utilities.utils import (
    make_dir,
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
        "--pred_ensemble_json",
        "-pe",
        help=(
            "Path to directory holding mmcif/pdb files which make up the"
            " submitted ensemble. Use one value per predicted ensemble "
            "directory."
        ),
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


def rmsf_to_bfac(rmsfs: list[float]):
    scalar_mult = 8.0 * np.pi * np.pi / 3.0
    bfac = scalar_mult * np.power(np.array(rmsfs, dtype=float), 2)
    bfac = bfac.tolist()
    # round the bfactor values to 2d.p.
    bfac = [round(val, 2) for val in bfac]
    return bfac


def rmsf_round(rmsfs: list[float]):
    rmsfs = [round(val, 2) for val in rmsfs]
    return rmsfs


def set_b_fac_atomlist(st: gemmi.Structure, bfacs: list[float]):
    for model in st:
        for chain in model:
            for res in chain:
                for bfac, atom in zip(bfacs, res):
                    atom.b_iso = bfac
    return st


def save_b_fac_pdb(
    dirname: str,
    pdbsavename: str,
    st: gemmi.Structure,
    subdir: str,
    bfacs: list[float],
):
    # check whether bfacs are passed
    # if not, don't update them
    if bfacs:
        bfac_st = set_b_fac_atomlist(
            st,
            bfacs,
        )
    else:
        bfac_st = st

    bfac_dir = os.path.join(
        os.path.dirname(dirname),
        subdir,
    )
    make_dir(bfac_dir)
    bfac_file = os.path.join(
        bfac_dir,
        pdbsavename,
    )

    bfac_st.write_pdb(bfac_file)


def plot_violin_noextras(
    labels: list[str],
    val_lists: list[list[float]],
    xlabel: str,
    ylabel: str,
    plot_dir: str,
    plot_name: str,
):
    # check that vals are not all zero
    mask = []
    for vals in val_lists:
        check_zeros = np.zeros(len(vals)).astype(float)
        if np.allclose(vals, check_zeros):
            mask.append(False)
        else:
            mask.append(True)

    # mask out labels
    labels = np.array(labels, dtype=str)[np.asarray(mask, dtype=bool)].tolist()
    # mask out val lists
    val_lists_masked = []
    for vals, mask_val in zip(val_lists, mask):
        if mask_val:
            val_lists_masked.append(vals)

    # create a fresh figure sized for the number of datasets
    plt.figure(figsize=(max(8, 0.8 * len(labels)), 6))
    # violin plot of the min-RMSD distribution for each dataset
    parts = plt.violinplot(
        val_lists_masked,
        positions=np.arange(len(labels)),
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body in parts["bodies"]:
        body.set_facecolor("steelblue")
        body.set_alpha(0.6)
        body.set_edgecolor("black")
    # overlaid box plot (median, mean, quartiles, whiskers, outliers) per dataset
    plt.boxplot(
        val_lists_masked,
        positions=np.arange(len(labels)),
        widths=0.15,
        patch_artist=True,
        showmeans=True,
        meanline=True,
        meanprops=dict(color="blue", linewidth=2, label="Mean"),
        boxprops=dict(facecolor="white", edgecolor="black"),
        medianprops=dict(color="red", linewidth=2, label="Median"),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
        flierprops=dict(marker="x", markerfacecolor="black", markersize=2),
    )
    plt.xticks(np.arange(len(labels)), labels, rotation=45, ha="right")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # grid with both major and minor ticks
    plt.grid(True, which="major", linestyle="--", alpha=0.5, axis="y")
    plt.grid(True, which="minor", linestyle=":", alpha=0.3, axis="y")
    # add legend for mean and median
    plt.legend(
        handles=[
            plt.Line2D([0], [0], color="blue", linewidth=2, label="Mean"),
            plt.Line2D([0], [0], color="red", linewidth=2, label="Median"),
        ],
        loc="upper right",
    )
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, plot_name + ".png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, plot_name + ".pdf"))
    plt.clf()


# compare b-factors between ensembles
"""
# compare backbone and sidechain b-factors for a given ensemble
def make_bfac_comp_plot(
    bfac_labels: list[str],
    bfacs: list[list[float]],
    xlabel: str,
    ylabel: str,
    plot_dir: str,
    plot_name: str,
    avg_over: int = 10,
    bb_mask: list[bool] = [],
):
    smoothed_bfacs = []
    for
    avgs = []
    for i, atom in enumerate(range(len(diff))):
        if i < avg_over:
            continue
        else:
            avgs.append(np.average(diff[i - avg_over : i]))

    avgs_arr = np.array(avgs)
    plt.plot(range(len(diff)), diff, color="red", label="GT-pred")
    plt.plot(
        range(len(diff))[
            round(avg_over / 2.0) : len(diff) - round(avg_over / 2.0)
        ],
        avgs_arr,
        color="black",
        label="100-atom-avg",
    )
    plt.xlabel("Atom number")
    plt.ylabel("RMSF ($\\AA$)")
    plt.axvspan(79, 2175, zorder=0, alpha=0.2, color="orange", label="D1")
    plt.axvspan(4346, 4366, zorder=0, alpha=0.2, color="blue", label="LNK")
    plt.axvspan(6465, 8085, zorder=0, alpha=0.2, color="green", label="D4")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "rmsf_diff.png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, "rmsf_diff.pdf"))
    plt.clf()
"""


def get_name():
    return "min_rmsd_displacement"


def main(args):
    # STEP 1 - load the json for each pred ensemble
    if args.verbose:
        print(args.pred_ensemble_json)

    # Create lists for keeping track of min-RMSD atomic displacements
    # One list per pred_ensemble. The sublist will be appended to be
    # len (n_atoms*mask) where mask is applicable for bb or sc
    min_rmsd_disps: list[list[float]] = []
    min_rmsd_labels: list[str] = []
    backbone_names = {"N", "CA", "C", "O", "OXT"}
    min_rmsd_bb_disps: list[list[float]] = []  # for backbone atoms
    min_rmsd_sc_disps: list[list[float]] = []  # for sidechain atoms

    # Create lists for keeping track of b-factor values for each
    # pred ensemble and for GT (if applicable)
    bfacs: list[list[float]] = []  # one list per pred ens plus gt ens
    bfacs_labels: list[str] = []  # gt and pred labels
    bfacs_diff: list[list[float]] = []  # one list per pred ens
    bfacs_abs_diff: list[list[float]] = []  # one list per pred ens
    bfacs_diff_labels: list[str] = []  # pred labels
    bb_mask: list[bool] = (
        []
    )  # fill mask to allow filtering of bfacs by bb or sc
    residue_mask: list[int] = []  # fill with residue index for each atom

    # Loop over the pred ensembles
    for i, (json_file, label) in enumerate(
        zip(
            args.pred_ensemble_json,
            args.pred_ensemble_label,
        )
    ):
        # STEP 2 - Open the json file
        if args.verbose:
            print("Loading json: {}".format(json_file))
            print("Label: {}".format(label))
        with open(json_file, "r") as f:
            results = json.load(f)

        # collect all min_rmsd_displacements for all models in ensemble
        # in a single list to make violinplot
        min_rmsd_disp: list[float] = []
        min_rmsd_bb_disp: list[float] = []
        min_rmsd_sc_disp: list[float] = []
        min_rmsd_labels.append(label)

        # STEP 3 - Get PDB version with per-atom-distances as b-factor col
        # Loop through pred_filenames
        for idx, pred_pdb in enumerate(results["pred_filenames"]):
            # Load the pred_pdb with Gemmi
            pred_st = gemmi.read_structure(pred_pdb)

            # Get min-RMSD GT model index
            min_rmsd_idx = results["ensemble"]["conf_stats"]["gt_index"][idx]

            # Load the gt_pdb with Gemmi
            gt_st = gemmi.read_structure(results["gt_filenames"][min_rmsd_idx])

            # and the min-RMSD GT model
            # There will be 1 model in each pdb
            atom_distances = []
            residue_rmsfs = []
            for pred_model, gt_model in zip(pred_st, gt_st):
                for pred_chain, gt_chain in zip(pred_model, gt_model):
                    for pred_res, gt_res in zip(pred_chain, gt_chain):
                        residue_rmsf = []
                        for pred_atom, gt_atom in zip(pred_res, gt_res):
                            pred_pos = np.array(
                                [
                                    pred_atom.pos.x,
                                    pred_atom.pos.y,
                                    pred_atom.pos.z,
                                ]
                            )
                            gt_pos = np.array(
                                [
                                    gt_atom.pos.x,
                                    gt_atom.pos.y,
                                    gt_atom.pos.z,
                                ]
                            )
                            # not that we only keep 2d.p. for distance
                            atom_distance = round(
                                float(
                                    np.absolute(
                                        np.linalg.norm(pred_pos - gt_pos)
                                    )
                                ),
                                2,
                            )
                            atom_distances.append(atom_distance)
                            min_rmsd_disp.append(atom_distance)
                            #
                            if (
                                pred_atom.name not in backbone_names
                                and not pred_atom.is_hydrogen()
                            ):
                                min_rmsd_sc_disp.append(atom_distance)
                                # get mask for bb/sc if first label+model
                                if i == 0 and idx == 0:
                                    bb_mask.append(False)
                            else:
                                min_rmsd_bb_disp.append(atom_distance)
                                # get mask for bb/sc if first label+model
                                if i == 0 and idx == 0:
                                    bb_mask.append(True)
                            # also get residue index mask if first label+model
                            if i == 0 and idx == 0:
                                residue_mask.append(pred_res.seqid.num)
                            # Add the distance for each atom (A)
                            # to the b-factor column for pred model
                            pred_atom.b_iso = atom_distance

                        residue_rmsfs.append(np.average(residue_rmsf))

            # Make per-residue LFC plot

            # TODO Add per model min-RMSD model atom distances to json?

            # TODO Save updated json?

            # Save pred_pdb atomic model with b-factor col filled
            # formulate save directory and fname from pred_pdb filename
            save_b_fac_pdb(
                os.path.dirname(pred_pdb),
                os.path.basename(pred_pdb),
                pred_st,
                "atom_dist",
                [],
            )

            # Look at b-factors only once per pred ensemble
            # STEP 4 - Get PDB version with
            # Load RMSF for GT models and convert to b-factor
            gt_bfacs = rmsf_to_bfac(
                results["ensemble"]["ensemble_stats"]["gt_rmsf"]
            )
            gt_rmsfs = rmsf_round(
                results["ensemble"]["ensemble_stats"]["gt_rmsf"]
            )

            # Load RMSF for pred models and convert to b-factor
            pred_bfacs = rmsf_to_bfac(
                results["ensemble"]["ensemble_stats"]["pred_rmsf"]
            )
            pred_rmsfs = rmsf_round(
                results["ensemble"]["ensemble_stats"]["pred_rmsf"]
            )

            # Compute difference in b-factor (pred-GT)
            diff_bfacs = np.array(pred_bfacs, dtype=float) - np.array(
                gt_bfacs, dtype=float
            )
            diff_rmsfs = np.array(pred_rmsfs, dtype=float) - np.array(
                gt_rmsfs, dtype=float
            )

            # Compute abs difference in b-factor (pred-GT)
            abs_diff_bfacs = np.absolute(diff_bfacs).tolist()
            abs_diff_rmsfs = np.absolute(diff_rmsfs).tolist()

            # Negative vals can be used for plots but not for pdbs
            diff_bfacs = diff_bfacs.tolist()
            diff_rmsfs = diff_rmsfs.tolist()

            # If first iteration, save an updated PDB file for GT
            # b-factor additions
            if i == 0 and idx == 0:
                bfacs.append(gt_bfacs)
                bfacs_labels.append("GT")
                save_b_fac_pdb(
                    os.path.dirname(results["gt_filenames"][min_rmsd_idx]),
                    os.path.basename(results["gt_filenames"][min_rmsd_idx]),
                    gt_st,
                    "b_facs",
                    gt_bfacs,
                )
                # save gt rmsfs values into pdb
                save_b_fac_pdb(
                    os.path.dirname(results["gt_filenames"][min_rmsd_idx]),
                    os.path.basename(results["gt_filenames"][min_rmsd_idx]),
                    gt_st,
                    "rmsfs",
                    gt_rmsfs,
                )

            # If first iteration of pred_ensemble
            # Save the GT model with b-factor abs-diff from pred model
            if idx == 0:
                bfacs.append(pred_bfacs)
                bfacs_labels.append(label)
                bfacs_diff_labels.append(label)
                bfacs_diff.append(diff_bfacs)
                bfacs_abs_diff.append(abs_diff_bfacs)
                # abs b-factor differences
                # This depends on pred_ens, so put in pred subdir
                save_b_fac_pdb(
                    os.path.dirname(pred_pdb),
                    os.path.basename(results["gt_filenames"][min_rmsd_idx]),
                    gt_st,
                    "b_facs_absdiff",
                    abs_diff_bfacs,
                )
                # save abs diff rmsfs in gt pdb
                save_b_fac_pdb(
                    os.path.dirname(pred_pdb),
                    os.path.basename(results["gt_filenames"][min_rmsd_idx]),
                    gt_st,
                    "rmsfs_absdiff",
                    abs_diff_rmsfs,
                )

                # Can't save negative b-factors to pdb file?
                # so can't also save b-factor differences to pdbs

                # If first iteration also save the pred model
                # with pred ensembles b-factors
                # This will be closest pred model to gt model
                # but will not always be super close for qual
                # eval, so don't rely on that
                save_b_fac_pdb(
                    os.path.dirname(pred_pdb),
                    os.path.basename(pred_pdb),
                    pred_st,
                    "b_facs",
                    pred_bfacs,
                )
                # save rmsfs of pred ensemble
                save_b_fac_pdb(
                    os.path.dirname(pred_pdb),
                    os.path.basename(pred_pdb),
                    pred_st,
                    "rmsfs",
                    pred_rmsfs,
                )

        # Extend list of atomic displacements for this pred_ensemble
        min_rmsd_disps.append(min_rmsd_disp)
        min_rmsd_bb_disps.append(min_rmsd_bb_disp)
        min_rmsd_sc_disps.append(min_rmsd_sc_disp)

    # Plot violinplots of per-atom min_rmsd_displacements
    # One violin per pred ensemble
    # for all atoms
    plot_violin_noextras(
        shorten_labels(min_rmsd_labels),
        min_rmsd_disps,
        "Ensemble",
        "Atom distance to min-RMSD GT model (Å)",
        results["pred_confs_dir"],
        "min_rmsd_atom_dists",
    )
    # for backbone
    plot_violin_noextras(
        shorten_labels(min_rmsd_labels),
        min_rmsd_bb_disps,
        "Ensemble",
        "Backbone atom distance to min-RMSD GT model (Å)",
        results["pred_confs_dir"],
        "min_rmsd_bb_atom_dists",
    )
    # for sidechain
    plot_violin_noextras(
        shorten_labels(min_rmsd_labels),
        min_rmsd_sc_disps,
        "Ensemble",
        "Sidechain atom distance to min-RMSD GT model (Å)",
        results["pred_confs_dir"],
        "min_rmsd_sc_atom_dists",
    )
    # The backbone and sidechain plots let you know relative spread

    # Plot violinplots of bfac distributions for each ensemble
    # (including GT). This shows range of flexibility globally.
    # for all atoms
    plot_violin_noextras(
        shorten_labels(bfacs_labels),
        bfacs,
        "Ensemble",
        "B-factor (atom) (Å^2)",
        results["pred_confs_dir"],
        "bfacs_all_atoms",
    )
    # for backbone
    bfacs_bb_only = []
    for bfac in bfacs:
        bfacs_bb_only.append([bf for mask, bf in zip(bb_mask, bfac) if mask])
    plot_violin_noextras(
        shorten_labels(bfacs_labels),
        bfacs_bb_only,
        "Ensemble",
        "B-factor backbone (atom) (Å^2)",
        results["pred_confs_dir"],
        "bfacs_bb_atoms",
    )
    # for sidechain
    bfacs_sc_only = []
    for bfac in bfacs:
        bfacs_sc_only.append(
            [bf for mask, bf in zip(bb_mask, bfac) if not mask]
        )
    plot_violin_noextras(
        shorten_labels(bfacs_labels),
        bfacs_sc_only,
        "Ensemble",
        "B-factor sidechain (atom) (Å^2)",
        results["pred_confs_dir"],
        "bfacs_sc_atoms",
    )
    # The backbone and sidechain plots let you know relative spread

    # Think about how you want to plot the b-factors clearly.
    # Do you want to use CA atoms only?
    # Which regions of the 6xou do you want to highlight?
    # D1, D4 lnk, but also switches?

    # TODO Load the saved pdb files with RMSF values
    # GT

    # pred ensembles

    # TODO Plot smoothed scatterplots (lines) showing b-factors
    # along the atoms in the atomic model. This shows local
    # flexibility
    # TODO for all atoms
    # TODO for backbone
    # TODO for sidechain

    # TODO Plot smoothed scatterplots (lines) showing b-factor
    # diff (pred-gt) along the atoms in the atomic model.
    # This shows where the pred ensemble is more flexible
    # than GT ensemble (+ve) and where it is less flex (-ve)
    # TODO for all atoms
    # TODO for backbone
    # TODO for sidechain
