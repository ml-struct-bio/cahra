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

# TODO fill docstrings

import argparse
import os
import copy
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from cahra_metrics.utilities.utils import (
    get_ensemble,
    replace_file_ext,
    Ensemble_Stats,
    save_dict_to_json,
    extract_min_rmsd_values,
)


def add_arguments(parser: argparse.ArgumentParser):
    """Set up arguments for parsing

    Args:
        parser (argparse.ArgumentParser): Waymarking argument parser

    Returns:
        argparse.ArgumentParser:  Waymarking argument parser with arguments
    """

    parser.add_argument(
        "--selection",
        "-s",
        help=(
            "Atom selection to make for alignment AND min-RMSD calculation."
            " Models aligned using this min-RMSD alignment selection will be"
            " used as input for PCV calculation. Defaults to 'protein'"
        ),
        type=str,
        choices=["protein", "backbone", "name CA"],
        default="protein",
    )

    parser.add_argument(
        "--exclude_hist_tag",
        "-eht",
        help="Do not include hist-tag in alignment of atomic models and SAC/LFC"
        " calculation",
        action="store_true",
    )

    parser.add_argument(
        "--pred_ensemble_dir",
        "-pe",
        help=(
            "Path to directory holding mmcif/pdb files which make up the"
            " submitted ensemble."
        ),
        type=str,
    )

    parser.add_argument(
        "--pred_ensemble_label",
        "-pl",
        help=(
            "Label for predicted ensemble. Used to create output directory."
        ),
        type=str,
    )

    parser.add_argument(
        "--truth_ensemble_dir",
        "-te",
        help=(
            "Path to directory holding mmcif/pdb files which make up the"
            " truth ensemble."
        ),
        type=str,
        default=(
            "/media/jg/My Book/JGreer/dinaMISMO/dinaMISMO-systems-prep/PreP/"
            "pilot_studies/ConformationSampling/extract_frames_6xou_unsteered_0-200ns/"
            "extract_frames_6xou_unsteered_0-200ns"
        ),
    )

    parser.add_argument(
        "--output_dir",
        "-o",
        help=("Directory to output results in."),
        type=str,
        default="calc_min_rmsd_metrics_outputs",
    )

    parser.add_argument(
        "--extension",
        "-e",
        help=("File type (extension) of each conformation."),
        choices=["cif", "pdb", "json"],
        default="cif",
    )

    parser.add_argument(
        "--write_gt",
        help="Write/overwrite GT conformations after self-alignment.",
        action="store_true",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        help="Increase output verbosity",
        action="store_true",
    )


def get_name():
    return "calc_min_rmsd_metrics"


def align_min_rmsd(
    pred_traj: mda.Universe,
    GT_traj: mda.Universe,
    pred_traj_copy: mda.Universe,
    selection: str,
    verbose: bool = False,
) -> tuple[dict, mda.Universe]:

    # loop through each predicted structure and find GT structure with min RMSD
    rmsd_min: dict = {}
    rmsd_min["pred_index"] = []
    rmsd_min["rmsds"] = []
    rmsd_min["gt_index"] = []
    # pred_traj_aligned is not yet aligned and is at this point just a copy of pred_traj
    for i, pred in enumerate(
        zip(pred_traj.trajectory, pred_traj_copy.trajectory)
    ):
        rmsd_vals = []
        for GT in GT_traj.trajectory:
            # calc rmsd
            # align the pred model to the GT model
            # rms.rmsd docs: https://docs.mdanalysis.org/1.1.1/
            # documentation_pages/analysis/rms.html#module-MDAnalysis.analysis.rms
            # it aligns set of coords (arg1) to ref (arg2)
            # superpose implies center=True

            # here we are not excluding hist tag in rmsd calc or in min rmsd alignment
            rmsd_vals.append(
                rms.rmsd(
                    pred_traj.select_atoms(
                        selection,
                        # TODO decide on whether to use all (heavy) atoms for
                        # rmsd calculation
                    ).positions,
                    GT_traj.select_atoms(
                        selection,
                    ).positions,  # TODO update as above
                    # ensuring QCP superposition seems to make rms.rmsd and
                    # align.alignto np.close to each other
                    superposition=True,
                )
            )
        rmsd_vals_arr: np.ndarray[tuple] = np.array(rmsd_vals, dtype=float)
        min_index = np.argmin(rmsd_vals_arr)
        # NOTE: We store pairwise min-rmsds in rmsds array to dictionary.
        # May be useful for later analyses.
        # Casting to int/float as numpy dtypes not save-able as json.
        rmsd_min["pred_index"].append(int(i))
        # save the rmsd_min array and use gt_index (ie: min_index) to recover min-RMSD
        # for each pred conformation. This is useful for later analyses when
        # rmsd_min_dict is saved to file
        rmsd_min["rmsds"].append(rmsd_vals_arr.astype(float).tolist())
        rmsd_min["gt_index"].append(int(min_index))

    # now that we have the indices of the best rmsd fit for each pred_traj
    # we now align them in the copied pred_traj
    # and fill in the new atomic positions
    for pred_idx, gt_idx, rmsds in zip(
        rmsd_min["pred_index"],
        rmsd_min["gt_index"],
        rmsd_min["rmsds"],
    ):
        # change the frame of the pred_traj and pred_traj_copy to correct idx
        pred_traj.trajectory[pred_idx]
        pred_traj_copy.trajectory[pred_idx]
        # change the frame of GT_traj to correct idx
        GT_traj.trajectory[gt_idx]
        # do the alignment between the two frames using GT as ref
        rmsd_unaligned, rmsd_aligned = align.alignto(
            pred_traj,
            GT_traj,
            select=selection,
            # NOTE: Aligning using protein assuming H is
            # already removed!
        )
        if verbose:
            print(
                "{}: {} -> {}".format(
                    pred_idx,
                    rmsd_unaligned,
                    rmsd_aligned,
                )
            )
        # check the alignment is approx the same as reported by rms.rmsd above
        assert np.isclose(
            rmsds[gt_idx],
            rmsd_aligned,
        ), "rmsd_min:{}\nrmsd_aligned:{}".format(rmsds[gt_idx], rmsd_aligned)

        # input aligned atomic positions to copy of pred_traj
        # Have to do here because pred_traj will reset to unaligned after
        # iteration through loop
        pred_traj_copy.atoms.positions = copy.deepcopy(
            pred_traj.atoms.positions
        )

    return rmsd_min, pred_traj_copy


def calc_lfc(rmsf_gt, rmsf_pred):
    abs_diff_rmsf = np.abs(rmsf_gt - rmsf_pred)
    lfc = np.average(abs_diff_rmsf)
    lfc_stddev = np.std(abs_diff_rmsf)
    return lfc, lfc_stddev, abs_diff_rmsf


def make_plot_lfc(ad_rmsf, plot_dir: str):
    # TODO add docstring with info about non-H atom indices
    plt.plot(
        range(len(ad_rmsf)), ad_rmsf, color="blue", label="|$\\Delta$RMSF|"
    )
    plt.xlabel("Atom index")
    plt.ylabel("RMSF ($\\AA$)")
    plt.axvspan(79, 2175, zorder=0, alpha=0.2, color="orange", label="D1")
    plt.axvspan(4346, 4366, zorder=0, alpha=0.2, color="blue", label="LNK")
    plt.axvspan(6465, 8085, zorder=0, alpha=0.2, color="green", label="D4")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "rmsf.png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, "rmsf.pdf"))
    plt.clf()


def make_over_under_represent_het_plot(
    rmsf_gt,
    rmsf_pred,
    plot_dir: str,
    avg_over: int = 100,
):
    # +ve value of diff means heterogeneity under-estimated
    # -ve value of diff means heterogeneity over-estimated
    diff = rmsf_gt - rmsf_pred
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


def make_plot_rmsd_min(rmsd_min_vals, plot_dir: str):
    """Plot a histogram of the per-conformation min-RMSD values with a KDE overlay.

    Shows the distribution of how well each predicted conformation matches the
    ground-truth ensemble, complementing the single SAC (mean) value.

    Args:
        rmsd_min_vals: Array-like of min-RMSD values (in Angstroms) for each
            predicted conformation.
        plot_dir (str): Directory to save the plot in.
    """
    # TODO improve error bar plotting style
    rmsd_min_arr = np.asarray(rmsd_min_vals, dtype=float)
    # histogram
    plt.hist(
        rmsd_min_arr,
        bins="auto",
        density=True,
        alpha=0.6,
        color="steelblue",
        label="Histogram",
    )
    # KDE overlay
    kde = scipy.stats.gaussian_kde(rmsd_min_arr)
    x_vals = np.linspace(rmsd_min_arr.min(), rmsd_min_arr.max(), 200)
    plt.plot(
        x_vals,
        kde(x_vals),
        color="black",
        label="KDE",
    )
    # vertical line at the mean (= SAC)
    plt.axvline(
        np.mean(rmsd_min_arr),
        color="red",
        linestyle="--",
        label="Mean (SAC) = {:.2f} $\\AA$".format(np.mean(rmsd_min_arr)),
    )
    plt.xlabel("Min-RMSD to GT ($\\AA$)")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "rmsd_min.png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, "rmsd_min.pdf"))
    plt.clf()


def make_plot_rmsd_min_ecdf(rmsd_min_vals, plot_dir: str):
    """Plot the empirical cumulative distribution function (ECDF) of min-RMSD.

    Shows the fraction of predicted conformations whose min-RMSD to the
    ground-truth ensemble is at or below a given value, complementing the
    single SAC (mean) value and the A1 histogram/KDE plot.

    Args:
        rmsd_min_vals: Array-like of min-RMSD values (in Angstroms) for each
            predicted conformation.
        plot_dir (str): Directory to save the plot in.
    """
    rmsd_min_arr = np.asarray(rmsd_min_vals, dtype=float)
    sorted_vals = np.sort(rmsd_min_arr)
    # ECDF: fraction of values <= x, as a step function
    y_vals = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
    plt.step(
        sorted_vals, y_vals, where="post", color="steelblue", label="ECDF"
    )
    # reference lines at the median and the mean (= SAC)
    median = np.median(rmsd_min_arr)
    mean = np.mean(rmsd_min_arr)
    plt.axvline(
        median,
        color="green",
        linestyle="--",
        label="Median = {:.2f} $\\AA$".format(median),
    )
    plt.axvline(
        mean,
        color="red",
        linestyle="--",
        label="Mean (SAC) = {:.2f} $\\AA$".format(mean),
    )
    plt.xlabel("Min-RMSD to GT ($\\AA$)")
    plt.ylabel("Fraction of conformations $\\leq$ value")
    plt.ylim(0, 1.05)
    # grid with both major and minor ticks
    plt.grid(True, which="major", linestyle="--", alpha=0.5)
    plt.grid(True, which="minor", linestyle=":", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "rmsd_min_ecdf.png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, "rmsd_min_ecdf.pdf"))
    plt.clf()


def make_plot_rmsd_min_box_violin(rmsd_min_vals, plot_dir: str):
    """Plot a box-and-violin summary of the per-conformation min-RMSD values.

    Shows the distribution of how well each predicted conformation matches the
    ground-truth ensemble as a violin (kernel density) with an overlaid box
    plot (median, quartiles, whiskers, outliers), complementing the A1
    histogram/KDE and A2 ECDF plots.

    Args:
        rmsd_min_vals: Array-like of min-RMSD values (in Angstroms) for each
            predicted conformation.
        plot_dir (str): Directory to save the plot in.
    """
    rmsd_min_arr = np.asarray(rmsd_min_vals, dtype=float)
    """
    The KDE used in the violin plots is poor and does not represent the
    distribution well. The violin plot is not very informative and is not
    used in the final report.
    # violin plot (kernel density estimate of the distribution)
    parts = plt.violinplot(
        rmsd_min_arr,
        positions=[1],
        showmeans=True,
        showmedians=True,
        showextrema=True,
        bw_method="silverman",
    )
    for body in parts["bodies"]:
        body.set_facecolor("steelblue")
        body.set_alpha(0.6)
        body.set_edgecolor("black")
    """
    # box plot (median, quartiles, whiskers, outliers)
    plt.boxplot(
        rmsd_min_arr,
        positions=[1],
        widths=0.15,
        patch_artist=True,
        showmeans=True,
        meanline=True,
        showfliers=True,
        showbox=True,
        showcaps=True,
        boxprops=dict(facecolor="steelblue", edgecolor="black"),
        meanprops=dict(color="red", linewidth=2, linestyle="-", label="Mean"),
        medianprops=dict(
            color="blue", linewidth=2, linestyle="-", label="Median"
        ),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
        flierprops=dict(
            marker="x", markerfacecolor="black", markersize=2, label="Outliers"
        ),
    )
    """
    # horizontal line at the mean (= SAC)
    mean = np.mean(rmsd_min_arr)
    plt.axhline(
        mean,
        color="red",
        linestyle="--",
        label="Mean (SAC) = {:.2f} $\\AA$".format(mean),
    )
    median = np.median(rmsd_min_arr)
    plt.axhline(
        median,
        color="blue",
        linestyle="--",
        label="Median (SAC) = {:.2f} $\\AA$".format(median),
    )
    """
    plt.xlabel("Predicted ensemble")
    plt.ylabel("Min-RMSD to GT ($\\AA$)")
    plt.xticks([1], ["Pred"])
    # grid with both major and minor ticks
    plt.grid(True, which="major", linestyle="--", alpha=0.5, axis="y")
    plt.grid(True, which="minor", linestyle=":", alpha=0.3, axis="y")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "rmsd_min_box_violin.png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, "rmsd_min_box_violin.pdf"))
    plt.clf()


def make_plot_assignment(gt_index, plot_dir: str):
    """Plot which GT conformation each predicted conformation is assigned to.

    Each predicted conformation is assigned to the GT conformation with the
    minimum RMSD (the ``gt_index``). This plot shows that assignment as a
    scatter of predicted-conformation index against assigned GT-conformation
    index, so it is easy to see how the predicted conformations map onto the
    ground-truth ensemble.

    Args:
        gt_index: Array-like of GT-conformation indices (one per predicted
            conformation).
        plot_dir (str): Directory to save the plot in.
    """
    gt_index_arr = np.asarray(gt_index, dtype=int)
    pred_index = np.arange(len(gt_index_arr))
    plt.scatter(
        pred_index,
        gt_index_arr,
        color="steelblue",
        s=20,
        alpha=0.7,
        edgecolor="black",
        linewidth=0.5,
    )
    plt.xlabel("Predicted conformation index")
    plt.ylabel("Assigned GT conformation index")
    # grid with both major and minor ticks
    plt.grid(True, which="major", linestyle="--", alpha=0.5)
    plt.grid(True, which="minor", linestyle=":", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "assignment.png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, "assignment.pdf"))
    plt.clf()


def make_plot_occupancy(gt_index, plot_dir: str):
    """Plot how many predicted conformations occupy each GT conformation.

    Counts how many predicted conformations are assigned to each GT
    conformation (via ``gt_index``) and plots this as a histogram, so it is
    easy to see which GT conformations are well-sampled by the predicted
    ensemble and which are not sampled at all.

    Args:
        gt_index: Array-like of GT-conformation indices (one per predicted
            conformation).
        plot_dir (str): Directory to save the plot in.
    """
    gt_index_arr = np.asarray(gt_index, dtype=int)
    n_gt = gt_index_arr.max() + 1
    counts = np.bincount(gt_index_arr, minlength=n_gt)
    plt.bar(
        np.arange(n_gt),
        counts,
        color="steelblue",
        edgecolor="black",
        linewidth=0.5,
    )
    plt.xlabel("GT conformation index")
    plt.ylabel("Number of predicted conformations")
    # grid with both major and minor ticks
    plt.grid(True, which="major", linestyle="--", alpha=0.5, axis="y")
    plt.grid(True, which="minor", linestyle=":", alpha=0.3, axis="y")
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "occupancy.png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, "occupancy.pdf"))
    plt.clf()


def make_plot_coverage(gt_index, plot_dir: str):
    """Plot the cumulative coverage of the GT ensemble by the predicted one.

    Shows, as a function of the number of GT conformations considered, the
    fraction of GT conformations that have at least one predicted conformation
    assigned to them. A steep curve that reaches 1.0 quickly means the
    predicted ensemble covers the GT ensemble well.

    Args:
        gt_index: Array-like of GT-conformation indices (one per predicted
            conformation).
        plot_dir (str): Directory to save the plot in.
    """
    gt_index_arr = np.asarray(gt_index, dtype=int)
    n_gt = gt_index_arr.max() + 1
    # unique GT conformations that have at least one predicted conformation
    covered = np.unique(gt_index_arr)
    # cumulative fraction of GT conformations covered as we consider more of
    # the GT ensemble (ordered by GT index)
    covered_set = set(covered.tolist())
    cumulative = np.array(
        [
            sum(1 for j in range(i + 1) if j in covered_set) / n_gt
            for i in range(n_gt)
        ]
    )
    plt.plot(
        np.arange(n_gt),
        cumulative,
        color="steelblue",
        linewidth=2,
        marker="o",
        markersize=3,
    )
    plt.xlabel("Number of GT conformations considered")
    plt.ylabel("Fraction of GT conformations covered")
    plt.ylim(0, 1.05)
    # grid with both major and minor ticks
    plt.grid(True, which="major", linestyle="--", alpha=0.5)
    plt.grid(True, which="minor", linestyle=":", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "coverage.png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, "coverage.pdf"))
    plt.clf()


def create_selection_string(
    select_str: str,
    include_tag: bool,
    verbose: bool,
) -> str:
    selection = select_str
    if include_tag:
        selection = selection + " and resnum 10:1014"
    if verbose:
        print("Selection string is: '{}'".format(selection))
    return selection


def main(args):
    # create atom selection string used for alignment and RMSD/RMSF calc
    selection = create_selection_string(
        args.selection,
        args.exclude_hist_tag,
        args.verbose,
    )
    # STEP 1 - load predicted and ground truth ensembles
    # ground truth ensemble
    gt_ensemble = get_ensemble(
        args.truth_ensemble_dir, file_ext=args.extension, verbose=args.verbose
    )
    if args.verbose:
        for conf in gt_ensemble:
            print(os.path.basename(conf))

    # pred ensemble
    pred_ensemble = get_ensemble(
        args.pred_ensemble_dir, file_ext=args.extension, verbose=args.verbose
    )
    if args.verbose:
        for conf in pred_ensemble:
            print(os.path.basename(conf))

    # ------------------------------------------------------------------------------------
    # STEP 2 - create mda universes to hold the ensembles
    gt_u = mda.Universe(
        gt_ensemble[0],  # topology
        gt_ensemble,
        in_memory=True,
    )

    pred_u = mda.Universe(
        gt_ensemble[0],
        # using GT for topology
        pred_ensemble,
        in_memory=True,
    )

    # make a copy of pred_u which we later fill with the conformations once they are
    # optimally aligned to ground truth ensemble
    pred_u_copy = mda.Universe(
        gt_ensemble[0],
        # using GT for topology
        pred_ensemble,
        in_memory=True,
    )

    # ------------------------------------------------------------------------------------
    # STEP 3 - self-align GT ensemble
    # NOTE: gt universe should already be globally aligned...
    # ie: aligner = align.AlignTraj(u, u, select="backbone",
    # in_memory=True).run() was applied to all 2000 conformations before
    # subsampling for dataset creation.
    # Re-aligning may change the alignment slightly but is done to allow
    # calculation of RMSF metrics against a truth ensemble aligned against
    # itself. This is done below

    # align to reference
    # NOTE: we do NOT remove hist-tag from ground-truth alignment by default as this
    # may slightly affect the GT RMSF distribution
    # It does not affect RMSDs calculated between GT and pred confs though
    # This should not affect individual min-RMSD alignments of pred
    # conformations to the aligned GT conformations
    # NOTE: we select all protein rather than backbone only
    # The hydrogens (and any other non-heavy atoms) must have) already been
    # removed!
    aligner_gt_u = align.AlignTraj(  # noqa: F841
        gt_u,
        gt_u,
        select=selection,  # and resnum 10:1014',
    ).run()
    # NOTE: use in_memory so the alignments persist!!!!!

    # NOTE: if using mass-weighting in RMSD/RMSF calculations, must also use
    # it to do the GT alignment (as this is via min RMSD!)
    # We DO NOT currently use mass-weighting.
    # If using mass-weighting keep this alignment AND add weight="mass"
    # as an argument!)
    # NOTE: If using weighting, SAVE GT mass/charge weighted aligned
    # conformations to file (mmcif) for reuse in PCV calculation

    # NOTE: we don't align pred_ensemble to itself.
    # It is only aligned on a per-conf basis to the (now aligned) GT confs

    # ------------------------------------------------------------------------------------
    # STEP 4 - min-RMSD align each conf of pred_ensemble to closest conf in GT ensemble
    # and get new min-RMSD-aligned pred_ensemble
    # dictionary holds min-RMSD value and index of min-RMSD GT conformation for each
    # predicted conformation
    # NOTE: Not using mass-weighting as RMSF calculation does not support this
    pred_u_rmsd_min_dict, pred_u_rmsd_min_aligned = align_min_rmsd(
        pred_u,
        gt_u,
        pred_u_copy,
        selection=selection,
    )
    if args.verbose:
        for key, val in pred_u_rmsd_min_dict.items():
            print("{}:{}".format(key, val))

    # ------------------------------------------------------------------------------------
    # STEP 5 - create directories to save aligned atomic models as pdb/xtc.
    # NOTE: Currently saving these aligned ensembles as dir of pdbs and a xtc trajectory
    # as writing mmcif is not supported by MDAnalysis from a mda.Universe
    # Also creating a directory which can be used as output to convert pdbs to
    # mmcifs via cahra-metric convert_ensemble <args>

    # create the directory structure and results object
    results = Ensemble_Stats(
        args.output_dir,
        args.pred_ensemble_label,
        verbose=args.verbose,
    )
    # create results directories
    results.make_results_dir_structure()

    # NOTE: SAVE self-aligned GT to disk in a new directory
    # (to be used as input for min-RMSD metric recalculation or PCV metric calculation)
    # gt_u is self-aligned GT ensemble mda.Universe
    # write xtc
    # TODO fix writing all conformations to xtc
    gt_u.atoms.write(
        os.path.join(
            results.gt_confs_dir_xtc, "{}.xtc".format(args.pred_ensemble_label)
        ),
        overwrite=args.write_gt,
    )
    # write pdbs
    gt_ensemble_filenames = []
    # loop through traj frames and write each one to file
    for gt_filename, gt_frame in zip(gt_ensemble, gt_u.trajectory):
        if args.verbose:
            print("Saving GT frame {} to pdb".format(gt_frame.frame))
        # keep same file labels but convert to pdb
        new_ensemble_filename = replace_file_ext(
            os.path.basename(gt_filename),
            "pdb",
            results.gt_confs_dir_pdbs,
            args.verbose,
        )
        # write pdbs in output directory
        gt_ensemble_filenames.append(new_ensemble_filename)
        if not os.path.isfile(new_ensemble_filename) or args.write_gt:
            print("Writing out self-aligned GT models")
            gt_u.atoms.write(
                new_ensemble_filename,
                # overwrite=args.write_gt,
            )
    # no way to write mmcifs from mda.Universe. Could use conversion tool from pdb

    # NOTE: SAVE min-RMSD-aligned pred models to disk in a new directory
    # pred_u_rmsd_min_aligned is min-RMSD aligned pred ensemble mda.Universe
    # write xtc
    # TODO fix writing all conformations to xtc
    pred_u_rmsd_min_aligned.atoms.write(
        os.path.join(
            results.pred_confs_dir_label_xtc,
            "{}.xtc".format(args.pred_ensemble_label),
        )
    )
    # write pdbs
    pred_ensemble_filenames = []
    for pred_filename, pred_u_rmsd_min_aligned_frame in zip(
        pred_ensemble, pred_u_rmsd_min_aligned.trajectory
    ):
        if args.verbose:
            print(
                "Saving pred frame {} to pdb".format(
                    pred_u_rmsd_min_aligned_frame.frame
                )
            )
        # keep same file labels but convert to pdb
        new_ensemble_filename = replace_file_ext(
            os.path.basename(pred_filename),
            "pdb",
            results.pred_confs_dir_label_pdbs,
            args.verbose,
        )
        # write pdbs in output directory
        pred_ensemble_filenames.append(new_ensemble_filename)
        pred_u_rmsd_min_aligned.atoms.write(
            new_ensemble_filename,
            # overwrite=True,
        )
    # no way to write mmcifs from mda.Universe. Could use conversion tool from pdb

    # add gt and pred saved filenames to results object now that the aligned models
    # are saved
    results.gt_filenames = gt_ensemble_filenames
    results.pred_filenames = pred_ensemble_filenames
    # also add these filenames to the rmsd_min dict as a list
    pred_u_rmsd_min_dict["pred_filename"] = results.pred_filenames
    # write this dict to the results object
    results.ensemble = {}
    results.ensemble["conf_stats"] = pred_u_rmsd_min_dict
    results.ensemble["ensemble_stats"] = {}
    results.ensemble["ensemble_stats"][
        "pred_ensemble_label"
    ] = args.pred_ensemble_label

    # ------------------------------------------------------------------------------------
    # STEP 6 - get SAC

    SAC = float(
        np.average(
            extract_min_rmsd_values(
                pred_u_rmsd_min_dict["rmsds"],
                pred_u_rmsd_min_dict["gt_index"],
            )
        )
    )
    results.ensemble["ensemble_stats"]["sac"] = SAC
    print("SAC: {}".format(SAC))
    # TODO Save pred_u_rmsd_min_dict and SAC/LFC values to file (use dict for latter?)
    # also could save atom-wise RMSF values for aligned-pred-ensemble (AND GT ensemble,
    # to check consistency of alignment/realignmnet if that is kept in the above code)

    # ------------------------------------------------------------------------------------
    # STEP 7 - calc RMSF on per-residue basis to allow calculation of LFC
    # NOTE: Not mass-weighting this currently. Currently just give all atoms same
    # weighting. The hydrogens should not be present in the submitted (or GT) mmcif
    # files but ensuring this is the case should be done by validate_ensembles code
    # It appears weights can be applied to an atomgroup object, ie provide rms.RMSF with
    # the argument 'weights="mass"', in similar way to RMSD calc above
    # see: https://docs.mdanalysis.org/stable/documentation_pages/core/groups.html#
    # MDAnalysis.core.groups.AtomGroup
    # HOWEVER, rms.RMSF reports that it does not use mass weighting...
    # see https://www.mdanalysis.org/pmda/api/rmsf.html
    # If this is the case, then do not mass/charge weight anything unless you are
    # sure it can be used in all atomic metrics. Do state that RMSDs are NOT MASS
    # WEIGHTED if you do not mass-weight them

    # so here get the RMSF of our aligned GT ensemble for each atom
    # TODO determine option to calc on per-atom basis
    # TODO choose whether to calc per-residue or per-atom RMSF
    rmsf_gt = rms.RMSF(gt_u.select_atoms(selection)).run()
    # and for the min-RMSD-aligned pred ensemble
    rmsf_pred = rms.RMSF(pred_u_rmsd_min_aligned.select_atoms(selection)).run()
    results.ensemble["ensemble_stats"][
        "gt_rmsf"
    ] = rmsf_gt.results.rmsf.tolist()
    results.ensemble["ensemble_stats"][
        "pred_rmsf"
    ] = rmsf_pred.results.rmsf.tolist()
    # ------------------------------------------------------------------------------------
    # STEP 8 - get LFC
    # make and save LFC plots
    lfc_all, lfc_stddev, abs_diff_rmsf = calc_lfc(
        rmsf_gt.results.rmsf,
        rmsf_pred.results.rmsf,
    )
    results.ensemble["ensemble_stats"]["lfc"] = lfc_all
    results.ensemble["ensemble_stats"]["lfc_abs_diff"] = abs_diff_rmsf.tolist()
    if args.verbose:
        print("LFC:{}\nLFC_stddev:{}".format(lfc_all, lfc_stddev))

    # save results dictionary to json
    save_dict_to_json(
        results.__dict__,
        os.path.join(results.pred_confs_dir_label_results, "rms_metrics.json"),
    )
    make_plot_lfc(
        abs_diff_rmsf,
        results.pred_confs_dir_label_figures,
    )
    make_over_under_represent_het_plot(
        rmsf_gt.results.rmsf,
        rmsf_pred.results.rmsf,
        results.pred_confs_dir_label_figures,
    )
    make_plot_rmsd_min(
        extract_min_rmsd_values(
            pred_u_rmsd_min_dict["rmsds"],
            pred_u_rmsd_min_dict["gt_index"],
        ),
        results.pred_confs_dir_label_figures,
    )
    make_plot_rmsd_min_ecdf(
        extract_min_rmsd_values(
            pred_u_rmsd_min_dict["rmsds"],
            pred_u_rmsd_min_dict["gt_index"],
        ),
        results.pred_confs_dir_label_figures,
    )
    make_plot_rmsd_min_box_violin(
        extract_min_rmsd_values(
            pred_u_rmsd_min_dict["rmsds"],
            pred_u_rmsd_min_dict["gt_index"],
        ),
        results.pred_confs_dir_label_figures,
    )
    make_plot_assignment(
        pred_u_rmsd_min_dict["gt_index"],
        results.pred_confs_dir_label_figures,
    )
    make_plot_occupancy(
        pred_u_rmsd_min_dict["gt_index"],
        results.pred_confs_dir_label_figures,
    )
    make_plot_coverage(
        pred_u_rmsd_min_dict["gt_index"],
        results.pred_confs_dir_label_figures,
    )

    # TODO option to also visualise/save rmsf via temp factor in (saved) min-RMSD
    # pred-ensemble mda universe
    # see https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/
    # rmsf.html#Visualising-RMSF-as-B-factors
    # NOTE I have done this previously - example in extract_frames_6xou_unsteered_0-
    # 200ns_alignment_metric_tests.ipynb
    # (ie: the notebook which showed volume-ensemble-metrics weren't invariant)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    if args.verbose:
        for arg in vars(args):
            print("{}, {}".format(arg, getattr(args, arg)))
    main(args)
