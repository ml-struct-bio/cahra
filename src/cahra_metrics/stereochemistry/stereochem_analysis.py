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

# import csv
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
import numpy as np

from cahra_metrics.utilities.utils import (
    # load_dict_from_json,
    get_stereochem_json_files,
    Ensemble_Stats,
    save_dict_to_json,
)


def add_arguments(parser: argparse.ArgumentParser):
    """Set up arguments for parsing

    Args:
        parser (argparse.ArgumentParser): Waymarking argument parser

    Returns:
        argparse.ArgumentParser:  Waymarking argument parser with arguments
    """

    parser.add_argument(
        "--results_dir",
        "-o",
        help=(
            "Path to directory holding the min_rmsd and"
            " calc_pcv results calculated for all gt and pred ensembles."
        ),
        type=str,
        default="",
    )

    parser.add_argument(
        "--stereochem_dir",
        "-s",
        help=(
            "Path to directory holding rota-rama-clash files which were"
            " calculated for an ensemble."
        ),
        type=str,
        default="",
    )

    parser.add_argument(
        "--gt_subdir",
        "-gts",
        help=(
            "Relative path to subdirectory holding rota-rama-clash files"
            " which were calculated for gt ensemble. Should be same for"
            " results_dir and sterochem_dir"
        ),
        type=str,
        default="gt_confs",
    )

    parser.add_argument(
        "--orig_subdir",
        "-os",
        help=(
            "Relative path to subdirectory holding rota-rama-clash files"
            " which were calculated for pdb structure. Should be same for"
            " results_dir and sterochem_dir"
        ),
        type=str,
        default="orig_6xou",
    )

    parser.add_argument(
        "--pred_subdir",
        "-ps",
        help=(
            "Relative path to subdirectory holding rota-rama-clash files"
            " which were calculated for pred ensembles. Should be same for"
            " results_dir and sterochem_dir"
        ),
        type=str,
        default="pred_confs",
    )

    parser.add_argument(
        "--pred_ensemble_label",
        "-pl",
        help=(
            "Label for predicted ensemble. Used to set input and output"
            " directories. Use one value per predicted ensemble directory."
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


def check_directories_exist(args):
    """Check that all required directories exist"""
    stereochem_dir = args.stereochem_dir
    results_dir = args.results_dir

    gt_subdir = args.gt_subdir
    orig_subdir = args.orig_subdir
    pred_subdir = args.pred_subdir
    pred_ensemble_label = args.pred_ensemble_label

    assert os.path.isdir(stereochem_dir), "{} does not exist".format(
        stereochem_dir
    )
    assert os.path.isdir(results_dir), "{} does not exist".format(results_dir)

    # check gt, pdb and pred ensemb subdirectories exist
    # for stereochem dir
    dir_to_check = os.path.join(stereochem_dir, gt_subdir)
    assert os.path.isdir(
        dir_to_check,
    ), "{} does not exist".format(dir_to_check)
    dir_to_check = os.path.join(stereochem_dir, orig_subdir)
    assert os.path.isdir(
        dir_to_check,
    ), "{} does not exist".format(dir_to_check)
    for label in pred_ensemble_label:
        dir_to_check = os.path.join(stereochem_dir, pred_subdir, label)
        assert os.path.isdir(
            dir_to_check,
        ), "{} does not exist".format(dir_to_check)

    # checking gt and pred ensemb subdirectories exist
    # for results_dir
    # not checking results_dir/orig_6xou as manually created
    dir_to_check = os.path.join(stereochem_dir, gt_subdir)
    assert os.path.isdir(
        dir_to_check,
    ), "{} does not exist".format(dir_to_check)
    for label in pred_ensemble_label:
        dir_to_check = os.path.join(stereochem_dir, pred_subdir, label)
        assert os.path.isdir(
            dir_to_check,
        ), "{} does not exist".format(dir_to_check)
    return


def check_basenames_identical(jsonpath: str, pdbpath: str):
    json_base = os.path.splitext(os.path.basename(jsonpath))[0]
    pdb_base = os.path.splitext(os.path.basename(pdbpath))[0]
    assert pdb_base in json_base, "{} not matching {}".format(
        json_base,
        pdb_base,
    )
    return


def extract_summary_rota(geom_jsonfile: str):
    with open(geom_jsonfile, "r") as f:
        results = json.load(f)
        summary_results = results["summary_results"][""]
    return summary_results


def extract_summary_rama(geom_jsonfile: str):
    with open(geom_jsonfile, "r") as f:
        results = json.load(f)
        summary_results = results["summary_results"][""]
    return summary_results


def extract_summary_clash(geom_jsonfile: str):
    with open(geom_jsonfile, "r") as f:
        results = json.load(f)
        summary_results = results["summary_results"][""]
    return summary_results


def add_geom_to_ensemble_stats(
    metrics_jsonfile: str,
    rota_jsonfiles: list[str],
    rama_jsonfiles: list[str],
    clash_jsonfiles: list[str],
):
    """Takes an ensemble (Ensemble_Stats class) and lists of
    molprobity calc'ed geometry files for each model in ensemble.
    Calculates averages over ensemble and saves updated
    all_geom_metrics.json file with per model molprobity summary
    info and averages of these values over the ensemble.

    Args:
        metrics_jsonfile (str): _description_
        rota_jsonfiles (list[str]): _description_
        rama_jsonfiles (list[str]): _description_
        clash_jsonfiles (list[str]): _description_
    """
    # Open all_metrics json
    with open(metrics_jsonfile, "r") as f:
        results = json.load(f)

    # Extract rotamer analysis scores for each model
    rota_num_favored = []
    rota_num_allowed = []
    rota_num_outliers = []
    rota_num_residues = []
    rota_outlier_percentage = []
    # loop whilst check natsort order correct
    for (
        jsonfile,
        pred_filename,
    ) in zip(rota_jsonfiles, results.pred_filename):
        check_basenames_identical(
            jsonfile,
            pred_filename,
        )
        summary_results = extract_summary_rota(jsonfile)
        rota_num_favored.append(summary_results["num_favored"])
        rota_num_allowed.append(summary_results["num_allowed"])
        rota_num_outliers.append(summary_results["num_outliers"])
        rota_num_residues.append(summary_results["num_residues"])
        rota_outlier_percentage.append(summary_results["outlier_percentage"])
    # Calculate average values for each
    rota_num_favored_avg = float(np.average(rota_num_favored))
    rota_num_allowed_avg = float(np.average(rota_num_allowed))
    rota_num_outliers_avg = float(np.average(rota_num_outliers))
    rota_num_residues_avg = float(np.average(rota_num_residues))
    rota_outlier_percentage_avg = float(np.average(rota_outlier_percentage))
    # Update json
    # Save model values for each pred model (to json)
    results.conf_stats["rota_num_favored"] = rota_num_favored
    results.conf_stats["rota_num_allowed"] = rota_num_allowed
    results.conf_stats["rota_num_outliers"] = rota_num_outliers
    results.conf_stats["rota_num_residues"] = rota_num_residues
    results.conf_stats["rota_outlier_percentage"] = rota_outlier_percentage
    # Save average values for each ensemble (to json)
    results.ensemble_stats["rota_num_favored_avg"] = rota_num_favored_avg
    results.ensemble_stats["rota_num_allowed_avg"] = rota_num_allowed_avg
    results.ensemble_stats["rota_num_outliers_avg"] = rota_num_outliers_avg
    results.ensemble_stats["rota_num_residues_avg"] = rota_num_residues_avg
    results.ensemble_stats["rota_outlier_percentage_avg"] = (
        rota_outlier_percentage_avg
    )

    # Extract rama analysis scores for each model
    rama_num_favored = []
    rama_num_allowed = []
    rama_num_outliers = []
    rama_num_residues = []
    rama_outlier_percentage = []
    rama_favored_percentage = []
    # loop whilst check natsort order correct
    for (
        jsonfile,
        pred_filename,
    ) in zip(rama_jsonfiles, results.pred_filename):
        check_basenames_identical(
            jsonfile,
            pred_filename,
        )
        summary_results = extract_summary_rama(jsonfile)
        rama_num_favored.append(summary_results["num_favored"])
        rama_num_allowed.append(summary_results["num_allowed"])
        rama_num_outliers.append(summary_results["num_outliers"])
        rama_num_residues.append(summary_results["num_residues"])
        rama_outlier_percentage.append(summary_results["outlier_percentage"])
        rama_favored_percentage.append(summary_results["favored_percentage"])
    # Calculate average values for each
    rama_num_favored_avg = float(np.average(rama_num_favored))
    rama_num_allowed_avg = float(np.average(rama_num_allowed))
    rama_num_outliers_avg = float(np.average(rama_num_outliers))
    rama_num_residues_avg = float(np.average(rama_num_residues))
    rama_outlier_percentage_avg = float(np.average(rama_outlier_percentage))
    rama_favored_percentage_avg = float(np.average(rama_favored_percentage))
    # Update json
    # Save model values for each pred model (to json)
    results.conf_stats["rama_num_favored"] = rama_num_favored
    results.conf_stats["rama_num_allowed"] = rama_num_allowed
    results.conf_stats["rama_num_outliers"] = rama_num_outliers
    results.conf_stats["rama_num_residues"] = rama_num_residues
    results.conf_stats["rama_outlier_percentage"] = rama_outlier_percentage
    results.conf_stats["rama_favored_percentage"] = rama_favored_percentage
    # Save average values for each ensemble (to json)
    results.ensemble_stats["rama_num_favored_avg"] = rama_num_favored_avg
    results.ensemble_stats["rama_num_allowed_avg"] = rama_num_allowed_avg
    results.ensemble_stats["rama_num_outliers_avg"] = rama_num_outliers_avg
    results.ensemble_stats["rama_num_residues_avg"] = rama_num_residues_avg
    results.ensemble_stats["rama_outlier_percentage_avg"] = (
        rama_outlier_percentage_avg
    )
    results.ensemble_stats["rama_favored_percentage_avg"] = (
        rama_favored_percentage_avg
    )

    # Extract clashscore for each model
    clashscore = []
    num_clashes = []
    for (
        jsonfile,
        pred_filename,
    ) in zip(clash_jsonfiles, results.pred_filename):
        check_basenames_identical(
            jsonfile,
            pred_filename,
        )
        summary_results = extract_summary_clash(jsonfile)
        clashscore.append(summary_results["clashscore"])
        num_clashes.append(summary_results["num_clashes"])
    # Calculate average values for each
    clashscore_avg = float(np.average(clashscore))
    num_clashes_avg = float(np.average(num_clashes))
    # Update json
    # Save model values for each pred model (to json)
    results.conf_stats["clashscore"] = clashscore
    results.conf_stats["num_clashes"] = num_clashes
    # Save average values for each ensemble (to json)
    results.ensemble_stats["clashscore_avg"] = clashscore_avg
    results.ensemble_stats["num_clashes_avg"] = num_clashes_avg

    # Save updated json with updated name
    updated_json_file = os.path.join(
        results.pred_confs_dir_label_results,
        "all_geom_metrics.json",
    )
    save_dict_to_json(results, updated_json_file)
    return


def plot_hist(vals, xlabel: str, plot_dir: str, plot_name: str):
    """_summary_

    Args:
        vals (_type_): _description_
        xlabel (str): _description_
        plot_dir (str): _description_
        plot_name (str): _description_
    """
    vals = np.asarray(vals, dtype=float)
    # histogram
    plt.hist(
        vals,
        bins="auto",
        density=True,
        alpha=0.6,
        color="steelblue",
        label="Histogram",
    )
    # KDE overlay
    kde = scipy.stats.gaussian_kde(vals)
    x_vals = np.linspace(vals.min(), vals.max(), 200)
    plt.plot(
        x_vals,
        kde(x_vals),
        color="black",
        label="KDE",
    )
    # vertical line at the mean
    plt.axvline(
        np.mean(vals),
        color="red",
        linestyle="--",
        label="Mean",
    )
    plt.xlabel(xlabel)
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, plot_name + ".png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, plot_name + ".pdf"))
    plt.clf()


def plot_scatter(
    vals_x,
    vals_y,
    xlabel: str,
    ylabel: str,
    plot_dir: str,
    plot_name: str,
):

    vals_x = np.asarray(vals_x, dtype=float)
    vals_y = np.asarray(vals_y, dtype=float)
    # 2d scatterplot of CV angle (x-axis) and assigned GT conformation index
    # (y-axis) with no lines connecting points:
    plt.scatter(
        vals_x,
        vals_y,
        s=50,
        c="steelblue",
        edgecolors="black",
        linewidth=0.5,
    )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # grid with both major and minor ticks
    plt.grid(True, which="major", linestyle="--", alpha=0.5)
    plt.grid(True, which="minor", linestyle=":", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, plot_name + ".png"), dpi=300)
    plt.savefig(os.path.join(plot_dir, plot_name + ".pdf"))
    plt.clf()


def plot_violin(
    labels: list[str],
    val_lists: list[list[float]],
    xlabel: str,
    ylabel: str,
    hline: float,
    hline_label: str,
    plot_dir: str,
    plot_name: str,
):

    # create a fresh figure sized for the number of datasets
    plt.figure(figsize=(max(8, 0.8 * len(labels)), 6))
    # violin plot of the min-RMSD distribution for each dataset
    parts = plt.violinplot(
        val_lists,
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
        val_lists,
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
    # hline
    plt.axvline(
        hline,
        color="red",
        linestyle="--",
        label=hline_label,
    )
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


def main(args):
    # STEP 1 - check that all required directories exist
    check_directories_exist(args)

    # ------------------------------------------------------------------------------------
    # STEP 2 - get lists of json files with sterochemical info

    # get lists of gt rota, rama and clash json files
    print("Grabbing gt json files")
    gt_rota_files = get_stereochem_json_files(
        os.path.join(
            args.stereochem_dir,
            args.gt_subdir,
            "pdb",
        ),
        pattern_match="rota.json",
        verbose=args.verbose,
    )
    gt_rama_files = get_stereochem_json_files(
        os.path.join(
            args.stereochem_dir,
            args.gt_subdir,
            "pdb",
        ),
        pattern_match="rama.json",
        verbose=args.verbose,
    )
    gt_clash_files = get_stereochem_json_files(
        os.path.join(
            args.stereochem_dir,
            args.gt_subdir,
            "pdb",
        ),
        pattern_match="clash.json",
        verbose=args.verbose,
    )
    if args.verbose:
        print("gt_rota_files:{}".format(gt_rota_files))
        print("gt_rama_files:{}".format(gt_rama_files))
        print("gt_clash_files:{}".format(gt_clash_files))

    # get lists of pred ensemble rota, rama and clash json files
    # one list per label
    pred_stereo_rama_files = {}
    pred_stereo_rota_files = {}
    pred_stereo_clash_files = {}
    for label in args.pred_ensemble_label:
        print("Grabbing {} json files".format(label))
        pred_stereo_rota_files[label] = get_stereochem_json_files(
            os.path.join(
                args.stereochem_dir,
                args.pred_subdir,
                label,
                "pdb",
            ),
            pattern_match="rota.json",
            verbose=args.verbose,
        )
        pred_stereo_rama_files[label] = get_stereochem_json_files(
            os.path.join(
                args.stereochem_dir,
                args.pred_subdir,
                label,
                "pdb",
            ),
            pattern_match="rama.json",
            verbose=args.verbose,
        )
        pred_stereo_clash_files[label] = get_stereochem_json_files(
            os.path.join(
                args.stereochem_dir,
                args.pred_subdir,
                label,
                "pdb",
            ),
            pattern_match="clash.json",
            verbose=args.verbose,
        )
        if args.verbose:
            print(
                "{}_rota_files:{}".format(label, pred_stereo_rota_files[label])
            )
            print(
                "{}_rama_files:{}".format(label, pred_stereo_rama_files[label])
            )
            print(
                "{}_clash_files:{}".format(
                    label, pred_stereo_clash_files[label]
                )
            )

    # get (single value) lists of rota, rama and clash json file for 6xou
    print("Grabbing orig_6xou json files")
    orig_6xou_rota_files = get_stereochem_json_files(
        os.path.join(
            args.stereochem_dir,
            args.orig_subdir,
            "pdb",
        ),
        pattern_match="rota.json",
        verbose=args.verbose,
    )
    orig_6xou_rama_files = get_stereochem_json_files(
        os.path.join(
            args.stereochem_dir,
            args.orig_subdir,
            "pdb",
        ),
        pattern_match="rama.json",
        verbose=args.verbose,
    )
    orig_6xou_clash_files = get_stereochem_json_files(
        os.path.join(
            args.stereochem_dir,
            args.orig_subdir,
            "pdb",
        ),
        pattern_match="clash.json",
        verbose=args.verbose,
    )
    if args.verbose:
        print("orig_6xou_rota_files:{}".format(orig_6xou_rota_files))
        print("orig_6xou_rama_files:{}".format(orig_6xou_rama_files))
        print("orig_6xou_clash_files:{}".format(orig_6xou_clash_files))
    # ------------------------------------------------------------------------------------
    # STEP 3 - get json results files (post calc_pcv) for each pred ensemble
    metrics_jsons = {}
    for label in args.pred_ensemble_label:
        # get fpath for all_metrics.json (includes pcv)
        fpath = os.path.join(
            arg.results_dir, args.pred_dir, label, "results/all_metrics.json"
        )
        metrics_jsons[label] = fpath

    # ------------------------------------------------------------------------------------
    # STEP 4 - make simplified stereochem jsons for gt and orig_6xou calculations
    # Note that the pred ensemble stereochemical calculations are added to
    # the pred ensemble json files
    # So here we write these as they don't exist previously for GT or 6xou
    gt_ens = Ensemble_Stats(
        output_dir=args.results_dir,
        pred_ensemble_label="gt",
        verbose=args.verbose,
    )
    gt_ens_jsonfile = os.path.join(
        gt_ens.pred_confs_dir_label_results,
        "all_metrics.json",
    )
    save_dict_to_json(
        gt_ens.__dict__,
        gt_ens_jsonfile,
    )

    orig_6xou_ens = Ensemble_Stats(
        output_dir=args.results_dir,
        pred_ensemble_label="6xou",
        verbose=args.verbose,
    )
    orig_6xou_ens_jsonfile = os.path.join(
        orig_6xou_ens.pred_confs_dir_label_results,
        "all_metrics.json",
    )
    save_dict_to_json(
        orig_6xou_ens.__dict__,
        orig_6xou_ens_jsonfile,
    )

    # ------------------------------------------------------------------------------------
    # STEP 5 - EXTRACT, CALCULATE AVERAGES AND MAKE PLOTS FOR GEOM ANALYSIS

    # Per pred ensemble
    ens_geom_lists = {
        "rota_outlier_percentage": [],
        "rama_outlier_percentage": [],
        "rama_favored_percentage": [],
        "clashscore": [],
        "labels": [],
    }
    for label in args.pred_ensemble_label:
        # Add molprobity summary info to json
        add_geom_to_ensemble_stats(
            metrics_jsons[label],
            pred_stereo_rota_files[label],
            pred_stereo_rama_files[label],
            pred_stereo_clash_files[label],
        )
        geom_jsonfile = os.path.join(
            os.path.dirname(metrics_jsons[label]),
            "all_geom_metrics.json",
        )
        with open(geom_jsonfile, "r") as f:
            results = json.load(f)

        # keep track to make comparative violinplot
        ens_geom_lists.rota_outlier_percentage.append(
            results.conf_stats.rota_outlier_percentage
        )
        ens_geom_lists.rama_outlier_percentage.append(
            results.conf_stats.rama_outlier_percentage
        )
        ens_geom_lists.rama_favored_percentage.append(
            results.conf_stats.rama_favored_percentage
        )
        ens_geom_lists.clashscore.append(results.conf_stats.clashscore)
        ens_geom_lists.labels.append(label)
        # Make 2d hists
        # rota_outlier_percentage
        plot_hist(
            results.conf_stats.rota_outlier_percentage,
            "Rotamer outlier percentage",
            results.pred_confs_dir_label_figures,
            "rota_outlier_percentage",
        )
        # rama_outlier_percentage
        plot_hist(
            results.conf_stats.rama_outlier_percentage,
            "Ramachandran outlier percentage",
            results.pred_confs_dir_label_figures,
            "rama_outlier_percentage",
        )
        # rama_favored percentage
        plot_hist(
            results.conf_stats.rama_favored_percentage,
            "Ramachandran favored percentage",
            results.pred_confs_dir_label_figures,
            "rama_favored_percentage",
        )
        # clashscore
        plot_hist(
            results.conf_stats.clashscore,
            "Clash Score",
            results.pred_confs_dir_label_figures,
            "clashscore",
        )

    # GT ensemble
    # Add molprobity summary info to json
    add_geom_to_ensemble_stats(
        gt_ens_jsonfile,
        gt_rota_files,
        gt_rama_files,
        gt_clash_files,
    )

    geom_jsonfile = os.path.join(
        os.path.dirname(gt_ens_jsonfile),
        "all_geom_metrics.json",
    )
    with open(geom_jsonfile, "r") as f:
        results = json.load(f)
    # keep track to make comparative violinplot
    ens_geom_lists.rota_outlier_percentage.append(
        results.conf_stats.rota_outlier_percentage
    )
    ens_geom_lists.rama_outlier_percentage.append(
        results.conf_stats.rama_outlier_percentage
    )
    ens_geom_lists.rama_favored_percentage.append(
        results.conf_stats.rama_favored_percentage
    )
    ens_geom_lists.clashscore.append(results.conf_stats.clashscore)
    ens_geom_lists.labels.append("GT")
    # Make 2d hists
    # rota_outlier_percentage
    plot_hist(
        results.conf_stats.rota_outlier_percentage,
        "Rotamer outlier percentage",
        results.pred_confs_dir_label_figures,
        "rota_outlier_percentage",
    )
    # rama_outlier_percentage
    plot_hist(
        results.conf_stats.rama_outlier_percentage,
        "Ramachandran outlier percentage",
        results.pred_confs_dir_label_figures,
        "rama_outlier_percentage",
    )
    # rama_favored percentage
    plot_hist(
        results.conf_stats.rama_favored_percentage,
        "Ramachandran favored percentage",
        results.pred_confs_dir_label_figures,
        "rama_favored_percentage",
    )
    # clashscore
    plot_hist(
        results.conf_stats.clashscore,
        "Clash Score",
        results.pred_confs_dir_label_figures,
        "clashscore",
    )

    # Make scatterplots of gt model rama, rota and clash scores
    # which needs to be ordered by time (rather than CV)
    # rota_outlier_percentage
    plot_scatter(
        np.arange(len(results.conf_stats.rota_outlier_percentage)),
        results.conf_stats.rota_outlier_percentage,
        "GT model index",
        "Rotamer outlier percentage",
        results.pred_confs_dir_label_figures,
        "rota_outlier_percentage_progression",
    )
    # rama_outlier_percentage
    plot_scatter(
        np.arange(len(results.conf_stats.rama_outlier_percentage)),
        results.conf_stats.rama_outlier_percentage,
        "GT model index",
        "Ramachandran outlier percentage",
        results.pred_confs_dir_label_figures,
        "rama_outlier_percentage_progression",
    )
    # rama_favored percentage
    plot_scatter(
        np.arange(len(results.conf_stats.rama_favored_percentage)),
        results.conf_stats.rama_favored_percentage,
        "GT model index",
        "Ramachandran favored percentage",
        results.pred_confs_dir_label_figures,
        "rama_favored_percentage_progression",
    )
    # clashscore
    plot_scatter(
        np.arange(len(results.conf_stats.clashscore)),
        results.conf_stats.clashscore,
        "GT model index",
        "Clash Score",
        results.pred_confs_dir_label_figures,
        "clashscore_progression",
    )

    # orig_6xou model
    # Add molprobity summary info to json
    add_geom_to_ensemble_stats(
        orig_6xou_ens_jsonfile,
        orig_6xou_rota_files,
        orig_6xou_rama_files,
        orig_6xou_clash_files,
    )

    geom_jsonfile = os.path.join(
        os.path.dirname(orig_6xou_ens_jsonfile),
        "all_geom_metrics.json",
    )
    with open(geom_jsonfile, "r") as f:
        results = json.load(f)
    # Make 2d hists
    # rota_outlier_percentage
    plot_hist(
        results.conf_stats.rota_outlier_percentage,
        "Rotamer outlier percentage",
        results.pred_confs_dir_label_figures,
        "rota_outlier_percentage",
    )
    # rama_outlier_percentage
    plot_hist(
        results.conf_stats.rama_outlier_percentage,
        "Ramachandran outlier percentage",
        results.pred_confs_dir_label_figures,
        "rama_outlier_percentage",
    )
    # rama_favored percentage
    plot_hist(
        results.conf_stats.rama_favored_percentage,
        "Ramachandran favored percentage",
        results.pred_confs_dir_label_figures,
        "rama_favored_percentage",
    )
    # clashscore
    plot_hist(
        results.conf_stats.clashscore,
        "Clash Score",
        results.pred_confs_dir_label_figures,
        "clashscore",
    )

    # Make violin-plot/whisker-plots for all ensembles, inc gt
    # Add hlines for recommended value and 6xou value
    # rota_outlier_percentage
    plot_violin(
        ens_geom_lists.labels,
        ens_geom_lists.rota_outlier_percentage,
        "Ensemble",
        "Rotamer Outlier Percentage",
        results.ensemble_stats.rota_outlier_percentage_avg,  # 6xou PDB value
        "6XOU",
        results.pred_confs_dir,  # should be okay as GT,6xou created above
        "rota_outlier_percentage",
    )
    # rama_outlier_percentage
    plot_violin(
        ens_geom_lists.labels,
        ens_geom_lists.rama_outlier_percentage,
        "Ensemble",
        "Ramuchandran Outlier Percentage",
        results.ensemble_stats.rama_outlier_percentage_avg,  # 6xou PDB value
        "6XOU",
        results.pred_confs_dir,  # should be okay as GT,6xou created above
        "rama_outlier_percentage",
    )
    # rama_favored_percentage
    plot_violin(
        ens_geom_lists.labels,
        ens_geom_lists.rama_favored_percentage,
        "Ensemble",
        "Ramuchandran Favored Percentage",
        results.ensemble_stats.rama_outlier_percentage_avg,  # 6xou PDB value
        "6XOU",
        results.pred_confs_dir,  # should be okay as GT,6xou created above
        "rama_favored_percentage",
    )
    # clashscore
    plot_violin(
        ens_geom_lists.labels,
        ens_geom_lists.clashscore,
        "Ensemble",
        "Clash Score",
        results.ensemble_stats.clashscore_avg,  # 6xou PDB value
        "6XOU",
        results.pred_confs_dir,  # should be okay as GT,6xou created above
        "clashscore",
    )

    # Make and save a csv table of all geom stats
    # add 6xou to dict and then grab average vals to csv
    ens_geom_lists.rota_outlier_percentage.append(
        results.conf_stats.rota_outlier_percentage
    )
    ens_geom_lists.rama_outlier_percentage.append(
        results.conf_stats.rama_outlier_percentage
    )
    ens_geom_lists.rama_favored_percentage.append(
        results.conf_stats.rama_favored_percentage
    )
    ens_geom_lists.clashscore.append(results.conf_stats.clashscore)
    ens_geom_lists.labels.append("6XOU")

    # get average values to put in CSV file table
    ens_geom_lists_csv = {
        "labels": ens_geom_lists.labels,
        "rota_outlier_percentage": [
            float(np.average(ens))
            for ens in ens_geom_lists.rota_outlier_percentage
        ],
        "rama_outlier_percentage": [
            float(np.average(ens))
            for ens in ens_geom_lists.rama_outlier_percentage
        ],
        "rama_favored_percentage": [
            float(np.average(ens))
            for ens in ens_geom_lists.rama_favored_percentage
        ],
        "clashscore": [
            float(np.average(ens)) for ens in ens_geom_lists.clashscore
        ],
    }
    # Save csv file table (averaged over ensemble if applicable)
    csv_filename = os.path.join(
        results.pred_confs_dir,
        "all_geom_metrics.csv",
    )
    pd.DataFrame.from_dict(data=ens_geom_lists_csv).to_csv(csv_filename)
    """
    with open(csv_filename, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ens_geom_lists_csv.keys())
        w.writeheader()
        w.writerows()
    """
    # ------------------------------------------------------------------------------------

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    if args.verbose:
        for arg in vars(args):
            print("{}, {}".format(arg, getattr(args, arg)))
    main(args)
