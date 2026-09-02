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

# import numpy as np
# import json
from cahra_metrics.utilities.utils import (
    # load_dict_from_json,
    get_stereochem_json_files,
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

    # ------------------------------------------------------------------------------------
    # STEP 4 - make simplified stereochem jsons for gt and orig_6xou calculations
    # Note that the pred ensemble stereochemical calculations are added to
    # the pred ensemble json files

    # ------------------------------------------------------------------------------------
    # STEP 5 - DETERMINE WHAT TO CALCULATE HERE... AVERAGES? ETC?

    # ------------------------------------------------------------------------------------
    # STEP 6 - make plot(s) for rota, rama, clash stats for each ensemble

    # ------------------------------------------------------------------------------------
    # STEP 7 - make visualisation of how different clash, rota, rama stats compare
    # histograms, violin plots etc
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    if args.verbose:
        for arg in vars(args):
            print("{}, {}".format(arg, getattr(args, arg)))
    main(args)
