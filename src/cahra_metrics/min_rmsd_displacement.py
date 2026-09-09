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

import numpy as np
import gemmi

from cahra_metrics.utilities.utils import (
    make_dir,
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
    return bfac.tolist()


def set_b_fac_atomlist(st: gemmi.Structure, bfacs: list[float]):
    for model in st:
        for chain in model:
            for res in chain:
                for bfac, atom in zip(bfacs, res):
                    atom.b_iso = bfac
    return st


def save_b_fac_pdb(
    pdbfile: str,
    st: gemmi.Structure,
    subdir: str,
    bfacs: list[float],
):
    dirpath, pdb_name = os.path.split(pdbfile)
    # check whether bfacs are passed
    # if not, don't update them
    if bfacs:
        bfac_st = set_b_fac_atomlist(
            st,
            bfacs,
        )

    bfac_dir = os.path.join(
        os.path.dirname(dirpath),
        subdir,
    )
    make_dir(bfac_dir)
    bfac_file = os.path.join(
        bfac_dir,
        pdb_name,
    )

    bfac_st.write_pdb(bfac_file)


def get_name():
    return "min_rmsd_displacement"


def main(args):
    # STEP 1 - load the json for each pred ensemble
    if args.verbose:
        print(args.pred_ensemble_json)

    # Create lists for keeping track of min-RMSD atomic displacements
    # One list per pred_ensemble. The sublist will be appended to be
    # len (n_pred_models*n_atoms*mask) where mask is applicable for bb or sc
    min_rmsd_disps: list[float] = []
    backbone_names = {"N", "CA", "C", "O", "OXT"}
    min_rmsd_bb_disps: list[float] = []  # for backbone atoms
    min_rmsd_sc_disps: list[float] = []  # for sidechain atoms

    # Create lists for keeping track of b-factor values for each
    # pred ensemble and for GT (if applicable)
    bfacs: list[list[float]] = []  # one list per pred ens plus gt ens
    bfacs_labels: list[str] = []  # gt and pred labels
    bfacs_diff: list[list[float]] = []  # one list per pred ens
    bfacs_abs_diff: list[list[float]] = []  # one list per pred ens
    bfacs_diff_labels: list[str] = []  # pred labels
    bb_mask: list[int] = (
        []
    )  # fill mask to allow filtering of bfacs by bb or sc

    # Loop over the pred ensembles
    for i, (json, label) in enumerate(
        zip(
            args.pred_ensemble_json,
            args.pred_ensemble_label,
        )
    ):
        # STEP 2 - Open the json file
        if args.verbose:
            print("Loading json: {}".format(json))
            print("Label: {}".format(label))
        with open(json, "r") as f:
            results = json.load(f)

        # TODO For each predicted model violin plot the distances
        # between pred model and min-RMSD GT model atoms
        # For (all, backbone, sidechain) categories

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
            for pred_model, gt_model in zip(pred_st, gt_st):
                for pred_chain, gt_chain in zip(pred_model, gt_model):
                    for pred_res, gt_res in zip(pred_chain, gt_chain):
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
                            atom_distance = round(
                                float(np.absolute(pred_pos - gt_pos)),
                                1,
                            )
                            atom_distances.append(atom_distance)
                            min_rmsd_disps.append(atom_distance)
                            if (
                                pred_atom.name not in backbone_names
                                and not pred_atom.is_hydrogen()
                            ):
                                min_rmsd_sc_disps.append(atom_distance)
                                bb_mask.append(0)
                            else:
                                min_rmsd_bb_disps.append(atom_distance)
                                bb_mask.append(1)
                            # Add the distance for each atom (A)
                            # to the b-factor column for pred model
                            pred_atom.b_iso = atom_distance

            # TODO Add per model min-RMSD model atom distances to json?

            # TODO Save updated json?

            # Save pred_pdb atomic model with b-factor col filled
            # formulate save directory and fname from pred_pdb filename
            save_b_fac_pdb(
                pred_pdb,
                pred_st,
                "atom_dist",
                [],
            )

            # STEP 4 - Get PDB version with
            # Load RMSF for GT models and convert to b-factor
            gt_bfacs = rmsf_to_bfac(
                results["ensemble"]["ensemble_stats"]["gt_rmsf"]
            )

            # Load RMSF for pred models and convert to b-factor
            pred_bfacs = rmsf_to_bfac(
                results["ensemble"]["ensemble_stats"]["pred_rmsf"]
            )

            # Compute difference in b-factor (pred-GT)
            diff_bfacs = np.array(pred_bfacs, dtype=float) - np.array(
                gt_bfacs, dtype=float
            )

            # Compute abs difference in b-factor (pred-GT)
            abs_diff_bfacs = np.absolute(diff_bfacs).tolist()

            # Negative vals can be used for plots but not for pdbs
            diff_bfacs = diff_bfacs.tolist()

            # If first iteration, save an updated PDB file for GT
            # b-factor additions
            if i == 0:
                bfacs.append(gt_bfacs)
                bfacs_labels.append("GT")
                save_b_fac_pdb(
                    results["gt_filenames"][min_rmsd_idx],
                    gt_st,
                    "b_facs",
                    gt_bfacs,
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
                    pred_pdb,
                    gt_st,
                    "b_facs_absdiff",
                    abs_diff_bfacs,
                )

                # Can't save negative b-factors to pdb file?
                # so can't also save b-factor differences to pdbs

                # If first iteration also save the pred model
                # with pred ensembles b-factors
                # This will be closest pred model to gt model
                # but will not always be super close for qual
                # eval, so don't rely on that
                save_b_fac_pdb(
                    pred_pdb,
                    pred_st,
                    "b_facs",
                    pred_bfacs,
                )

        # TODO Save updated PDB files for pred ensemble min-distance
        # and a single PDB file for each b-factor file

    # TODO Save per atom b-factor plots for each pred-ensemble
    # similar to RMSF plots
