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

import math
import argparse
import os

import numpy as np
import gemmi

from cahra_metrics.utilities.utils import (
    # load_dict_from_json,
    get_ensemble,
)


def add_arguments(parser: argparse.ArgumentParser):
    """Set up arguments for parsing

    Args:
        parser (argparse.ArgumentParser): Waymarking argument parser

    Returns:
        argparse.ArgumentParser:  Waymarking argument parser with arguments
    """

    parser.add_argument(
        "--pred_ensemble_dir",
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
        "--verbose",
        "-v",
        help="Increase output verbosity",
        action="store_true",
    )


def make_grid(
    pdblist: list[str], spacing: float, verbose: bool = False
) -> gemmi.FloatGrid:
    """Create a FloatGrid covering the bounding box of the structure."""
    min_pos = np.array([np.inf] * 3)
    max_pos = np.array([-np.inf] * 3)

    for pdb in pdblist:
        st = gemmi.read_structure(pdb)
        for model in st:
            for chain in model:
                for res in chain:
                    for atom in res:
                        pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                        min_pos = np.minimum(min_pos, pos)
                        max_pos = np.maximum(max_pos, pos)

    # set size only using max positions.
    # size = max_pos - min_pos
    size = max_pos
    nx = int(math.ceil(size[0] / spacing)) + 1
    ny = int(math.ceil(size[1] / spacing)) + 1
    nz = int(math.ceil(size[2] / spacing)) + 1

    grid = gemmi.FloatGrid(nx, ny, nz)
    grid.set_unit_cell(gemmi.UnitCell(size[0], size[1], size[2], 90, 90, 90))
    return grid


def ensemble_to_density_map(
    pdblist: list[str],
    grid: gemmi.FloatGrid,
    d_min: float = 3.0,
    rate: float = 1.5,
    blur: float = 0.0,
    verbose=False,
) -> gemmi.FloatGrid:
    """Load a list of PDB files and compute a density map averaged over all
    structures, on a grid sized from the largest extent across all files.

    Args:
        pdblist: List of PDB file paths (one structure per file).
        grid: Grid determined from GT ensemble.
        d_min: Resolution limit for the density calculation.
        rate: Oversampling rate (grid spacing = d_min / (2 * rate)).
        blur: Additional B-factor blur to apply (Angstrom^2).

    Returns:
        gemmi.FloatGrid: Density map averaged over all structures.
    """
    if not pdblist:
        raise ValueError("pdblist is empty")

    # STEP 1 - Read all structures and size the grid from the largest extent
    structures = [gemmi.read_structure(pdb) for pdb in pdblist]

    # STEP 2 - Set up the density calculator on the shared grid
    dencalc = gemmi.DensityCalculatorE()
    dencalc.d_min = d_min
    dencalc.rate = rate
    dencalc.blur = blur
    # Resize the calculator's grid in place (grid property is read-only)
    # Using unit cell derived from GT ensemble
    dencalc.grid.set_size(grid.nu, grid.nv, grid.nw)
    dencalc.grid.set_unit_cell(grid.unit_cell)

    # STEP 3 - Accumulate density from every structure
    acc = np.zeros(
        (dencalc.grid.nu, dencalc.grid.nv, dencalc.grid.nw), dtype=np.float64
    )

    for st, pdb in zip(structures, pdblist):
        if verbose:
            print(f"Adding {pdb} to map gen")
        model = st[0]
        # Zero the grid manually (avoid initialize_grid, which resizes
        # the grid based on d_min/rate and breaks the fixed grid size)
        dencalc.grid.array.fill(0.0)
        dencalc.add_model_density_to_grid(model)
        dencalc.grid.symmetrize_sum()
        acc += np.asarray(dencalc.grid.array, dtype=np.float64)

    # STEP 4 - Average over all structures
    acc /= len(structures)

    # STEP 5 - Build the averaged grid
    result = gemmi.FloatGrid(dencalc.grid.nu, dencalc.grid.nv, dencalc.grid.nw)
    result.set_unit_cell(grid.unit_cell)
    np.copyto(result.array, acc)  # copy into the existing array
    return result


def make_difference_map(
    gt_map: gemmi.FloatGrid,
    pred_map: gemmi.FloatGrid,
    verbose: bool = False,
) -> gemmi.FloatGrid:
    """Compute a difference map (pred - gt) between two density maps.

    Args:
        gt_map: Ground truth density map (FloatGrid).
        pred_map: Predicted density map (FloatGrid).
        verbose: Increase output verbosity.

    Returns:
        gemmi.FloatGrid: Difference map (pred - gt) on the same grid.
    """
    # Check the grids are compatible
    if (gt_map.nu, gt_map.nv, gt_map.nw) != (
        pred_map.nu,
        pred_map.nv,
        pred_map.nw,
    ):
        raise ValueError(
            f"Grid size mismatch: gt={gt_map.nu, gt_map.nv, gt_map.nw} "
            f"pred={pred_map.nu, pred_map.nv, pred_map.nw}"
        )

    # Compute the difference (pred - gt)
    diff = np.asarray(pred_map.array, dtype=np.float64) - np.asarray(
        gt_map.array, dtype=np.float64
    )

    # Build the difference grid on the same lattice as gt_map
    result = gemmi.FloatGrid(gt_map.nu, gt_map.nv, gt_map.nw)
    result.set_unit_cell(gt_map.unit_cell)
    np.copyto(result.array, diff)

    if verbose:
        print(
            f"Difference map: {result.nu}x{result.nv}x{result.nw}, "
            f"range [{diff.min():.3f}, {diff.max():.3f}]"
        )

    return result


def save_gemmi_grid(ensemble_density, output_path="ensemble_map.mrc"):
    """Save a gemmi FloatGrid as a CCP4/MRC map file."""
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = ensemble_density
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map(output_path)


def make_dir(dirname: str):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)


def get_name():
    return "ensemble_maps"


def main(args):
    # STEP 1 - Get GT and Pred ensemble PDB file lists
    if args.verbose:
        print("GT pdb dir: {}".format(args.truth_ensemble_dir))
    gt_ensemble_filenames: list[str] = get_ensemble(args.truth_ensemble_dir)
    if args.verbose:
        print("GT models: {}".format(gt_ensemble_filenames))

    # We have a list of predicted ensembles
    if args.verbose:
        for pred_dir in args.pred_ensemble_dir:
            print("Pred pdb dir: {}".format(pred_dir))
    pred_ensembles_filenames: list[list[str]] = [
        get_ensemble(pred_ensemble) for pred_ensemble in args.pred_ensemble_dir
    ]
    if args.verbose:
        for pred_ens in pred_ensembles_filenames:
            print("Pred models: {}".format(pred_ens))

    # Make the grid from GT ensemble
    grid = make_grid(gt_ensemble_filenames, spacing=1.0)
    print(grid)

    # STEP 2 - Get density map for GT model
    gt_map = ensemble_to_density_map(
        gt_ensemble_filenames,
        grid,
        verbose=args.verbose,
    )

    # STEP 3 - Get density map for each Pred model
    pred_maps = []
    for pred_ens in pred_ensembles_filenames:
        pred_maps.append(
            ensemble_to_density_map(
                pred_ens,
                grid,
                verbose=args.verbose,
            )
        )
    # STEP 4 - Get difference map for each Pred model
    diff_maps = []
    for pred_map in pred_maps:
        diff_maps.append(
            make_difference_map(gt_map, pred_map, verbose=args.verbose)
        )

    # STEP 5 - Create directories for ensemble maps
    # and save ensemble and difference maps to file
    for (
        pred_ensemble_dir,
        label,
        pred_map,
        diff_map,
    ) in zip(
        args.pred_ensemble_dir,
        args.pred_ensemble_label,
        pred_maps,
        diff_maps,
    ):
        # create map directories for pred ensembles
        map_dir = os.path.join(
            os.path.dirname(pred_ensemble_dir),
            "ensemble_map",
        )
        make_dir(map_dir)
        save_gemmi_grid(
            pred_map,
            os.path.join(
                map_dir,
                "{}_map.mrc".format(label),
            ),
        )
        # for pred diffmaps
        map_dir = os.path.join(
            os.path.dirname(pred_ensemble_dir),
            "ensemble_diffmap",
        )
        make_dir(map_dir)
        save_gemmi_grid(
            diff_map,
            os.path.join(
                map_dir,
                "{}_diffmap.mrc".format(label),
            ),
        )
    # for gt map
    map_dir = os.path.join(
        os.path.dirname(args.truth_ensemble_dir),
        "ensemble_map",
    )
    make_dir(map_dir)
    save_gemmi_grid(
        gt_map,
        os.path.join(
            map_dir,
            "gtmap.mrc",
        ),
    )
