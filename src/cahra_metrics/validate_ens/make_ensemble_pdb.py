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
import gemmi
import MDAnalysis as mda
import os
import numpy as np
import matplotlib.pyplot as plt

from cahra_metrics.utilities.utils import get_ensemble, make_dir


def add_arguments(parser: argparse.ArgumentParser):
    """Set up arguments for parsing

    Args:
        parser (argparse.ArgumentParser): Waymarking argument parser

    Returns:
        argparse.ArgumentParser:  Waymarking argument parser with arguments
    """

    parser.add_argument(
        "--ensemble_dir",
        "-d",
        help=(
            "Path to directory holding mmcif/pdb files which make up the"
            " an ensemble."
        ),
        type=str,
        default="",
    )

    parser.add_argument(
        "--ensemble_pdb_dir",
        "-o",
        help=(
            "Path to directory to create ensemble pdb file within. Creates"
            " new directory if doesn't already exist."
        ),
        type=str,
        default="",
    )

    parser.add_argument(
        "--ensemble_pdb_filename",
        "-f",
        help="Ensemble pdb file name",
        type=str,
        default="",
    )

    parser.add_argument(
        "--orig_extension",
        "-e",
        help=("File type (extension)"),
        choices=["cif", "pdb", "json"],
        default="pdb",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        help="Increase output verbosity",
        action="store_true",
    )


def unit_vec(vec: np.ndarray) -> float:
    return vec / np.linalg.norm(vec)


def angle(vec1: np.ndarray, vec2: np.ndarray):
    vec1_unit = unit_vec(vec1)
    vec2_unit = unit_vec(vec2)
    return float(np.arccos(np.clip(np.dot(vec1_unit, vec2_unit), -1.0, 1.0)))


def rad_to_deg(rad):
    return float(180.0 * rad / np.pi)


def calc_cv(u: mda.Universe) -> tuple[list[int], list[float], list[float]]:
    # determine the CV for each MD frame
    # residue ranges
    # D1: 33-288###10-265
    # linker: 562-564###539-541
    # D4: 855-1037###832-1014

    iteration = []
    CVrad = []
    CVangle = []

    for ts in u.trajectory:
        iteration.append(ts.frame)
        # Make atomgroup for each residue range
        d1 = u.select_atoms("resnum 10:265")  # and name CA")
        linker = u.select_atoms("resnum 539:541")  # and name CA")
        d4 = u.select_atoms("resnum 832:1014")  # and name CA")

        # Get COM of each
        d1_com = d1.center_of_mass(compound="group")
        linker_com = linker.center_of_mass(compound="group")
        d4_com = d4.center_of_mass(compound="group")

        # determine angle between COMs using linker in center
        ld1 = d1_com - linker_com
        ld4 = d4_com - linker_com
        radians = angle(ld1, ld4)
        CVrad.append(radians)
        CVangle.append(rad_to_deg(radians))

    return iteration, CVrad, CVangle


def get_name():
    return "make_ensemble_pdb"


def main(args):
    # create output dir if not exists
    make_dir(args.ensemble_pdb_dir, args.verbose)

    # Load ensemble
    ensemble = get_ensemble(
        args.ensemble_dir,
        file_ext=args.orig_extension,
        verbose=args.verbose,
    )
    if args.verbose:
        for conf in ensemble:
            print(os.path.basename(conf))

    # read each ensemble into MDAnalysis Universe
    u = mda.Universe(
        ensemble[0],  # topology
        ensemble,
        in_memory=True,
    )

    if args.verbose:
        backbone = u.select_atoms("backbone")
        print("There are {} backbone atoms".format(len(backbone)))

        calpha = u.select_atoms("name CA")
        print("There are {} CA atoms".format(len(calpha)))

    # compute CV (collective variable)
    iteration, cvrad, cvangle = calc_cv(u)
    if args.verbose:
        for it, ang in zip(iteration, cvangle):
            print("Iter {}: {} deg".format(it, ang))

    # plot CV versus frame index
    # TODO update to save plot in output directory
    plt.plot(iteration, cvangle, marker="o")
    plt.xlabel("Frame index")
    plt.ylabel("Collective variable (deg)")
    plt.savefig(
        os.path.join(
            args.ensemble_pdb_dir,
            "CV_modelindex.png",
        ),
        dpi=300,
    )
    plt.clf()
    # TODO add CV to json dict output

    # order models fom open (high CV angle) to closed
    # find indices needed to reorder CV list
    order_by_CV_idxs = np.array(
        np.flip(np.argsort(np.array(cvangle, dtype=float))),
        dtype=int,
    )
    # save CV ordered from open to closed plot to output dir
    plt.plot(
        np.arange(len(order_by_CV_idxs)),
        np.array(cvangle, dtype=float)[order_by_CV_idxs],
    )
    plt.xlabel("Frame index")
    plt.ylabel("Collective variable (deg)")
    plt.savefig(
        os.path.join(
            args.ensemble_pdb_dir,
            "CV_ordered.png",
        ),
        dpi=300,
    )
    plt.clf()

    if args.verbose:
        print("index order open->closed:{}".format(order_by_CV_idxs))
        print("type:{}".format(type(order_by_CV_idxs)))

    # create ensemble pdb containing all models from most open to most closed
    st = gemmi.Structure()
    for i, (pdbfile, cv, rad) in enumerate(
        zip(
            np.array(ensemble, dtype=str)[order_by_CV_idxs],
            np.array(cvangle, dtype=float)[order_by_CV_idxs],
            np.array(cvrad, dtype=float)[order_by_CV_idxs],
        )
    ):
        if args.verbose:
            print(
                "Model {}/{} is {} with CV {} deg ({} rad)".format(
                    i + 1, len(ensemble), os.path.basename(pdbfile), cv, rad
                )
            )
        # append the structure
        conf_st = gemmi.read_structure(pdbfile)
        conf_st.setup_entities()  # Important for proper entity handling

        model = conf_st[0]
        # TODO ensure model names (from filename) are faithfully retained as name labels
        # model.name=pdbfile
        st.add_model(model)

    # Renumber models sequentially
    st.renumber_models()

    # write ensemble pdb to output dir
    st.write_pdb(
        os.path.join(
            args.ensemble_pdb_dir,
            args.ensemble_pdb_filename,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    if args.verbose:
        for arg in vars(args):
            print("{}, {}".format(arg, getattr(args, arg)))
    main(args)
