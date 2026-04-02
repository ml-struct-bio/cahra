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

from cahra_metrics.utilities.utils import get_ensemble


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
            " submitted ensemble."
        ),
        type=str,
        default="",
    )

    parser.add_argument(
        "--orig_extension",
        "-e",
        help=("File type (extension) to convert to "),
        choices=["cif", "pdb", "json"],
        default="cif",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        help="Increase output verbosity",
        action="store_true",
    )


def get_name():
    return "add_missing_oxt"


def add_missing_oxt_atom(pdb: gemmi.Structure) -> gemmi.Structure:
    oxt = gemmi.Atom()
    oxt.name = "OXT"
    oxt.element = gemmi.Element("O")
    oxt.occ = 1.00
    oxt.b_iso = 0.00

    # grab the existing O atom position in terminal residue (for duplication)
    # and assign it as the position of the OXT atom also
    model = pdb[0]
    chain = model[0]
    residue = chain[-1]
    copy_atom = residue.find_atom(
        "O", altloc=" "
    )  # space instead of empty str for altloc
    x = copy_atom.pos.x
    y = copy_atom.pos.y
    z = copy_atom.pos.z

    oxt.pos = gemmi.Position(x, y, z)

    residue.add_atom(oxt).pos
    return pdb


def main(args):

    # Load predicted ensemble
    pred_ensemble_files = get_ensemble(
        args.pred_ensemble_dir,
        file_ext=args.orig_extension,
        verbose=args.verbose,
    )

    pred_ensemble = []
    for f in pred_ensemble_files:
        pred_ensemble.append(gemmi.read_structure(f))
    if args.verbose:
        print(
            "Loaded {} predicted models".format(
                len(pred_ensemble),
            )
        )

    print(
        (
            "Adding missing OXT atom. Duplicating O position in terminal residue"
            " and using name OXT."
            " This will degrade metric score, but probably not in a meaningful way"
            " (~1%? at a guess)."
        )
    )

    # duplicate the terminal O atom to create OXT
    for st in pred_ensemble:
        add_missing_oxt_atom(st)

    # overwrite the files with the added OXT
    for f, st in zip(pred_ensemble_files, pred_ensemble):
        st.write_pdb(f)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    if args.verbose:
        for arg in vars(args):
            print("{}, {}".format(arg, getattr(args, arg)))
    main(args)
