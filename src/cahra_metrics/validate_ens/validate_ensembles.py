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

def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
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
        "--truth_ensemble_dir",
        "-te",
        help=(
            "Path to directory holding mmcif/pdb files which make up the"
            " truth ensemble."
        ),
        type=str,
        default="/media/jg/My Book/JGreer/dinaMISMO/dinaMISMO-systems-prep/PreP/pilot_studies/ConformationSampling/extract_frames_6xou_unsteered_0-200ns/extract_frames_6xou_unsteered_0-200ns",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        help="Increase output verbosity",
        action="store_true",
    )

def main(args):
    # TODO: ADD ensemble validation tools
    # gemmi can be used to convert PDB to MMCIF
    # see gemmi convert: https://gemmi.readthedocs.io/en/latest/program.html
    # TODO: potentially make any alterations required to standardise mmcif file(s) here and resave ensemble(s)
    # can check for hydrogens with gemmi
    # atom.is_hydrogen(): see https://gemmi.readthedocs.io/en/latest/mol.html
    # can remove hydrogens if necessary
    # see structure.remove_hydrogens() from same link
    # TODO check gemmi serial numbers are the same for every atom (ie: element, residue etc etc)
    # for both the GT ensemble and the pred ensemble relative to a single conf of the GT ensemble
    # (which can also be used as the topology file)
    # If gemmi serial numbers are not the same, or the molecule needs hydrogens stripping/similar,
    # recalc the serial numbers after the edit(s) and check they are now identical to the reference.
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    if args.verbose:
        for arg in vars(args):
            print("{}, {}".format(arg, getattr(args, arg)))
    main(args)