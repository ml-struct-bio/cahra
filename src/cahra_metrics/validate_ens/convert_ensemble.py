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
import gemmi
import os
import subprocess

from cahra_metrics.utilities.utils import get_ensemble, replace_file_ext, make_dir

def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Set up arguments for parsing

    Args:
        parser (argparse.ArgumentParser): Waymarking argument parser

    Returns:
        argparse.ArgumentParser:  Waymarking argument parser with arguments
    """

    parser.add_argument(
        "--input_ensemble_dir",
        "-i",
        help=(
            "Path to directory holding files which make up the"
            " ensemble to convert the format of."
        ),
        type=str,
        default="",
    )

    parser.add_argument(
        "--output_ensemble_dir",
        "-o",
        help=(
            "Path to directory to hold files which make up the"
            " converted ensemble."
        ),
        type=str,
        default="",
    )

    parser.add_argument(
        "--update-chainid",
        "-c",
        help=(
            "Update the chain ID, in gemmi format, ie: X:A"
        ),
        type=str,
        default="A:A",
    )

    parser.add_argument(
        "--orig-extension",
        "-oe",
        help=(
            "File type (extension) to convert to "
        ),
        choices=["cif", "pdb", "json"],
        default="pdb"
    )

    parser.add_argument(
        "--extension",
        "-e",
        help=(
            "File type (extension) to convert to "
        ),
        choices=["cif", "pdb", "json"],
        default="cif"
    )

    parser.add_argument(
        "--verbose",
        "-v",
        help="Increase output verbosity",
        action="store_true",
    )

def get_name():
    return "convert_ensemble"

# gemmi can be used to convert between PDB, MMCIF and MMJSON
# see gemmi convert: https://gemmi.readthedocs.io/en/latest/program.html
def main(args):
    # get ordered list of ensemble conformations from directory
    original_ens = get_ensemble(
        args.input_ensemble_dir,
        file_ext=args.orig_extension,
        verbose=args.verbose,
    )

    make_dir(args.output_ensemble_dir, verbose=args.verbose)
    
    # determine updated file name
    for f in original_ens:
        new_ensemble_filename=replace_file_ext(
            f,
            args.extension,
            args.output_ensemble_dir,
            args.verbose,
        )
        
        # change chain name or fill empty chain before doing rest of conversion
        subprocess.run(
            [
                "gemmi",
                "convert",
                "--rename-chain={}".format(args.update_chainid),
                "-v",
                f,
                f,
            ]
        )

        # Use gemmi cmd line to convert and save to new location
        subprocess.run(
            [
                "gemmi",
                "convert",
                "--verbose",
                "--remove-h", # remove hydrogens if present
                "--remove-waters", # remove water if present
                "--rename-chain={}".format(args.update_chainid),
                "-v",
                f,
                new_ensemble_filename,
            ],
            capture_output=args.verbose,
        )
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    if args.verbose:
        for arg in vars(args):
            print("{}, {}".format(arg, getattr(args, arg)))
    main(args)