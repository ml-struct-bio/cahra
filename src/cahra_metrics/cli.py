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

import os
import argparse
import pkg_resources
from pathlib import Path

import cahra_metrics.minrmsd.calc_min_rmsd_metrics
import cahra_metrics.pcv.calc_pcv

def get_cahra_metrics_parent() -> Path:
    return Path(__file__).parent.parent.parent

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--version",
        action="version",
        version="CAHRA-metrics version {}, license is found in {}".format(
            pkg_resources.get_distribution("cahra-metrics").version,
            os.path.join(get_cahra_metrics_parent(), "LICENSE"),
        ),
    )

    subparsers = parser.add_subparsers(
        title="subcommands",
        description="valid subcommands",
        help="",
    )
    subparsers.required = True

    modules = [
        cahra_metrics.minrmsd.calc_min_rmsd_metrics,
        cahra_metrics.pcv.calc_pcv,
    ]

    module_helptext = [
        "Sampling a molecular dynamics trajectory and saving the"
        + " conformations to PDB files.",
        "Simulation of micrograph/tomogram dataset using Parakeet software.",
    ]

    for helptext, module in zip(module_helptext, modules):
        this_parser = subparsers.add_parser(
            module.get_name(),
            help=helptext,
        )
        module.add_arguments(this_parser)
        this_parser.set_defaults(func=module.main)

    args = parser.parse_args()

    args.func(args)


if __name__ == "__main__":
    main()
