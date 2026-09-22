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
from importlib.metadata import version
from pathlib import Path

import cahra_metrics.validate_ens.convert_ensemble
import cahra_metrics.validate_ens.add_missing_oxt
import cahra_metrics.validate_ens.validate_ensembles
import cahra_metrics.validate_ens.make_ensemble_pdb
import cahra_metrics.minrmsd.calc_min_rmsd_metrics
import cahra_metrics.pcv.calc_pcv
import cahra_metrics.stereochemistry.stereochem_analysis
import cahra_metrics.analyse_aligned.ensemble_maps
import cahra_metrics.analyse_aligned.min_rmsd_displacement
import cahra_metrics.analyse_aligned.plot_rmsf


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
            version("cahra-metrics"),
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
        cahra_metrics.validate_ens.convert_ensemble,
        cahra_metrics.validate_ens.add_missing_oxt,
        cahra_metrics.validate_ens.validate_ensembles,
        cahra_metrics.validate_ens.make_ensemble_pdb,
        cahra_metrics.minrmsd.calc_min_rmsd_metrics,
        cahra_metrics.pcv.calc_pcv,
        cahra_metrics.stereochemistry.stereochem_analysis,
        cahra_metrics.analyse_aligned.ensemble_maps,
        cahra_metrics.analyse_aligned.min_rmsd_displacement,
        cahra_metrics.analyse_aligned.plot_rmsf,
    ]

    module_helptext = [
        (
            "Convert an ensemble of conformations from one file-type to another via"
            "gemmi"
        ),
        "Add missing OXT atom (duplicating position of terminal residue O atom)",
        "Validate an ensemble against itself and/or a reference ensemble.",
        "Create ensemble pdb files with all conformations and labelled models",
        (
            "All to all optimal alignment of ensemble(s) to reference ensemble via"
            " min-RMSD heavy-atom alignment and calculation of resulting RMSD"
            " and RMSF-based metrics."
        ),
        (
            "Calculate proportion of captured variance between ensemble(s)"
            " and reference ensemble. Ensure ensembles are pre-aligned (which"
            " can be done by running RMSD/RMSF metric calculation)."
        ),
        (
            "Take in molprobity rota, rama and clashscore jsons and add summary info"
            " for each model to the metrics jsons. Creates dummy metric jsons for GT"
            " and 6xou so that plotting of stereochemistry is uniform."
            " Includes making the stereochemistry plots."
        ),
        "Create density maps and difference maps from aligned ensembles.",
        (
            "Determine the distances between atoms in the min-RMSD-aligned "
            " pred ensemble"
            " and their min-RMSD GT model. Make violinplots for overall,"
            " backbone only and"
            " sidechain only. Also make violinplots for each pred ensemble"
            " (including GT)"
            " showing the b-factor values for each atom, as well as for"
            " backbone only and"
            "sidechain only. Outputs a pdb file for each model in each"
            " predicted ensemble"
            " where the b-factor column is filled by the distance to the"
            " same atom in the"
            " min-RMSD GT model. Also outputs 1 pdb file for the GT model"
            " where b-factor"
            " column is filled with the b-factor of each atom in the ensemble."
            " Also outputs 1 pdb file for each pred ensemble where the"
            " b-factor column is"
            " filled with the b-factor of each atom in the predicted"
            " ensemble."
            " Also outputs 1 pdb file for each pred ensemble where"
            " the GT model with"
            " the min-RMSD to the pred model (ie: the one saved with"
            " b-factors)"
            " is saved with the b-factor column filled by the absolute"
            " difference in b-factor."
            " These pdbs may be useful for visualisation."
        ),
        "Plot the rmsf",
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
