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
        "--truth_ensemble_dir",
        "-te",
        help=(
            "Path to directory holding mmcif/pdb files which make up the"
            " truth ensemble."
        ),
        type=str,
        default=(
            "/media/jg/My Book/JGreer/dinaMISMO/dinaMISMO-systems-prep/PreP/"
            "pilot_studies/ConformationSampling/extract_frames_6xou_unsteered_0-"
            "200ns/extract_frames_6xou_unsteered_0-200ns",
        ),
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


def ensemble_validation(
    fnames: list[str],
    ref_st: gemmi.Structure,
    ensemble: list[gemmi.Structure],
    verbose: bool = False,
):
    # Can be applied using reference model as first model in this ensemble
    # for self-validation or the first model in the GT ensemble for validation
    # of pred ensemble against GT ensemble
    # GT ensemble must be self-validated before being used to validate
    # pred ensemble.
    # Pred ensemble should also be self-validated
    # Use first ensemble as reference to compare others against

    # returns the pred_conf with (if any are made) the updates

    # check the correct number of models are in the list
    assert len(fnames) == len(
        ensemble
    ), "# fnames ({}) and ensembles ({}) do not match!".format(
        len(fnames),
        len(ensemble),
    )
    # loop through ensemble conformations
    for ensemble_idx, (fname, st) in enumerate(zip(fnames, ensemble)):
        if verbose:
            print("Validate conformation {}: {}".format(ensemble_idx, fname))
        # compare reference and this model from the ensemble
        for model_idx, (ref_model, model) in enumerate(zip(ref_st, st)):
            assert (
                ref_model.num == model.num
            ), "Ensemble {}, model {} numbers do not match: {} vs {}".format(
                ensemble_idx,
                model_idx,
                ref_model.num,
                model.num,
            )
            # check all chains are the same
            for chain_idx, (ref_chain, chain) in enumerate(
                zip(ref_model, model)
            ):
                assert ref_chain.name == chain.name, (
                    "Ensemble {}, model {}, chain {} names do not match: "
                    "{} vs {}"
                ).format(
                    ensemble_idx,
                    model_idx,
                    chain_idx,
                    ref_chain.name,
                    chain.name,
                )
                assert (
                    ref_chain.count_atom_sites() == chain.count_atom_sites()
                ), (
                    "Ensemble {}, model {}, chain {} # atoms in chain {} vs {}"
                    " not the same"
                ).format(
                    ensemble_idx,
                    model_idx,
                    chain_idx,
                    ref_chain.count_atom_sites(),
                    chain.count_atom_sites(),
                )
                # could compare mass, but no need if comparing all atoms are the same
                # check all residues are the same
                for residue_idx, (ref_residue, residue) in enumerate(
                    zip(ref_chain, chain)
                ):
                    assert ref_residue.name == residue.name, (
                        "Ensemble {}, model {}, chain {}, residue {} names do not"
                        " match: {} vs {}"
                    ).format(
                        ensemble_idx,
                        model_idx,
                        chain_idx,
                        residue_idx,
                        ref_residue.name,
                        residue.name,
                    )
                    # check all atoms are the same
                    for atom_idx, (ref_atom, atom) in enumerate(
                        zip(ref_residue, residue)
                    ):
                        # as atoms in a residue don't have a defined order in a cif file
                        # we instead need to check that all atoms in both
                        # are present in the other
                        # Currently not counting # times each is present

                        # check pred_atoms in reference residue
                        assert atom.name in [
                            residue_atom.name for residue_atom in ref_residue
                        ], (
                            "Ensemble {}, model {}, chain {}, residue {}, atom {}"
                            " name {} not present in reference!"
                        ).format(
                            ensemble_idx,
                            model_idx,
                            chain_idx,
                            residue_idx,
                            atom_idx,
                            atom.name,
                        )
                        assert ref_atom.element in [
                            residue_atom.element
                            for residue_atom in ref_residue
                        ], (
                            "Ensemble {}, model {}, chain {}, residue {}, atom {}"
                            " element {} not present in reference!"
                        ).format(
                            ensemble_idx,
                            model_idx,
                            chain_idx,
                            residue_idx,
                            atom_idx,
                            atom.element,
                        )

                        # check ref_atoms in pred residue
                        assert ref_atom.name in [
                            residue_atom.name for residue_atom in residue
                        ], (
                            "Ensemble {}, model {}, chain {}, residue {}, atom {}"
                            " name {} not present in predicted!"
                        ).format(
                            ensemble_idx,
                            model_idx,
                            chain_idx,
                            residue_idx,
                            atom_idx,
                            ref_atom.name,
                        )
                        assert ref_atom.element in [
                            residue_atom.element for residue_atom in residue
                        ], (
                            "Ensemble {}, model {}, chain {}, residue {}, atom {}"
                            " element {} not present in predicted!"
                        ).format(
                            ensemble_idx,
                            model_idx,
                            chain_idx,
                            residue_idx,
                            atom_idx,
                            ref_atom.element,
                        )

                        # check atom serial number is identical

                        assert ref_atom.serial == atom.serial, (
                            "Ensemble {}, model {}, chain {}, residue {},"
                            "atom {} serial numbers do not match: {} vs {}"
                        ).format(
                            ensemble_idx,
                            model_idx,
                            chain_idx,
                            residue_idx,
                            atom_idx,
                            ref_atom.serial,
                            atom.serial,
                        )

                        # check atom is not hydrogen
                        # see https://gemmi.readthedocs.io/en/latest/mol.html
                        assert ~atom.is_hydrogen(), (
                            "Ensemble {}, model {}, chain {}, residue {},"
                            " atom {} is hydrogen"
                        ).format(
                            ensemble_idx,
                            model_idx,
                            chain_idx,
                            residue_idx,
                            atom_idx,
                        )
    return


def ensemble_cross_validation(
    reference: gemmi.Structure, ensemble: list[gemmi.Structure]
):
    pass


def get_name():
    return "validate_ensembles"


def main(args):
    # Load GT ensemble
    gt_ensemble_files = get_ensemble(
        args.truth_ensemble_dir,
        file_ext=args.orig_extension,
        verbose=args.verbose,
    )
    # Load predicted ensemble
    pred_ensemble_files = get_ensemble(
        args.pred_ensemble_dir,
        file_ext=args.orig_extension,
        verbose=args.verbose,
    )

    # TODO: read_structure only rather than have read cif block to also access metadata
    # st = gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)

    # read GT and pred ensembles into lists
    gt_ensemble = []
    for f in gt_ensemble_files:
        gt_ensemble.append(gemmi.read_structure(f))
    pred_ensemble = []
    for f in pred_ensemble_files:
        pred_ensemble.append(gemmi.read_structure(f))
    if args.verbose:
        print(
            "Loaded {} GT models and {} predicted models".format(
                len(gt_ensemble),
                len(pred_ensemble),
            )
        )

    print(
        "GT ensemble self-validation using reference:\n{}\n".format(
            gt_ensemble_files[0]
        )
    )
    # GT self-validation
    ensemble_validation(
        fnames=gt_ensemble_files,
        ref_st=gt_ensemble[0],
        ensemble=gt_ensemble,
        verbose=True,
    )

    print(
        "Pred ensemble self-validation using reference:\n{}\n".format(
            pred_ensemble_files[0]
        )
    )
    # Pred self-validation
    ensemble_validation(
        fnames=pred_ensemble_files,
        ref_st=pred_ensemble[0],
        ensemble=pred_ensemble,
        verbose=True,
    )

    print(
        "Pred ensemble validation against GT using reference:\n{}\n".format(
            gt_ensemble_files[0]
        )
    )
    # Pred validation against GT ref
    ensemble_validation(
        fnames=pred_ensemble_files,
        ref_st=gt_ensemble[0],
        ensemble=pred_ensemble,
        verbose=True,
    )

    # TODO Could also use the supported mmcif validation tool to ensure that
    # all mmcif files.
    # TODO Could also just use gemmi validate command...
    # Would be useful if validate_ensembles modifies mmcif file, but it currently
    # doesn't
    # (ie: can be pip installed via https://pypi.org/project/pdbe-mmcif-validator/
    # and is on github at https://github.com/PDBeurope/mmcif-validator)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    if args.verbose:
        for arg in vars(args):
            print("{}, {}".format(arg, getattr(args, arg)))
    main(args)
