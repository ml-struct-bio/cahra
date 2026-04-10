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

# TODO: add docstrings

import os
import re
import glob
import json
from natsort import natsorted


class Ensemble_Stats(object):
    def __init__(
        self,
        output_dir: str,
        pred_ensemble_label: str,
        gt_filenames: list[str] = [],
        pred_filenames: list[str] = [],
        ensemble: dict = {},
        verbose: bool = False,
    ):
        self.output_dir = output_dir
        self.pred_ensemble_label = pred_ensemble_label
        if gt_filenames:
            self.gt_filenames = gt_filenames
        self.pred_filenames: list[str] = pred_filenames

        # ensemble dictionary has sub-dicts accessible via keys:

        # CONF_STATS
        # "conf_stats" dict - holds per conformation statistics/information
        # keys: pred_index: list, pred_filename: list, rmsd_min: list, gt_index: list]

        # ENSEMBLE_STATS
        # "ensemble_stats" dict - holds per ensemble statistics/information
        # keys:
        # pred_ensemble_label: str,
        # sac: float,
        # lfc: float,
        # lfc_abs_diff: list[float],
        # pcv: float,
        # gt_rmsf: list[float],
        # pred_rmsf: list[float]
        # In "rmsf" we have per-atom rmsf values for the gt_ensemble and
        # pred_ensemble (in Angstroms)

        self.verbose = verbose

        self.gt_confs_dir = os.path.join(self.output_dir, "gt_confs")
        self.gt_confs_dir_cifs = os.path.join(self.gt_confs_dir, "cifs")
        self.gt_confs_dir_pdbs = os.path.join(self.gt_confs_dir, "pdbs")
        self.gt_confs_dir_xtc = os.path.join(self.gt_confs_dir, "xtc")

        self.pred_confs_dir = os.path.join(self.output_dir, "pred_confs")
        self.pred_confs_dir_label = os.path.join(
            self.pred_confs_dir, self.pred_ensemble_label
        )
        self.pred_confs_dir_label_cifs = os.path.join(
            self.pred_confs_dir_label, "cifs"
        )
        self.pred_confs_dir_label_pdbs = os.path.join(
            self.pred_confs_dir_label, "pdbs"
        )
        self.pred_confs_dir_label_xtc = os.path.join(
            self.pred_confs_dir_label, "xtc"
        )
        self.pred_confs_dir_label_figures = os.path.join(
            self.pred_confs_dir_label, "figures"
        )
        self.pred_confs_dir_label_results = os.path.join(
            self.pred_confs_dir_label, "results"
        )

        print(self.__dict__)

    def make_results_dir_structure(self):
        # top level
        make_dir(self.output_dir, verbose=self.verbose)

        # gt ensemble self-aligned
        make_dir(self.gt_confs_dir)
        make_dir(self.gt_confs_dir_cifs)
        make_dir(self.gt_confs_dir_pdbs)
        make_dir(self.gt_confs_dir_xtc)

        # pred ensemble min-RMSD-aligned
        make_dir(self.pred_confs_dir)
        make_dir(self.pred_confs_dir_label)
        make_dir(self.pred_confs_dir_label_cifs)
        make_dir(self.pred_confs_dir_label_pdbs)
        make_dir(self.pred_confs_dir_label_xtc)
        make_dir(self.pred_confs_dir_label_figures)
        make_dir(self.pred_confs_dir_label_results)


def natural_sort_key(s, _nsre=re.compile(r"(\\d+)")):
    """Sort strings alphanumberically without requiring leading zeros
        Now replaced by natsorted # TODO: remove?
    Args:
        s _type_: _description_
        _nsre (, optional): Regular expression formulation.
        Defaults to re.compile(r'(\\d+)').

    Returns:
        : _description_
    """
    return [
        int(text) if text.isdigit() else text.lower()
        for text in _nsre.split(s)
    ]


def get_ensemble(
    conformations_dir: str,
    file_ext: str = ".pdb",
    verbose: bool = False,
) -> list[str]:
    """_summary_

    Args:
        conformations_dir (str): _description_
        file_ext (str, optional): _description_. Defaults to ".pdb".
        verbose (bool, optional): _description_. Defaults to False.

    Returns:
        list[str]: _description_
    """
    confs_path = os.path.join(conformations_dir, "*{}".format(file_ext))
    confs = natsorted(glob.glob(confs_path))
    if verbose:
        print(confs)
    return confs


def replace_file_ext(
    fname: str, new_ext: str, new_dir: str, verbose: bool = False
) -> str:
    # find original extension
    dir_path, base_name = os.path.split(fname)
    orig_basename, orig_ext = os.path.splitext(base_name)
    # check new extension is different to original extension
    assert (
        orig_ext != new_ext
    ), "New file extension ({}) is same as original ({})!".format(
        new_ext, orig_ext
    )
    # replace extension
    new_fname = os.path.join(
        new_dir,
        ".".join([orig_basename, new_ext]),
    )
    if verbose:
        print("Will convert {} to {}".format(fname, new_fname))
    return new_fname


def make_dir(directory: str, verbose: bool = False):
    if not os.path.isdir(directory):
        os.makedirs(directory)
        if verbose:
            print("Making {} results directory".format(directory))


def save_dict_to_json(dictionary: dict, json_file: str, indent: int = 4):
    with open(json_file, "w") as f:
        json.dump(dictionary, f, indent=indent)


# def load_dict_from_json(json_file: str):
#     return json.load(json_file)
