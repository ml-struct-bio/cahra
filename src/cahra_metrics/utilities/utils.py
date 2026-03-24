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

class Ensemble_Stats(object):
    def __init__(
        self,
        filenames: list[str] = [],
        ensemble: dict = {},
    ):
        if filenames:
            self.filenames=filenames
        else:
            self.filenames=None
        
        if ensemble:
            self.ensemble = ensemble
        else:
            self.ensemble = None
        # ensemble dictionary has structure sub-dicts accessible via keys
        # "conf_stats" dict - holds per conformation statistics/information
        # "ensemble_stats" dict - holds per ensemble statistics/information
        
        # CONF_STATS
        # "conf_stats" dict keys are the filenames of each conformation
        # each holding a sub-dict for that conformation
        
        # ENSEMBLE_STATS
        # "ensemble_stats" dict keys are...

    
    @classmethod
    def get_filenames_from_dict(self):
        """Filenames are keys in ensemble 
        """
        if self.ensemble["conf_stats"]:
            self.filenames = list(self.ensemble["filenames"].keys())

    def setup_conf_stats_dicts(self):
        if len(self.filenames)>0:
            self.ensemble["conf_stats"] = {}
            for f in self.filenames:
                self.ensemble["conf_stats"][f] = {}
        else:
            raise ValueError("Filenames list empty!")

def natural_sort_key(s, _nsre=re.compile(r'(\d+)')):
    """Sort strings alphanumberically without requiring leading zeros

    Args:
        s _type_: _description_
        _nsre (, optional): Regular expression formulation. Defaults to re.compile(r'(\d+)').

    Returns:
        : _description_
    """
    return [
        int(text) if text.isdigit() else text.lower()for text in _nsre.split(s)
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
    confs = sorted(glob.glob(confs_path), key=natural_sort_key)
    if verbose:
        print(confs)
    return confs