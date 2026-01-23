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
import torch
import mdtraj
import os
import glob
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import random

# ensure reproducibility
# see https://docs.pytorch.org/docs/stable/notes/randomness.html
torch.backends.cudnn.benchmark = False
torch.manual_seed(5)
random.seed(6)
np.random.seed(7)

# TODO update metrics website to note that PCV is calculated using atom coords, not volumes!!!


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

def get_name():
    return "calc_min_rmsd_metrics"

def compute_pca(atomic_models, max_rank=50):
    """
    Given a set of atomic structures, compute the principal components using SVD.
    
    ** Arguments **
    atomic_models : torch.Tensor
        Tensor of shape (num_frames, num_atoms, 3)
    max_rank : int
        Maximum number of principal components to compute

    ** Returns **
    u : torch.Tensor
        Left singular vectors
    s : torch.Tensor
        Singular values
    vh : torch.Tensor
        Right singular vectors

    - u @ torch.diag(s) gives the projection of the atomic_models onto the principal components.
    - vh contains the eigenvectors of the covariance matrix, and serves as a projector 
        from the original space to the PCA space.
    """
    atomic_models_pca = atomic_models.reshape(atomic_models.shape[0], -1)
    atomic_models_pca -= atomic_models_pca.mean(axis=0, keepdim=True)
    u, s, vh = torch.svd_lowrank(atomic_models_pca, q=max_rank)
    # NOTE: In general, use the full-rank SVD implementation torch.linalg.svd()
    # for dense matrices due to its 10x higher performance characteristics.
    # The low-rank SVD will be useful for huge sparse matrices that
    # torch.linalg.svd() cannot handle.
    return u, s, vh

def project_frames_to_pca(atomic_models, v):
    """
    Given a set of atomic structures and the right singular vectors from PCA,
    project the structures onto the PCA space.

    ** Arguments **
    atomic_models : torch.Tensor
        Tensor of shape (num_frames, num_atoms, 3)
    v : torch.Tensor
        Right singular vectors from PCA
    """
    atomic_models_pca = atomic_models.reshape(atomic_models.shape[0], -1)
    atomic_models_pca -= atomic_models_pca.mean(axis=0, keepdim=True)
    coeffs = atomic_models_pca @ v
    return coeffs

def compute_captured_variance(
    eigvecs_test_subspace,
    eigvecs_ref_subspace,
    singvals_ref_subspace,
):
    """
    Compute the captured variance of the test subspace with respect to the reference subspace.
    
    ** Arguments **
    eigvecs_test_subspace : torch.Tensor
        Eigenvectors of the test subspace
    eigvecs_ref_subspace : torch.Tensor
        Eigenvectors of the reference subspace
    singvals_ref_subspace : torch.Tensor
        Singular values of the reference subspace
    """
    US = eigvecs_ref_subspace @ torch.diag(singvals_ref_subspace)
    return torch.norm(torch.adjoint(eigvecs_test_subspace) @ US) ** 2 / torch.sum(
        singvals_ref_subspace**2
    )

def load_pdb_models(path_to_models):
    # TODO convert to mmcif compatible
    pdb_files = natsorted(glob.glob(os.path.join(path_to_models, "*.pdb")))
    models = np.asarray([mdtraj.load(f).xyz[0] * 10.0 for f in pdb_files])
    return torch.from_numpy(models)

# TODO convert the argparse input
reference = mdtraj.load_pdb("0_6xou_0-200ns_115.pdb")
topology = reference.topology

def prepare_trajectory_for_pca(path_to_pdbs, filename_prefix, topology, verbose=False):
    # TODO sort code to make this conversion once for each (new) pred-ensemble.
    # it uses prepare_trajectory_for_pca() function
    # It should do it once for the truth ensemble but then can refind it each time via
    # a deterministic filename/path

    # TODO HERE is where to check if deterministic .xtc trajectory already exists
    # TODO If so, load that and torchify it instead of making trajectory from mmcif files
    models = load_pdb_models(path_to_pdbs)
    
    traj = mdtraj.Trajectory(
        models.numpy() / 10.0, topology=topology
    )
    
    traj.superpose(reference, frame=0)
    traj.save(f"{filename_prefix}_traj.xtc")
    models_as_numpy = torch.from_numpy(traj.xyz * 10.0)

    if verbose:
        print("{} has shape {}".format(os.path.basename(filename_prefix), models_as_numpy.shape))
    return models_as_numpy


def main(args):

    # STEP 1 - convert pred and truth ensembles to .xtc format
    # AND load xtc trajectories from file and convert to torch
    # TODO use prepare_trajectory_for_pca(), see above

    # ------------------------------------------------------------------------------------
    # STEP 2 - compute PCA on trajectories

    # ------------------------------------------------------------------------------------
    # STEP 3 - compute the PCV
    # TODO determine whether average is correct way to get PCV (PCV(ref to GT)+PCV(GT to ref))/2  

    # ------------------------------------------------------------------------------------
    # STEP 4 - save PCV info to file (use dict format?)

    # ------------------------------------------------------------------------------------
    # STEP 5 - visualise the # TODO what and why?

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    if args.verbose:
        for arg in vars(args):
            print("{}, {}".format(arg, getattr(args, arg)))
    main(args)