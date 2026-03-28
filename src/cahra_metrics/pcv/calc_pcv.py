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
import json
import glob
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
import random
from tabulate import tabulate

from cahra_metrics.utilities.utils import get_ensemble, save_dict_to_json

# ensure reproducibility
# see https://docs.pytorch.org/docs/stable/notes/randomness.html
torch.backends.cudnn.benchmark = False
torch.manual_seed(5)
random.seed(6)
np.random.seed(7)

# TODO update metrics website to note that PCV is calculated using atom coords, not volumes!!!

# TODO update to do PCV on mass/charge weighted aligned GT ensemble (or however it is self-aligned before the pred-ensemble is aligned to it!)

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
        nargs="+", # load multiple ensembles to compare against each other
    )

    parser.add_argument(
        "--pred_ensemble_label",
        "-pl",
        help=(
            "Label for predicted ensemble. Used to set output directories."
            " Use one value per predicted ensemble directory."
        ),
        type=str,
        nargs="+", # one for each pred ensemble
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
        "--max_rank",
        "-r",
        help="Max rank to use in pca computation",
        type=int,
        default=50,
    )

    parser.add_argument(
        "--verbose",
        "-v",
        help="Increase output verbosity",
        action="store_true",
    )

def get_name():
    return "calc_pcv"

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

# TODO: Add function to create and save plot of frames projected into PCA space
# using the right singular vectors of the ground truth PCA computation

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

def load_models_as_tensor(model_filenames: list[str])->torch.Tensor:
    """Grab each conformation from mmcif/pdb files and convert the resulting
     list of mdtraj.Trajectory objects into a torch.Tensor holding the whole
     ensemble

    Args:
        model_filenames (list[str]): list of filenames of conformations

    Returns:
        torch.Tensor: ensemble of atomic conformations, shape=(n_conformations, n_atoms, 3)
    """
    # NOTE: Model_filenames must already be naturally sorted! ie:
    # use utilities.utils.get_ensemble()
    # NOTE: Scaling of 10.0 is used because:
    # 1. MDtraj loads file into mdtraj.Trajectory which uses nanometres as the
    #  unit for coordinates/unit cell lengths
    # 2. To convert to numpy array in Angstroms we multiply coords by 10
    # and return the subsequent conversion to a torch.Tensor
    # NOTE: [0] index is used because each mmcif/pdb has been defined to only
    # include a single conformation
    models = np.asarray([mdtraj.load(f).xyz[0] * 10.0 for f in model_filenames])
    return torch.from_numpy(models)

"""
def prepare_trajectory_for_pca(path_to_pdbs, filename_prefix, topology, verbose=False):
    # TODO COULD REMOVE THIS AND JUST USE THE load_models_as_tensor FUNCTION
    # TO PREPARE TENSORS FOR PCV
    # TODO: IF YOU REMOVE THIS FUNCTION, ADD .xtc WRITING OPTION TO load_models_as_tensor
    # FUNCTION SO YOU CAN VISUALISE THE ENSEMBLES IN MD CODE (such as VMD)

    # TODO sort code to make this conversion once for each (new) pred-ensemble.
    # it uses prepare_trajectory_for_pca() function
    # It should do it once for the truth ensemble but then can refind it each time via
    # a deterministic filename/path

    # TODO HERE is where to check if deterministic .xtc trajectory already exists
    # TODO If so, load that and torchify it instead of making trajectory from mmcif files
    models = get_ensemble(path_to_pdbs)
    
    traj = mdtraj.Trajectory(
        models.numpy() / 10.0, topology=topology
    )
    
    traj.superpose(reference, frame=0)
    traj.save(f"{filename_prefix}_traj.xtc")
    models_as_numpy = torch.from_numpy(traj.xyz * 10.0)

    if verbose:
        print("{} has shape {}".format(os.path.basename(filename_prefix), models_as_numpy.shape))
    return models_as_numpy
"""
    
def main(args):
    # check for CUDA
    if args.verbose:
        print("CUDA available: {}".format(torch.cuda.is_available()))
    # check each pred ensemble has a label
    assert len(
        args.pred_ensemble_dir
    )==len(
        args.pred_ensemble_label
    ), "Supply same number of predicted ensemble directories ({}) and labels ({})".format(
        len(args.pred_ensemble_dir),
        len(args.pred_ensemble_label),
    )


    # STEP 1 - load the pred and truth ensembles
    # NOTE: THESE MUST BE THE SELF-ALIGNED GT AND min-RMSD ALIGNED
    # NOTE: Already saved xtc in calc_min_rmsd_metrics, no need to do it
    # again here
    # NOTE: cahra-metrics calc_min_rmsd_metrics only wrote pdbs instead of
    # mmcifs. Just load the pdbs and use those here
    gt_ensemble_filenames: list[str] = get_ensemble(args.truth_ensemble_dir)
    # We have a list of predicted ensembles
    pred_ensembles_filenames: list[list[str]] = [
        get_ensemble(pred_ensemble) for pred_ensemble in args.pred_ensemble_dir
    ]

    # convert GT mmcifs to tensor
    gt_ensemble_t: torch.Tensor = load_models_as_tensor(gt_ensemble_filenames)
    # convert pred mmcifs of each ensemble to list[torch.Tensor]
    pred_ensembles_t: list[torch.Tensor] = [
        load_models_as_tensor(ensemble) for ensemble in pred_ensembles_filenames
    ]

    if args.verbose:
        print("gt_ensemble shape: {}".format(gt_ensemble_t.shape))
        for (label, ensemble) in zip(args.pred_ensemble_label, pred_ensembles_t):
            print("pred_ensemble {} shape: {}\n".format(label, ensemble.shape))
    
    # NOTE: As long as we are using min-RMSD aligned pred_ensemble mmcif/pdb
    # directories as inputs to args.pred_ensemble_dir we don't need to align
    # against a reference
    # This could otherwise be done with:
    """
    topology = mdtraj.load_pdb(gt_ensemble_filenames[0]).topology
    pred_ensembles_realigned = [
        gt.superpose(
            mdtraj.Trajectory(
                pred_ensemble.numpy() / 10.0,
                topology=topology,
            )
        ) for pred_ensemble in pred_ensembles_t
    ]
    pred_ensembles_t = [
        torch.from_numpy(
            np.asarray(
                pred_ensemble.xyz * 10.0 for pred_ensemble in pred_ensembles_realigned
            )
        )
    ]
    """

    # ------------------------------------------------------------------------------------
    # STEP 3 - compute PCA on trajectories and compute the PCV
    u_gt, s_gt, vh_gt = compute_pca(gt_ensemble_t, max_rank=args.max_rank)
    
    u_pred = []
    s_pred = []
    vh_pred = []
    pcv_pred = []

    for (label, ensemble) in zip(args.pred_ensemble_label, pred_ensembles_t):
        if args.verbose:
            print("Computing captured variance for {}".format(label))
        u, s, vh = compute_pca(ensemble, max_rank=args.max_rank)

        # NOTE: We are calculating PCV relative to GT ensemble (rather than
        # (PCV(ref to GT)+PCV(GT to ref))/2
        # convert to float to store in json
        pcv = float(compute_captured_variance(vh, vh_gt, s_gt))

        u_pred.append(u)
        s_pred.append(s)
        vh_pred.append(vh)
        pcv_pred.append(pcv)

    # ------------------------------------------------------------------------------------
    # STEP 4 - save PCV info to file (use dict format?)
    # TODO: tidy this up to have set filenames in class
    # Also, set pred_confs dir from any of the loaded jsons
    pred_confs_dir = ""
    for (label, pred_ensemble_dir, pcv) in zip(args.pred_ensemble_label, args.pred_ensemble_dir, pcv_pred):
        # print of PCV for each pred ensemble
        print(
            "{}: PCV={}".format(label, pcv)
        )

        if args.verbose:
            print("Adding PCV to json for {}".format(label))
        json_file=os.path.join(
            os.path.dirname(os.path.dirname(pred_ensemble_dir)),
            label,
            "results",
            "rms_metrics.json",
        )
        # in this case results variable is a dict rather than a class
        with open(json_file, "r") as f:
            results=json.load(f)
            results["ensemble"]["ensemble_stats"]["pcv"]=pcv
            # save dict to json with updated name
            updated_json_file = os.path.join(
                os.path.dirname(json_file),
                "all_metrics.json",
            )
            save_dict_to_json(results, updated_json_file)
            # also set dir for collective plots
            if not pred_confs_dir:
                pred_confs_dir=results["pred_confs_dir"]

    # ------------------------------------------------------------------------------------
    # STEP 5 - visualise the projection onto GT and ...
    # plot_s_viz 
    plt.plot(s_gt.numpy(), label="GT")
    for (label, s) in zip(args.pred_ensemble_label, s_pred):
        plt.plot(s.numpy(), label=label)
    plt.yscale("log")
    plt.xlabel("Principal Component")
    plt.legend()
    print("Writing plots to {}".format(pred_confs_dir)) # re-using from last loaded json
    plt.savefig(os.path.join(pred_confs_dir, "s_plot.png"), dpi=300)
    plt.savefig(os.path.join(pred_confs_dir, "s_plot.pdf"))
    plt.clf()

    # plot_s_viz_normalised
    s_gt_cumul = torch.cumsum(s_gt, dim=0) / torch.sum(s_gt)
    plt.plot(s_gt_cumul.numpy(), label="GT")
    for (label, s) in zip(args.pred_ensemble_label, s_pred):
        s_cumul = torch.cumsum(s, dim=0) / torch.sum(s)
        plt.plot(s_cumul.numpy(), label=label)
    plt.yscale("log")
    plt.xlabel("Principal Component")
    plt.legend()
    print("Writing plots to {}".format(pred_confs_dir)) # re-using from last loaded json
    plt.savefig(os.path.join(pred_confs_dir, "s_plot_norm.png"), dpi=300)
    plt.savefig(os.path.join(pred_confs_dir, "s_plot_norm.pdf"))
    plt.clf()

    # plot_projection_onto_GT
    fig, ax = plt.subplots(
        len(args.pred_ensemble_label)+1, # +1 for GT projection
        1,
        figsize=(10, 10*(len(args.pred_ensemble_label)+1)),
        layout="compressed",
        sharex=True,
        sharey=True
    )
    pc0 = 0
    pc1 = 1

    # gt
    coeffs = u_gt @ np.diag(s_gt)
    ax[0].scatter(coeffs[:, pc0], coeffs[:, pc1], c=np.arange(coeffs.shape[0]))
    ax[0].set_title("PCA of GT frames")
    # pred ensembles
    for idx, (label, pred_ensemble, s) in enumerate(
        zip(args.pred_ensemble_label, pred_ensembles_t, s_pred)
    ):
        coeffs_pred_ensemble = project_frames_to_pca(pred_ensemble, vh_gt)
        sc = ax[idx+1].scatter(
            coeffs_pred_ensemble[:, pc0],
            coeffs_pred_ensemble[:, pc1],
            c=np.arange(coeffs_pred_ensemble.shape[0]),
        )
        ax[idx+1].set_title("PCA of {} models".format(label))
        fig.colorbar(sc, ax=ax[idx+1], label="Frame index")
    plt.savefig(os.path.join(pred_confs_dir, "pcv_projection.png"), dpi=300)
    plt.savefig(os.path.join(pred_confs_dir, "pcv_projection.pdf"))
    plt.clf()

    # STEP 6 - for now print out the summary statistics at end of this script
    # TODO: make standalone script to make pd.DataFrame/csv from the jsons
    table = [
        ["label", "sac", "lfc", "pcv", "sac_rank", "lfc_rank", "pcv_rank", "summed_score", "relative_ranking"]
    ]
    lines = []
    for (label, pred_ensemble_dir) in zip(args.pred_ensemble_label, args.pred_ensemble_dir):
        json_file=os.path.join(
            os.path.dirname(os.path.dirname(pred_ensemble_dir)),
            label,
            "results",
            "all_metrics.json",
        )
        with open(json_file, "r") as f:
            results=json.load(f)
            line = [
                label,
                results["ensemble"]["ensemble_stats"]["sac"],
                results["ensemble"]["ensemble_stats"]["lfc"],
                results["ensemble"]["ensemble_stats"]["pcv"],
            ]
            lines.append(line)
    
    # now we have each line, find the rank for each metric
    # NOTE: low is good for SAC and LFC
    # high is good for PCV
    sacs = np.asarray([item[1] for item in lines])
    lfcs = np.asarray([item[2] for item in lines])
    pcvs = np.asarray([item[3] for item in lines])

    sac_ranks = np.argsort(
        np.argsort(sacs)
    )
    lfc_ranks = np.argsort(
        np.argsort(lfcs)
    )
    # higher value is better so flip the index ordering to be negative before 2nd argsort
    pcv_ranks = np.argsort(
        np.argsort(pcvs)*(-1)
    )

    # low is good, high bad 
    summed_scores = sac_ranks+lfc_ranks+pcv_ranks

    # NOTE: low is good, but beware of summed_scores which are identical!
    # Any identical summed scores will not have the same relative_rank
    relative_ranks = np.argsort(
        np.argsort(summed_scores)
    )

    for (
        line,
        sac_rank,
        lfc_rank,
        pcv_rank,
        summed_score,
        relative_rank
    ) in zip(
        lines,
        sac_ranks.tolist(),
        lfc_ranks.tolist(),
        pcv_ranks.tolist(),
        summed_scores.tolist(),
        relative_ranks.tolist()
    ):
        line.append(sac_rank)
        line.append(lfc_rank)
        line.append(pcv_rank)
        line.append(summed_score)
        line.append(relative_rank)

        table.append(line)

    # print the formatted table
    print(
        tabulate(
            table,
            headers="firstrow",
            tablefmt="fancy_grid",
        )
    )
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    if args.verbose:
        for arg in vars(args):
            print("{}, {}".format(arg, getattr(args, arg)))
    main(args)