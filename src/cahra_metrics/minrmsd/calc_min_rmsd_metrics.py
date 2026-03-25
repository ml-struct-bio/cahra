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
import os
import re
import copy
import glob2 as glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align, rms
from cahra_metrics.utilities.utils import get_ensemble

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

def align_min_rmsd(
        pred_traj: mda.Universe,
        GT_traj: mda.Universe,
        pred_traj_copy: mda.Universe,
        verbose: bool = False,
) -> tuple[dict, mda.Universe]:

    # loop through each predicted structure and find GT structure with min RMSD
    rmsd_min = {}
    rmsd_min["pred_index"]=[]
    rmsd_min["rmsd_min"]=[]
    rmsd_min["gt_index"]=[]
    # pred_traj_aligned is not yet aligned and is at this point just a copy of pred_traj
    for i, pred in enumerate(zip(pred_traj.trajectory, pred_traj_copy.trajectory)):
        rmsd_vals = []
        for GT in GT_traj.trajectory:
            # calc rmsd
            # align the pred model to the GT model
            # rms.rmsd docs: https://docs.mdanalysis.org/1.1.1/documentation_pages/analysis/rms.html#module-MDAnalysis.analysis.rms
            # it aligns set of coords (arg1) to ref (arg2)
            # superpose implies center=True

            # here we are not excluding hist tag in rmsd calc or in min rmsd alignment
            rmsd_vals.append(
                rms.rmsd(
                    pred_traj.select_atoms("backbone").positions, # TODO decide on whether to use all (heavy) atoms for rmsd calculation
                    GT_traj.select_atoms("backbone").positions, # TODO update as above
                    # ensuring QCP superposition seems to make rms.rmsd and align.alignto
                    # np.close to each other
                    superposition=True,
                )
            )
        rmsd_vals=np.array(rmsd_vals,dtype=float)
        min_index = np.argmin(rmsd_vals)
        # TODO: decide whether to also store pairwise min-rmsds to dictionary. May be useful for later analyses
        # when rmsd_min dict is saved to file?
        rmsd_min["pred_index"].append(i)
        rmsd_min["rmsd_min"].append(rmsd_vals[min_index])
        rmsd_min["gt_index"].append(min_index)
    
    # now that we have the indices of the best rmsd fit for each pred_traj
    # we now align them in the copied pred_traj
    # and fill in the new atomic positions
    for (pred_idx,gt_idx, rmsd_min_val) in zip(
        rmsd_min["pred_index"],
        rmsd_min["gt_index"],
        rmsd_min["rmsd_min"],
    ):
        # change the frame of the pred_traj and pred_traj_copy to correct idx
        pred_traj.trajectory[pred_idx]
        pred_traj_copy.trajectory[pred_idx]
        # change the frame of GT_traj to correct idx 
        GT_traj.trajectory[gt_idx]
        # do the alignment between the two frames using GT as ref
        rmsd_unaligned, rmsd_aligned = align.alignto(
            pred_traj,
            GT_traj,
            select="backbone" # not ignoring hist-tag here
        )
        if verbose:
            print(
                "{}: {} -> {}".format(
                    pred_idx,
                    rmsd_unaligned,
                    rmsd_aligned,
                )
            )
        # check the alignment is approx the same as reported by rms.rmsd above
        assert(
            np.isclose(
            rmsd_min_val,
            rmsd_aligned,
            )
        ), "rmsd_min:{}\nrmsd_aligned:{}".format(
            rmsd_min_val, rmsd_aligned
        )

        # input aligned atomic positions to copy of pred_traj
        # have to do here because pred_traj will reset to unaligned after iteration through loop
        pred_traj_copy.atoms.positions = copy.deepcopy(
            pred_traj.atoms.positions
        )

    return rmsd_min, pred_traj_copy

def calc_lfc(rmsf_gt, rmsf_pred):
    abs_diff_rmsf = np.abs(rmsf_gt-rmsf_pred)
    lfc = np.average(
        abs_diff_rmsf
    )
    lfc_stddev= np.std(
        abs_diff_rmsf
    )
    return lfc, lfc_stddev, abs_diff_rmsf

def make_plot_lfc(ad_rmsf):
    # TODO add docstring with info about non-H atom indices
    plt.plot(range(len(ad_rmsf)), ad_rmsf, color="blue", label="|$\Delta$RMSF|")
    plt.xlabel('Atom index')
    plt.ylabel('RMSF ($\AA$)')
    plt.axvspan(142, 4239, zorder=0, alpha=0.2, color='orange', label='D1')
    plt.axvspan(8565, 8603, zorder=0, alpha=0.2, color='blue', label='LNK')
    plt.axvspan(13183, 16025, zorder=0, alpha=0.2, color='green', label='D4')
    plt.legend()
    plt.plot()
    plt.show()
    plt.clf()
    # TODO convert this to save instead of show plot 

def make_over_under_represent_het_plot(
        rmsf_gt, rmsf_pred
):
    # +ve value of diff means heterogeneity under-estimated
    # -ve value of diff means heterogeneity over-estimated
    diff = rmsf_gt-rmsf_pred
    avg_over = 100 # 100 atoms TODO make configurable
    avgs = []
    for i, atom in enumerate(range(len(diff))):
        if i<avg_over:
            continue
        else:
            avgs.append(np.average(diff[i-avg_over:i]))
        
    avgs = np.array(avgs)
    plt.plot(range(len(diff)), diff, color="red", label="GT-pred")
    plt.plot(range(len(diff))[100:], avgs, color="black", label="100-atom-avg")
    plt.xlabel('Atom number')
    plt.ylabel('RMSF ($\AA$)')
    plt.axvspan(142, 4239, zorder=0, alpha=0.2, color='orange', label='D1')
    plt.axvspan(8565, 8603, zorder=0, alpha=0.2, color='blue', label='LNK')
    plt.axvspan(13183, 16025, zorder=0, alpha=0.2, color='green', label='D4')
    plt.legend()
    plt.plot()
    plt.show()
    plt.clf()
    # TODO convert this to save instead of show plot 


def main(args):
    
    # STEP 1 - load predicted and ground truth ensembles
    # ground truth ensemble
    # TODO check with mmcif format instead of pdb
    gt_ensemble = get_ensemble(
        args.truth_ensemble_dir, file_ext=".pdb", verbose=args.verbose
    )
    if args.verbose:
        for conf in gt_ensemble:
            print(os.path.basename(conf))

    # pred ensemble
    pred_ensemble = get_ensemble(
        args.pred_ensemble_dir, file_ext=".pdb", verbose=args.verbose
    )
    if args.verbose:
        for conf in pred_ensemble:
            print(os.path.basename(conf))
    
    # TODO check all confs in each ensemble have same atoms present and same number of atoms
    # and no hydrogens... ie: add the code to validate all the mmcifs are valid (validate_ensembles)!
    
    # ------------------------------------------------------------------------------------
    # STEP 2 - create mda universes to hold the ensembles
    gt_u = mda.Universe(
        gt_ensemble[0], # topology
        gt_ensemble,
        in_memory=True,
    )

    pred_u = mda.Universe(
        pred_ensemble[0], # topology
        pred_ensemble,
        in_memory=True,
    )

    # make a copy of pred_u which we later fill with the conformations once they are
    # optimally aligned to ground truth ensemble
    pred_u_copy = mda.Universe(
        pred_ensemble[0], # topology
        pred_ensemble,
        in_memory=True,
    )

    # ------------------------------------------------------------------------------------
    # STEP 3 - self-align GT ensemble
    # NOTE: gt universe should already be globally aligned...
    # ie: aligner = align.AlignTraj(u, u, select="backbone",
    # in_memory=True).run() was applied to all 2000 conformations before
    # subsampling for dataset creation.
    # Re-aligning may change the alignment slightly but is done to allow
    # calculation of RMSF metrics against a truth ensemble aligned against
    # itself
    
    # align to reference
    # NOTE: we do NOT remove hist-tag from ground-truth alignment as this
    # may slightly affect the GT RMSF distribution
    # It does not affect RMSDs calculated between GT and pred confs though
    # This should not affect individual min-RMSD alignments of pred
    # conformations to the aligned GT conformations
    # align gt_u
    aligner_gt_u = align.AlignTraj(
        gt_u,
        gt_u,
        select='backbone', # and resnum 10:1014',
    ).run()
    # NOTE: use in_memory so the alignments persist!!!!!
    
    # NOTE: if using mass-weighting in RMSD/RMSF calculations, must also use
    # it to do the GT alignment (as this is via min RMSD!
    # Hence if using mass-weighting keep this alignment AND add weight="mass"
    # as an argument!)
    # NOTE: If using weighting, SAVE GT mass/charge weighted aligned
    # conformations to file (mmcif) for reuse in PCV calculation

    # NOTE: we don't align pred_ensemble to itself.
    # It is only aligned on a per-conf basis to the (now aligned) GT confs

    # ------------------------------------------------------------------------------------
    # STEP 4 - min-RMSD align each conf of pred_ensemble to closest conf in GT ensemble
    # and get new min-RMSD-aligned pred_ensemble
    # dictionary holds min-RMSD value and index of min-RMSD GT conformation for each
    # predicted conformation
    # NOTE: not using mass weighting as RMSF calculation does not support this
    # related... TODO decide how to check hydrogens are not present in mmcif and where to do this.
    pred_u_rmsd_min_dict, pred_u_rmsd_min_aligned = align_min_rmsd(
        pred_u,
        gt_u,
        pred_u_copy,
    )
    if args.verbose:
        for key,val in pred_u_rmsd_min_dict.items():
            print("{}:{}".format(key,val))
    # TODO SAVE mass/charge-weighted min-RMSD-aligned pred-ensemble to file for reuse 
    # in PCV (and any other future comparison to mass/charge-weighted GT-ensemble!)

    # TODO: SAVE self-aligned GT as mmcifs to disk in a new directory
    # (to be used as input for min-RMSD metric recalculation (or
    # PCV metric calculation)???)
    # TODO: SAVE min-RMSD-aligned pred models to disk in a new directory 

    # ------------------------------------------------------------------------------------
    # STEP 5 - get SAC
    SAC = np.average(np.array(pred_u_rmsd_min_dict["rmsd_min"]))
    print(SAC)
    # TODO Save pred_u_rmsd_min_dict and SAC/LFC values to file (use dict for latter?)
    # also could save atom-wise RMSF values for aligned-pred-ensemble (AND GT ensemble,
    # to check consistency of alignment/realignmnet if that is kept in the above code)

    # ------------------------------------------------------------------------------------
    # STEP 6 - calc RMSF on per-residue basis to allow calculation of LFC
    # NOTE: Not mass-weighting this currently. Currently just give all atoms same
    # weighting. The hydrogens should not be present in the submitted (or GT) mmcif files
    # but ensuring this is the case should be done by validate_ensembles code
    # It appears weights can be applied to an atomgroup object, ie provide rms.RMSF with
    # the argument 'weights="mass"', in similar way to RMSD calc above
    # see: https://docs.mdanalysis.org/stable/documentation_pages/core/groups.html#
    # MDAnalysis.core.groups.AtomGroup
    # HOWEVER, rms.RMSF reports that it does not use mass weighting...
    # see https://www.mdanalysis.org/pmda/api/rmsf.html
    # If this is the case, then do not mass/charge weight anything unless you are
    # sure it can be used in all atomic metrics. Do state that RMSDs are NOT MASS
    # WEIGHTED if you do not mass-weight them

    # so here get the RMSF of our aligned GT ensemble for each atom
    rmsf_gt = rms.RMSF(gt_u.select_atoms("protein")).run()
    # and for the min-RMSD-aligned pred ensemble
    rmsf_pred = rms.RMSF(pred_u_rmsd_min_aligned.select_atoms("protein")).run()

    # ------------------------------------------------------------------------------------
    # STEP 7 - get LFC
    # make and save LFC plots
    lfc_all, lfc_stddev, abs_diff_rmsf= calc_lfc(
        rmsf_gt.results.rmsf,
        rmsf_pred.results.rmsf,
    )
    if args.verbose:
        print(
            "EDC:{}\nEDC_stddev:{}".format(
                lfc_all, lfc_stddev
            )
        )

    make_plot_lfc(abs_diff_rmsf)
    make_over_under_represent_het_plot(
        rmsf_gt.results.rmsf,
        rmsf_pred.results.rmsf
    )
    # TODO improve error bar plotting style

    # TODO save min-RMSD aligned predicted ensemble to file
    # TODO option to also visualise/save rmsf via temp factor in (saved) min-RMSD pred-ensemble mda universe
    # see https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsf.html#Visualising-RMSF-as-B-factors
    # NOTE I have done this previously - example in extract_frames_6xou_unsteered_0-200ns_alignment_metric_tests.ipynb
    # (ie: the notebook which showed volume-ensemble-metrics weren't invariant)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    if args.verbose:
        for arg in vars(args):
            print("{}, {}".format(arg, getattr(args, arg)))
    main(args)