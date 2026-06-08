#!/usr/bin/bash

TRUTH_ENSEMBLE_DIR="extract_frames_6xou_unsteered_0-200ns_pdb"
EXT="pdb"

label="lemon1"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
output_dir="metric_calculation/"
cahra-metrics calc_min_rmsd_metrics -o $output_dir -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT --write_gt -v > calc_rmsd_$label.log 2>&1 &

label="lemon2"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
output_dir="metric_calculation/"
cahra-metrics calc_min_rmsd_metrics -o $output_dir -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT --write_gt -v > calc_rmsd_$label.log 2>&1 &

label="lime1"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
output_dir="metric_calculation/"
cahra-metrics calc_min_rmsd_metrics -o $output_dir -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > calc_rmsd_$label.log 2>&1 &

label="orange1"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
output_dir="metric_calculation/"
cahra-metrics calc_min_rmsd_metrics -o $output_dir -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > calc_rmsd_$label.log 2>&1 &

label="kiwi1"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
output_dir="metric_calculation/"
cahra-metrics calc_min_rmsd_metrics -o $output_dir -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > calc_rmsd_$label.log 2>&1 &

label="kiwi2"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
output_dir="metric_calculation/"
cahra-metrics calc_min_rmsd_metrics -o $output_dir -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > calc_rmsd_$label.log 2>&1 &

label="kiwi3"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
output_dir="metric_calculation/"
cahra-metrics calc_min_rmsd_metrics -o $output_dir -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > calc_rmsd_$label.log 2>&1 &

label="kiwi4"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
output_dir="metric_calculation/"
cahra-metrics calc_min_rmsd_metrics -o $output_dir -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > calc_rmsd_$label.log 2>&1 &

# MDAnalysis cannot read banana's PDB files. Not sure exactly why, but suspect it 
# may be because they are missing ter OXT atom. Solved by manual addition.
label="banana1"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
output_dir="metric_calculation/"
cahra-metrics calc_min_rmsd_metrics -o $output_dir -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > calc_rmsd_$label.log 2>&1 &

