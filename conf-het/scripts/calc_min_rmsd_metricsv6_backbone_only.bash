#!/usr/bin/bash
METRICS_LABEL="_v6_backbone_only"
ATOM_SELECTION="backbone"

LOG_DIR="calc_min_rmsd_metrics${METRICS_LABEL}"
TRUTH_ENSEMBLE_DIR="extract_frames_6xou_unsteered_0-200ns_pdb"
EXT="pdb"
OUTPUT_DIR="metric_calculation${METRICS_LABEL}"

mkdir $LOG_DIR
mkdir $OUTPUT_DIR

label="lemon1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION --write_gt -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="lemon2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="lemon3"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="lime1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="lime2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="orange1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="orange2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &

label="kiwi1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="kiwi2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &

label="kiwi3"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="kiwi4"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="kiwi5"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

# MDAnalysis cannot read banana's PDB files. Not sure exactly why, but suspect it 
# was because they were missing ter OXT atom. Solved by manual addition.
label="banana1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="banana2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -s $ATOM_SELECTION  -v > "${LOG_DIR}/calc_rmsd_$label.log" 2>&1 &
wait $!

