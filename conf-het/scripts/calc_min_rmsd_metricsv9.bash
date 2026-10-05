#!/usr/bin/bash
METRICS_LABEL="_v9"

LOG_DIR="calc_min_rmsd_metrics${METRICS_LABEL}"
TRUTH_ENSEMBLE_DIR="extract_frames_6xou_unsteered_0-200ns_pdb_reordered"
EXT="pdb"
OUTPUT_DIR="metric_calculation${METRICS_LABEL}"

mkdir $LOG_DIR
mkdir $OUTPUT_DIR

label="poised-bobcat-1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT --write_gt -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="poised-bobcat-2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="poised-bobcat-3"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="exotic-stingray-1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="exotic-stingray-2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="curious-fulmar-1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="curious-fulmar-2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="persuasive-penguins-1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="persuasive-penguins-2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="persuasive-penguins-3"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="persuasive-penguins-4"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="persuasive-penguins-5"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

# MDAnalysis cannot read banana's PDB files. Not sure exactly why, but suspect it 
# was because they were missing ter OXT atom. Solved by manual addition.
label="kiharalab-1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_${label}.log" 2>&1 &
wait $!

label="kiharalab-2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb_reordered"
cahra-metrics calc_min_rmsd_metrics -o $OUTPUT_DIR -pl $label -pe $pred_ensemble_dir -te $TRUTH_ENSEMBLE_DIR  -e $EXT -v > "${LOG_DIR}/calc_rmsd_$label.log" 2>&1 &
wait $!

