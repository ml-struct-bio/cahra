#!/usr/bin/bash
LOG_DIR="validate_ensembles_logs"
TRUTH_PDB_DIR="submissions/extract_frames_6xou_unsteered_0-200ns"
FILE_EXT="pdb"

mkdir $LOG_DIR

label="poised-bobcat-1"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="poised-bobcat-2"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="poised-bobcat-3"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="exotic-stingray-1"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="exotic-stingray-2"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="curious-fulmar-1"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="curious-fulmar-2"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="persuasive-penguins-1"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="persuasive-penguins-2"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="persuasive-penguins-3"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="persuasive-penguins-4"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="persuasive-penguins-5"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="kiharalab-1"
pred_ensemble_dir="submissions/${label}/confs_pdb"
# I previously removed conf 13 and conf 29 as partial. But this appears to have resulted from
# conert_ensemblesv2.bash not having a "wait $!" after:
# ###"cahra-metrics add_missing_oxt -v -pe $outdir -e $NEW_EXT > add_missing_oxt_kiharalab-1.log 2>&1 &"
# This likely caused a read/write race condition and resulting partial models.
# This has been fixed in convert_ensembles_v3.bash and appears to keep all models complete
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="kiharalab-2"
pred_ensemble_dir="submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!