#!/usr/bin/bash
LOG_DIR="validate_ensembles_logs_v6"
TRUTH_PDB_DIR="extract_frames_6xou_unsteered_0-200ns_pdb"
FILE_EXT="pdb"

mkdir $LOG_DIR

label="lemon1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="lemon2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="lemon3"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="lime1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="lime2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="orange1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="orange2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="kiwi1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="kiwi2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="kiwi3"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="kiwi4"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="kiwi5"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="banana1"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
# I previously removed conf 13 and conf 29 as partial. But this appears to have resulted from
# conert_ensemblesv2.bash not having a "wait $!" after:
# ###"cahra-metrics add_missing_oxt -v -pe $outdir -e $NEW_EXT > add_missing_oxt_banana1.log 2>&1 &"
# This likely caused a read/write race condition and resulting partial models.
# This has been fixed in convert_ensembles_v3.bash and appears to keep all models complete
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!

label="banana2"
pred_ensemble_dir="prelim_submissions/${label}/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > "${LOG_DIR}/${label}.log" 2>&1 &
wait $!