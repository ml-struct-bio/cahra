#!/usr/bin/bash
DIR="validationv2"
TRUTH_PDB_DIR="extract_frames_6xou_unsteered_0-200ns_pdb"
FILE_EXT="pdb"

mkdir $DIR

label="lemon1"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > $DIR/$label.log
wait$!

label="lime1"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > $DIR/$label.log
wait $!

label="orange1"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > $DIR/$label.log
wait $!

label="kiwi1"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > $DIR/$label.log
wait $!

label="banana1"
pred_ensemble_dir="prelim_submissions/$label/confs_pdb"
cahra-metrics validate_ensembles -pe $pred_ensemble_dir -te $TRUTH_PDB_DIR -e $FILE_EXT -v > $DIR/$label.log
wait $!


