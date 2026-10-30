#!/usr/bin/bash
METRICS_LABEL="_label"

LOG_DIR="make_ensemble_pdb${METRICS_LABEL}"
OUTPUT_DIR="metric_calculation${METRICS_LABEL}"
EXT="pdb"

mkdir $LOG_DIR

# get ensemble pdb for the aligned GT ensemble 
ensemble_dir="${OUTPUT_DIR}/gt_confs/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/gt_confs/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_GT"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="poised-bobcat-1"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="poised-bobcat-2"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="poised-bobcat-3"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="exotic-stingray-1"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="exotic-stingray-2"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="curious-fulmar-1"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="curious-fulmar-2"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="persuasive-penguins-1"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="persuasive-penguins-2"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="persuasive-penguins-3"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="persuasive-penguins-4"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="persuasive-penguins-5"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiharalab-1"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiharalab-2"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

