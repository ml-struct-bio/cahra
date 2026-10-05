#!/usr/bin/bash
METRICS_LABEL="_v6"

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
label="lemon1"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="lemon2"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="lemon3"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="lime1"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="lime2"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="orange1"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="orange2"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiwi1"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiwi2"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiwi3"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiwi4"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiwi5"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="banana1"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="banana2"
ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/pdbs"
output_ensemble_dir="${OUTPUT_DIR}/pred_confs/${label}/ensemble_pdb"
ensemble_pdb_label="ensemble_pdb_${label}"
fname=$ensemble_pdb_label".pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v > "${LOG_DIR}/${ensemble_pdb_label}.log"
wait $!

