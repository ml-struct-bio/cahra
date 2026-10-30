#!/usr/bin/bash

source /opt/xtal/ccp4-9/bin/ccp4.setup-sh
wait $!

# setup env vars
METRICS_LABEL="_label"

LOG_DIR="run_rota_rama_clash${METRICS_LABEL}"
EXT="pdb"
OUTPUT_DIR="metric_calculation${METRICS_LABEL}"

mkdir $LOG_DIR

orig_6xou_dir="${OUTPUT_DIR}/orig_6xou/pdbs"
log_dir="${LOG_DIR}/orig_6xou/${EXT}"
mkdir -p $log_dir
# NOTE: could add hydrogens with reduce
# --- Run Ramachandran analysis on every PDB in the ensemble ---
# ramalyze doesn't use model=<model>
for pdb in ${orig_6xou_dir}/*.${EXT}; do
    base=$(basename "$pdb" .${EXT})
    echo "Running ramalyze on $base"
    molprobity.ramalyze \
        "$pdb" \
        json=True \
        --json-filename "${log_dir}/${base}_rama.json" \
        > "${log_dir}/${base}_rama.log" 2>&1
    wait $!
done
wait $!
# NOTE: could add hydrogens with reduce
# --- Run Rotamer analysis on every PDB in the ensemble ---
# rotalyze does use model=<model>
for pdb in ${orig_6xou_dir}/*.${EXT}; do
    base=$(basename "$pdb" .${EXT})
    echo "Running rotalyze on $base"
    molprobity.rotalyze \
        model="$pdb" \
        json=True \
        --json-filename "${log_dir}/${base}_rota.json" \
        > "${log_dir}/${base}_rota.log" 2>&1
    wait $!
done
wait $!
# NOTE: could add hydrogens with reduce if so, change keep_hydrogens=False to keep_hydrogens=True
# keep_hydrogens=False apparently regenerates them
# --- Run Clashscore analysis on every PDB in the ensemble ---
for pdb in ${orig_6xou_dir}/*.${EXT}; do
    base=$(basename "$pdb" .${EXT})
    echo "Running clashscore on $base"
    molprobity.clashscore \
        model="$pdb" \
        keep_hydrogens=False \
        json=True \
        --json-filename "${log_dir}/${base}_clash.json" \
        > "${log_dir}/${base}_clash.log" 2>&1
    wait $!
done
wait $!