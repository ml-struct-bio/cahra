#!/usr/bin/bash

TRUTH_ENSEMBLE_DIR="metric_calculation/gt_confs/pdbs"
MAX_RANK=100 # updated from 50 to 100

# min-RMSD aligned to GT ensemble PDB locations
lemon1_dir="metric_calculation/pred_confs/lemon1/pdbs"
lemon2_dir="metric_calculation/pred_confs/lemon2/pdbs"
lime1_dir="metric_calculation/pred_confs/lime1/pdbs"
orange1_dir="metric_calculation/pred_confs/orange1/pdbs"
kiwi1_dir="metric_calculation/pred_confs/kiwi1/pdbs"
kiwi2_dir="metric_calculation/pred_confs/kiwi2/pdbs"
kiwi3_dir="metric_calculation/pred_confs/kiwi3/pdbs"
kiwi4_dir="metric_calculation/pred_confs/kiwi4/pdbs"
banana1_dir="metric_calculation/pred_confs/banana1/pdbs"

lemon1_label="lemon1"
lemon2_label="lemon2"
lime1_labeL="lime1"
orange1_label="orange1"
kiwi1_label="kiwi1"
kiwi2_label="kiwi2"
kiwi3_label="kiwi3"
kiwi4_label="kiwi4"
banana1_label="banana1"

# note that banana currently cannot be processed in previous step via MDAnalysis and needs further debugging

cahra-metrics calc_pcv -v -r $MAX_RANK -te $TRUTH_ENSEMBLE_DIR \
-pe $lemon1_dir $lemon2_dir $lime1_dir $orange1_dir $kiwi1_dir $kiwi2_dir $kiwi3_dir $kiwi4_dir $banana1_dir \
-pl $lemon1_label $lemon2_label $lime1_labeL $orange1_label $kiwi1_label $kiwi2_label $kiwi3_label $kiwi4_label $banana1_label \
> calc_pcv4.log 2>&1 &
