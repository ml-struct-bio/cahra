#!/usr/bin/bash

TRUTH_ENSEMBLE_DIR="metric_calculation/gt_confs/pdbs"
MAX_RANK=50

# min-RMSD aligned to GT ensemble PDB locations
lemon1_dir="metric_calculation/pred_confs/lemon1/pdbs"
lime1_dir="metric_calculation/pred_confs/lime1/pdbs"
orange1_dir="metric_calculation/pred_confs/orange1/pdbs"
kiwi1_dir="metric_calculation/pred_confs/kiwi1/pdbs"
banana1_dir="metric_calculation/pred_confs/banana1/pdbs"

lemon1_label="lemon1"
lime1_labeL="lime1"
orange1_label="orange1"
kiwi1_label="kiwi1"
banana1_label="banana1"

# note that banana currently cannot be processed in previous step via MDAnalysis and needs further debugging

cahra-metrics calc_pcv -v -r $MAX_RANK -te $TRUTH_ENSEMBLE_DIR \
-pe $lemon1_dir $lime1_dir $orange1_dir $kiwi1_dir $banana1_dir \
-pl $lemon1_label $lime1_labeL $orange1_label $kiwi1_label $banana1_label \
> calc_pcv2.log 2>&1 &
