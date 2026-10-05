#!/usr/bin/bash

METRICS_LABEL="_v6"
LOG_DIR="ensemble_maps${METRICS_LABEL}"

mkdir -p $LOG_DIR

lemon1label="lemon1"
pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${lemon1label}/pdbs"

gtlabel="GT"
gt_ensemble_dir="metric_calculation${METRICS_LABEL}/gt_confs/pdbs"

cahra-metrics ensemble_maps \
-pe $pred_ensemble_dir \
-pl $lemon1label \
-te $gt_ensemble_dir \
-v > "${LOG_DIR}/ensemble_maps.log"

