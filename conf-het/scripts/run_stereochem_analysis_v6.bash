#!/usr/bin/bash

METRICS_LABEL="_v6"
LOG_DIR="stereochem_analysis${METRICS_LABEL}"
RESULTS_DIR="metric_calculation${METRICS_LABEL}"

mkdir $LOG_DIR

STEREOCHEM_DIR="run_rota_rama_clash${METRICS_LABEL}"

GT_SUBDIR="gt_confs"
ORIG_SUBDIR="orig_6xou"
PRED_SUBDIR="pred_confs"

lemon1_label="lemon1"
lemon2_label="lemon2"
lemon3_label="lemon3"
lime1_label="lime1"
lime2_label="lime2"
orange1_label="orange1"
orange2_label="orange2"
kiwi1_label="kiwi1"
kiwi2_label="kiwi2"
kiwi3_label="kiwi3"
kiwi4_label="kiwi4"
kiwi5_label="kiwi5"
banana1_label="banana1"
banana2_label="banana2"

cahra-metrics stereochem_analysis -v \
-o $RESULTS_DIR -s $STEREOCHEM_DIR -gts $GT_SUBDIR -os $ORIG_SUBDIR -ps $PRED_SUBDIR \
-pl $lemon1_label $lemon2_label $lemon3_label $lime1_label $lime2_label $orange1_label $orange2_label $kiwi1_label $kiwi2_label $kiwi3_label $kiwi4_label $kiwi5_label $banana1_label $banana2_label \
> $LOG_DIR/stereochem_analysis${METRICS_LABEL}.log 2>&1 &
