#!/usr/bin/bash

bash run_rota_rama_clash_gt_v6.sh > run_rota_rama_clash_gt_v6.log
wait $!

#labels=("lemon1" "lemon2" "lemon3" "lime1" "lime2" "orange1" "orange2" "kiwi1" "kiwi2" "kiwi3" "kiwi4" "kiwi5" "banana1" "banana2")
#
#for label in "${labels[@]}"; do
#    echo "Processing label: $label (length: ${#label})"
#    bash run_rota_rama_clash_pred_v6.sh $label > run_rota_rama_clash_pred_v6.log
#    wait $!
#done