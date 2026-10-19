#!/usr/bin/bash

bash run_rota_rama_clash_gt.sh > run_rota_rama_clash_gt.log
wait $!

labels=("poised-bobcat-1" "poised-bobcat-2" "poised-bobcat-3" "exotic-stingray-1" "exotic-stingray-2" "curious-fulmar-1" "curious-fulmar-2" "persuasive-penguins-1" "persuasive-penguins-2" "persuasive-penguins-3" "persuasive-penguins-4" "persuasive-penguins-5" "kiharalab-1" "kiharalab-2")

for label in "${labels[@]}"; do
    echo "Processing label: $label (length: ${#label})"
    bash run_rota_rama_clash_pred.sh $label > run_rota_rama_clash_pred.log
    wait $!
done

bash run_rota_rama_clash_6xou.sh > run_rota_rama_clash_6xou.log
wait $!