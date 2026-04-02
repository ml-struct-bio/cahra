# CAHRA-metrics
Automated metrics calculation for 2025 CAHRA Challenge conf-het dataset

## Running the conf-het metrics

1. Download the conformations
2. Put entries into directory `prelim_scoring/prelim_submissions/<entrylabel>/confs`
3. Copy scripts to `prelim_scoring` and run as below:
```
cp scripts/*.bash prelim_scoring
pushd prelim_scoring
convert_ensemblesv2.bash
wait $!
validate_ensemblesv2.bash
wait $!
calc_min_rmsd_metricsv2.bash
wait $!
calc_pcv2.bash
wait $!
popd
```

Look at the final table in `prelim_scoring/calc_pcv2.log`