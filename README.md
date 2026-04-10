# cahra-internal
The Community-Wide Assessment of Cryo-EM Heterogeneous Reconstruction Algorithms (CAHRA)


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


## Merging of CAHRA-metrics repository into cahra-internal
The simplest way to integrate the CAHRA-metrics (conf-het metrics) repository into the cahra-internal repo seems to be to merge the repositories (rather than using submodules).

This keeps the history of both and allows further development to be done in cahra-internal only. Following [this example](https://stackoverflow.com/questions/1425892/how-do-you-merge-two-git-repositories) and solving merge conflicts in the license and README.
```
git remote add CAHRA-metrics /home/jg/Software/CAHRA-metrics
git fetch CAHRA-metrics
git merge --allow-unrelated-histories CAHRA-metrics/main

```