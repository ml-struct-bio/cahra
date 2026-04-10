# cahra-internal
The Community-Wide Assessment of Cryo-EM Heterogeneous Reconstruction Algorithms (CAHRA)

## Conformational Heterogeneity Challenge Assessment

Code for running this via a set of command line tools is included in `src` and scripts for reproducing on preliminary entries are in `conf-het/scripts`.

The `pre-commit` is set out in `.pre-commit-config.yaml`. Currently `pytest` is manually commented out due to lack of any tests.

### Installation
Install [miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file#unix-like-platforms-macos-linux--wsl)

```
conda create -n cahra-internal python=3.14
conda activate cahra-internal
git clone git@github.com:ml-struct-bio/cahra-internal.git
pip install -e .[dev]
pre-commit install
```

###  Running the conf-het metrics

1. Download the conformations
2. Put entries into directory `prelim_scoring/prelim_submissions/<entrylabel>/confs`
3. Copy scripts to `prelim_scoring` and run as below:
```
cp conf-het/scripts/*.bash prelim_scoring
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


### Merging of CAHRA-metrics (conf-het) repository into cahra-internal
The simplest way to integrate the CAHRA-metrics (conf-het metrics) repository into the cahra-internal repo seems to be to merge the repositories (rather than using submodules).

This keeps the history of both and allows further development to be done in cahra-internal only. Following [this example](https://stackoverflow.com/questions/1425892/how-do-you-merge-two-git-repositories) and solving merge conflicts in the license and README.
```
git remote add CAHRA-metrics /home/jg/Software/CAHRA-metrics
git fetch CAHRA-metrics
git merge --allow-unrelated-histories CAHRA-metrics/main

```