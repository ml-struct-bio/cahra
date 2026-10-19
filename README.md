# The Community-Wide Assessment of Cryo-EM Heterogeneous Reconstruction Algorithms (CAHRA)
Please use the [Discussions page](https://github.com/ml-struct-bio/cahra/discussions) for challenge-related discourse. For more information, please see the [challenge website](https://heterogeneity.notion.site/challenge).

* ADD LOCATION OF DATASETS INCLUDING ZENODO CONF-HET DATA
* ADD LOCATION OF ENTRIES
* ADD LOCATION OF THIS REPO (ISSUE TRACKER)
* ADD LOCATION OF PREPRINT

## Installation

Install [miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file#unix-like-platforms-macos-linux--wsl)

```
conda create -n cahra python=3.14
conda activate cahra
git clone git@github.com:ml-struct-bio/cahra.git
cd cahra

# For users
pip install .

# For developers
pip install -e .[dev]
pre-commit install
```

## Compositional Heterogeneity Challenge Assessment

## Pose-Entanglement Heterogeneity Challenge Assessment

## Conformational Heterogeneity Challenge Assessment

Code for running this via a set of command line tools is included in `src` and scripts for reproducing on entries are in `conf-het/scripts`.

The `pre-commit` is set out in `.pre-commit-config.yaml`. Currently `pytest` is manually commented out due to lack of any tests.



### Conf-het metrics CLI

Once pip-installed, command line interface utilities and options available for conf-het can be accessed via:
`cahra-metrics -h`

###  Running the conf-het metrics

1. Download the conformations
2. Put entries into directory `scoring/prelim_submissions/<entrylabel>/confs`
3. Copy scripts to `scoring` and run as below:
```
mkdir -p scoring
ADD COPYING OF DATA FROM ZENODO
cp conf-het/scripts/*sh scoring/
conda activate cahra # if not already active
pushd scoring
bash convert_ensembles.bash
wait $!
bash validate_ensembles.bash
wait $!
bash calc_min_rmsd_metrics.bash
wait $!
bash calc_pcv.bash
wait $!
ADD PLOTTING OF RMSF
bash make_ensemble_pdb.bash
wait $!
bash run_ensemble_maps.bash
wait $!
bash run_min_rmsd_displacement.sh # optional
wait $!
popd
```

To run stereochemistry plot generation routine, need to first grab 6xou.pdb
and add a subdir to metrics directory (e.g. `mkdir -p scoring/metric_calculation_label/orig_6xou/pdbs/`, `cp 6xou.pdb scoring/metric_calculation_label/orig_6xou/pdbs/`). 
Then run:

```
cp conf-het/scripts/*.sh scoring/
pushd scoring
bash run_gt_pred_rota_rama_clash.sh
wait $!
bash run_stereochem_analysis.bash 
wait $!
popd
```
Look at the final table in `scoring/calc_pcv.log`

Metrics calculated (and plots), along with aligned ensembles will be stored in `metric_calculation` directory


