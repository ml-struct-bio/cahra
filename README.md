# The Community-Wide Assessment of Cryo-EM Heterogeneous Reconstruction Algorithms (CAHRA)
Please use the [Discussions page](https://github.com/ml-struct-bio/cahra/discussions) for challenge-related discourse. For more information, please see the [challenge website](https://heterogeneity.notion.site/challenge).

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
The `pre-commit` is set out in `.pre-commit-config.yaml`. `pytest` is manually commented out due to lack of any tests.

## Compositional Heterogeneity Challenge Assessment

## Pose-Entanglement Heterogeneity Challenge Assessment

## Conformational Heterogeneity Challenge Assessment

Code for running this via a set of command line tools is included in `src` and scripts for reproducing on entries are in `conf-het/scripts`.

Once pip-installed (see above), command line interface utilities and options available for conf-het can be accessed via:
`cahra-metrics -h`

### Set Up To Re-run Conf-Het Metrics

Firstly create a directory to calculate the metrics within and also download the predicted & ground truth ensembles.
```
mkdir -p scoring
pushd scoring
wget -O submissions.zip "https://zenodo.org/records/23497320/files/submissions.zip?download=1"
wait $!
unzip submissions.zip
wait $!
rm submissions.zip
popd
```

The copy the scripts used for metric calculation:
```
cp conf-het/scripts/*.bash scoring/
```

### Run The Conf-Het Metrics 
Remember to check the log files after executing each command.

The run SAC & LFC calculation, alongside outputting min-RMSD-aligned ensembles, run:
```
pushd scoring
bash calc_min_rmsd_metrics.bash
wait $!
popd
```


The outputs (by default in the `metric_calculation_label` directory) are the min-RMSD-aligned ensembles. The GT ensemble is in the `gt_confs` subdirectory and predicted ensembles are in the `pred_confs` subdirectory. In further subdirectories, a number of figures are saved for each predicted ensemble, along with a json tracking the results (and allowing reuse of intermediate results, such as min-RMSD values and indices).

To then calculate the PCV, we run another script as follows:
```
pushd scoring
bash calc_pcv.bash
wait $!
popd
```
Plots over all predicted ensembles will be added (by default to the `metric_calculation_label/pred_confs/` directory), with a `all_metrics.csv` file, containing the metric values and relative rankings for the entries. The `pred_conf/<label>/rms_metrics.json` is also updated to `pred_conf/<label>/all_metrics.json` for each predicted ensemble.

More concentrated RMSF plots overlapping a subset of predicted ensembles can be created using the scripts in `conf-het/scripts/plot_rmsf` directory, with example bash scripts therein. Before running these, copy the `plot_rmsf` directory to the `pred_confs directory` as so: `cp conf-het/scripts/plot_rmsf scoring/metric_calculation_label/pred_confs/`. This will allow the script to find the correct data.

### Create Ensemble PDBs, Ensemble Maps and Difference Maps 

To create an ensemble PDB file for each of the aligned ensembles run the corresponding script:
```
pushd scoring
make_ensemble_pdb.bash
popd
```
The ensemble PDB files go from the most open conformation (according to collective variable) to the most closed and are put into the `metric_calculation_label/pred_confs/<label>/ensemble_pdb/` directory by default. 

To create a density map from each of the aligned ensembles (as well as a difference map compared to the GT ensemble), run the corresponding script:
```
pushd scoring
bash run_ensemble_maps.bash
popd
```
The output density maps are output in `metric_calculation_label/pred_confs/<label>/ensemble_map` and the difference maps in `metric_calculation_label/pred_confs/<label>/ensemble_diffmap`. The difference maps are relative to the GT map, which will be generated at: `metric_calculation_label/gt_confs/ensemble_map`.


### Run Stereochemistry Analysis

To run stereochemistry plot generation routine, need to first grab 6xou.pdb from the PDB and add a subdir to metrics directory (e.g. `mkdir -p scoring/metric_calculation_label/orig_6xou/pdbs/ ; cp 6xou.pdb scoring/metric_calculation_label/orig_6xou/pdbs/`). These scripts require a CCP4 install and are set up to expect the sourcing script in the location: `/opt/xtal/ccp4-9/bin/ccp4.setup-sh`. Adjust as appropriate in `run_rota_rama_clash_pred.bash`, `run_rota_rama_clash_gt.bash` and `run_rota_rama_clash_6xou.bash`.
Then run:
```
cp conf-het/scripts/*.bash scoring/
pushd scoring
bash run_gt_pred_rota_rama_clash.bash
wait $!
bash run_stereochem_analysis.bash 
wait $!
popd
```


### (Extra) A Note On Validating Ensembles

**The atomic model ensembles downloaded above include PDB versions which have already had their formatting and atomic model ordering checked. This does not need to be repeated. Details of this processing are included in this section.**

Note that some software, particularly MD simulation programs do not enforce conventional atomic ordering in PDB files.
Note also that the scoring metrics utilise libraries which depend on PDB files and on topology files, such as MDAnalysis. This means that MMCIF files submitted needed to be converted to PDB before use.
It also means that any submitted (or MD-simulation-derived) PDB files need to have all the atoms in each residue in the same order.

To confirm that this is the case, the atomic models in each ensemble can be checked for self-similarity and for similarity to another ensemble using the `cahra-metrics validate_ensembles` CLI utility. This was done using the `validate_ensembles.bash` script.
If you want to rerun it, you can do so via:
```
pushd scoring
bash validate_ensembles.bash
popd
```
Remember to check the output log files after running!

The CLI also offers a command to update the chain label, remove hydrogens and transform MMCIF files to PDB: `cahra-metrics convert_ensemble`. An example script including some of the original processing done on the submitted ensembles is `convert_ensembles.bash`. An example snippet to rerun it is follows, but again note that the atomic models provided are already fully pre-processed & validated.
```
pushd scoring
bash convert_ensembles.bash
popd
```
Remember to check the output log files after running!

Also made available are scripts for checking that all PDB files have atoms in each residue in the same order as a reference file (`6XOU.pdb` in the example). These can be found in `conf-het/scripts/reorder-pdbs`, which an example of how to run them on the conf-het ensembles in `conf-het/scripts/reorder-pdbs/reorder-all-trajs.bash`.

### (Extra) The Entire Conf-Het Processing Pipeline

Below is how the scoring was run during competiton. Note that the `wait $!` commands are essential if running the whole scoring & analysis pipeline.

1. Download the conformations
2. Put (MMCIF) entries into directory `scoring/submissions/<entrylabel>/confs`
3. Copy scripts to `scoring` and run as below:
```
mkdir -p scoring
cp conf-het/scripts/*.bash scoring/
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
bash make_ensemble_pdb.bash
wait $!
bash run_ensemble_maps.bash
wait $!
popd
```