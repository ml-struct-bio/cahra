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

Note that if you only want to install dependencies, for comp-het or pose-entanglement assessment, there are `requirements.txt` files in the corresponding subdirectories.

# For developers
pip install -e .[dev]
pre-commit install
```
The `pre-commit` is set out in `.pre-commit-config.yaml`. `pytest` is manually commented out due to lack of any tests.

## Compositional Heterogeneity Challenge Assessment
Instructions to run comp-het challenge assessment can be found in `comp-het/README.md`.

## Conformational Heterogeneity Challenge Assessment

Code for running this via a set of command line tools is included in `src` and scripts for reproducing on entries are in `conf-het/scripts`.

Once pip-installed (see above), command line interface utilities and options available for conf-het can be accessed via:
`cahra-metrics -h`

Instructions to run conf-het challenge assessment can be found in `conf-het/README.md`.

## Pose-Entanglement Heterogeneity Challenge Assessment
Instructions to run pose-entanglement challenge assessment can be found in `comp-het/README.md`.
