# Pose Entanglement Challenge Evaluation

Create a directory to put output plots: `mkdir plots`.

Compute the pose entanglement scoring metrics with the command: `python scoring.py ground_truth_star_fname=... pre_wrangled_submissions_fname=...` 
The ground-truth STAR file can be download via `wget "https://ftp.ebi.ac.uk/empiar/world_availability/13991/data/source.star"`.

For example: `mkdir -p plots; python scoring.py ground_truth_star_fname="source.star" pre_wrangled_submissions_fname="wrangled_submissions_20261008.pkl"`


