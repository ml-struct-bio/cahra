#!usr/bin/bash

python3 batch_reorder_pdb.py \
  --reference /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/6XOU.pdb \
  --mirror-base /mnt/clathrin/CAHRA-results-analysis \
  --output-root /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/reorder-all-trajs/output \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/extract_frames_6xou_unsteered_0-200ns_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/curious-fulmar-1/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/curious-fulmar-2/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/poised-bobcat-1/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/poised-bobcat-2/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/poised-bobcat-3/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/exotic-stingray-1/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/exotic-stingray-2/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/persuasive-penguins-1/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/persuasive-penguins-2/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/persuasive-penguins-3/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/persuasive-penguins-4/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/persuasive-penguins-5/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/kiharalab-1/confs_pdb \
  /mnt/clathrin/CAHRA-results-analysis/pdbfiles-to-reformat/prelim_submissions/kiharalab-2/confs_pdb \
