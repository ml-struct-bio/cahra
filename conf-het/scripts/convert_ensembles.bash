#!/usr/bin/bash

# If needed, for the GT pdb files, get rid of non-standard residue names and update chain label,
# including removing hydrogens
# Note that this is already done for GT conformations provided as part of Zenodo/EMPIAR dataset(s)
# To rerun use l8-l24

LOG_DIR="convert_ensembles_logs"
GT_DIR="submissions/extract_frames_6xou_unsteered_0-200ns"
OUT_DIR="submissions/extract_frames_6xou_unsteered_0-200ns"
ORIG_EXT="pdb"
NEW_EXT="pdb"
UPDATE_CHAIN="X:A"

mkdir $LOG_DIR

pushd $GT_DIR
for pdbfile in *.pdb; do sed -i 's/GLH/GLU/g' ${pdbfile}; done;
wait $!
for pdbfile in *.pdb; do sed -i 's/HID/HIS/g' ${pdbfile}; done;
wait $!
popd
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $GT_DIR -o $OUT_DIR -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_gt_ensemble.log" 2>&1 &
wait $!


# Now iterate through each submitted ensemble and correct chains, atom labels etc which were found via validate_ensemble.bash
# Note that submitted and GT conformations provided as part of Zenodo/EMPIAR dataset(s) are already corrected using the following
# To rerun use l31 onwards

# Convert poised-bobcat-1 mmcifs to pdb so we can load using pdb file as mda.Universe topology
UPDATE_CHAIN=":A"
ORIG_EXT="cif"
NEW_EXT="pdb"
indir="submissions/poised-bobcat-1/confs_cif"
outdir="submissions/poised-bobcat-1/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_poised-bobcat-1_ensemble.log" 2>&1 &
wait $!

# convert poised-bobcat-2 mmcifs to pdb so we can load using pdb file as mda.Universe topology
UPDATE_CHAIN=":A"
indir="submissions/poised-bobcat-2/confs_cif"
outdir="submissions/poised-bobcat-2/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_poised-bobcat-2_ensemble.log" 2>&1 &
wait $!

# convert poised-bobcat-3 mmcifs to pdb so we can load using pdb file as mda.Universe topology
UPDATE_CHAIN=":A"
indir="submissions/poised-bobcat-3/confs_cif"
outdir="submissions/poised-bobcat-3/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_poised-bobcat-3_ensemble.log" 2>&1 &
wait $!

# exotic-stingray-1 submission has a model which appears to have the hydrogens left in it
# updating here
# UPDATED TO USE MODEL POST-H-REMOVAL (ie: confs directory instead of confs-1withH)
UPDATE_CHAIN=":A"
indir="submissions/exotic-stingray-1/confs_cif"
outdir="submissions/exotic-stingray-1/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_exotic-stingray-1_ensemble.log" 2>&1 &
wait $!

# exotic-stingray-2 submission has a couple of models which appears to have the hydrogens left in it
# updating here
UPDATE_CHAIN=":A"
indir="submissions/exotic-stingray-2/confs_cif"
outdir="submissions/exotic-stingray-2/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_exotic-stingray-2_ensemble.log" 2>&1 &
wait $!

# orange submission has models with author-assigned chain IDs of C instead of A, so change them here
indir="submissions/curious-fulmar-1/confs-C-labelled"
outdir="submissions/curious-fulmar-1/confs_pdb"
UPDATE_CHAIN="C:A"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_curious-fulmar-1_ensemble.log" 2>&1 &
wait $!
# these mmcifs have labels OT1 and OT2 in the C-terminus residue
pushd "submissions/curious-fulmar-1/confs_pdb" # $out_dir
wait $!
# replacing OT1 with O
for pdbfile in *.pdb; do sed -i 's/OT1/O  /g' ${pdbfile}; done; # for pdbs the spaces after O matter
wait $!
# replacing OT2 with OXT
for pdbfile in *.pdb; do sed -i 's/OT2/OXT/g' ${pdbfile}; done;
wait $!
popd

# convert curious-fulmar-2 to pdbs
indir="submissions/curious-fulmar-2/confs_cif"
outdir="submissions/curious-fulmar-2/confs_pdb"
UPDATE_CHAIN="C:A"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_curious-fulmar-2_ensemble.log" 2>&1 &
wait $!
# these mmcifs have labels OT1 and OT2 in the C-terminus residue
pushd "submissions/curious-fulmar-2/confs_pdb" # $out_dir
wait $!
# replacing OT1 with O
for pdbfile in *.pdb; do sed -i 's/OT1/O  /g' ${pdbfile}; done; # for pdbs the spaces after O matter
wait $!
# replacing OT2 with OXT
for pdbfile in *.pdb; do sed -i 's/OT2/OXT/g' ${pdbfile}; done;
wait $!
popd

# convert persuasive-penguins-1 to pdbs
UPDATE_CHAIN=":A"
indir="submissions/persuasive-penguins-1/confs_cif"
outdir="submissions/persuasive-penguins-1/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_persuasive-penguins-1_ensemble.log" 2>&1 &
wait $!

# convert persuasive-penguins-2 to pdbs
UPDATE_CHAIN=":A"
indir="submissions/persuasive-penguins-2/confs_cif"
outdir="submissions/persuasive-penguins-2/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_persuasive-penguins-2_ensemble.log" 2>&1 &
wait $!

# convert persuasive-penguins-3 to pdbs
UPDATE_CHAIN=":A"
indir="submissions/persuasive-penguins-3/confs_cif"
outdir="submissions/persuasive-penguins-3/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_persuasive-penguins-3_ensemble.log" 2>&1 &
wait $!

# convert persuasive-penguins-4 to pdbs
UPDATE_CHAIN=":A"
indir="submissions/persuasive-penguins-4/confs_cif"
outdir="submissions/persuasive-penguins-4/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_persuasive-penguins-4_ensemble.log" 2>&1 &
wait $!

# convert persuasive-penguins-4 to pdbs
UPDATE_CHAIN=":A"
indir="submissions/persuasive-penguins-5/confs_cif"
outdir="submissions/persuasive-penguins-5/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_persuasive-penguins-5_ensemble.log" 2>&1 &
wait $!

# convert kiharalab-1 to pdbs
UPDATE_CHAIN=":A"
indir="submissions/kiharalab-1/confs_cif"
outdir="submissions/kiharalab-1/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_kiharalab-1_ensemble.log" 2>&1 &
wait $!
# also need to add a OXT atom to the final residue in kiharalab-1
cahra-metrics add_missing_oxt -v -pe $outdir -e $NEW_EXT > "${LOG_DIR}/add_missing_oxt_kiharalab-1.log" 2>&1 &
wait $!
# add_missing_oxt writes out gemmi.Structure as pdb only, so convert these back to <fname>.$NEW_EXT files
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $outdir -o $outdir -oe "pdb" -e $NEW_EXT -v  > "${LOG_DIR}/convert_kiharalab-1_ensemble_afterOXTaddition.log" 2>&1 &

# convert kiharalab-2 to pdbs
UPDATE_CHAIN=":A"
indir="submissions/kiharalab-2/confs_cif"
outdir="submissions/kiharalab-2/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > "${LOG_DIR}/convert_kiharalab-2_ensemble.log" 2>&1 &
wait $!
# TODO: CHECK IF THE BELOW IS NEEDED FOR kiharalab-2!!!!
# also need to add a OXT atom to the final residue in kiharalab-2
cahra-metrics add_missing_oxt -v -pe $outdir -e $NEW_EXT > "${LOG_DIR}/add_missing_oxt_kiharalab-2.log" 2>&1 &
wait $!
# add_missing_oxt writes out gemmi.Structure as pdb only, so convert these back to <fname>.$NEW_EXT files
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $outdir -o $outdir -oe "pdb" -e $NEW_EXT -v  > "${LOG_DIR}/convert_kiharalab-2_ensemble_afterOXTaddition.log" 2>&1 &

