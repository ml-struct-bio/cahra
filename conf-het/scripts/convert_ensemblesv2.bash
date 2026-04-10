#!/usr/bin/bash

# convert gt pdbs to mmcif files
# including removing hydrogens

GT_DIR="extract_frames_6xou_unsteered_0-200ns"
OUT_DIR="extract_frames_6xou_unsteered_0-200ns_pdb"
ORIG_EXT="pdb"
NEW_EXT="pdb"
UPDATE_CHAIN="X:A"

pushd $GT_DIR
for pdbfile in *.pdb; do sed -i 's/GLH/GLU/g' ${pdbfile}; done;
wait $!
for pdbfile in *.pdb; do sed -i 's/HID/HIS/g' ${pdbfile}; done;
wait $!
popd
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $GT_DIR -o $OUT_DIR -oe $ORIG_EXT -e $NEW_EXT -v > convert_gt_ensemble.log 2>&1 &
wait $!

# convert lemon mmcifs to pdb so we can load using pdb file as mda.Universe topology
UPDATE_CHAIN=":A"
ORIG_EXT="cif"
NEW_EXT="pdb"
indir="prelim_submissions/lemon1/confs"
outdir="prelim_submissions/lemon1/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > convert_lemon1_ensemble.log 2>&1 &
wait $!

# lime submission has a model which appears to have the hydrogens left in it
# updating here
indir="prelim_submissions/lime1/confs-1withH"
outdir="prelim_submissions/lime1/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > convert_lime1_ensemble.log 2>&1 &
wait $!

# orange submission has models with author-assigned chain IDs of C instead of A, so change them here
indir="prelim_submissions/orange1/confs-C-labelled"
outdir="prelim_submissions/orange1/confs_pdb"
UPDATE_CHAIN="C:A"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > convert_orange1_ensemble.log 2>&1 &
wait $!
# these mmcifs have labels OT1 and OT2 in the C-terminus residue
pushd "prelim_submissions/orange1/confs_pdb" # $out_dir
wait $!
# replacing OT1 with O
for pdbfile in *.pdb; do sed -i 's/OT1/O  /g' ${pdbfile}; done; # for pdbs the spaces after O matter
wait $!
# replacing OT2 with OXT
for pdbfile in *.pdb; do sed -i 's/OT2/OXT/g' ${pdbfile}; done;
wait $!
popd

# convert kiwi1 to pdbs
UPDATE_CHAIN=":A"
indir="prelim_submissions/kiwi1/confs"
outdir="prelim_submissions/kiwi1/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > convert_kiwi1_ensemble.log 2>&1 &
wait $!

# convert banana1 to pdbs
indir="prelim_submissions/banana1/confs"
outdir="prelim_submissions/banana1/confs_pdb"
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $indir -o $outdir -oe $ORIG_EXT -e $NEW_EXT -v > convert_banana1_ensemble.log 2>&1 &
wait $!
# also need to add a OXT atom to the final residue in banana1
cahra-metrics add_missing_oxt -v -pe $outdir -e $NEW_EXT > add_missing_oxt_banana1.log 2>&1 &
# add_missing_oxt writes out gemmi.Structure as pdb only, so convert these back to <fname>.$NEW_EXT files
cahra-metrics convert_ensemble -c $UPDATE_CHAIN -i $outdir -o $outdir -oe "pdb" -e $NEW_EXT -v  > convert_banana1_ensemble_afterOXTaddition.log 2>&1 &

