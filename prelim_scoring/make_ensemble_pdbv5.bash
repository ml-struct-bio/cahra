#!/usr/bin/bash

EXT="pdb"

# get ensemble pdb for the aligned GT ensemble 
ensemble_dir="metric_calculation/gt_confs/pdbs"
output_ensemble_dir="metric_calculation/gt_confs/ensemble_pdb"
fname="ensemble_pdb_GT.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="lemon1"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="lemon2"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="lemon3"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="lime1"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="lime2"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="orange1"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="orange2"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiwi1"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiwi2"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiwi3"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiwi4"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="kiwi5"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="banana1"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

# get ensemble pdb for the min-rmsd-aligned pred ensembles
label="banana2"
ensemble_dir="metric_calculation/pred_confs/${label}/pdbs"
output_ensemble_dir="metric_calculation/pred_confs/${label}/ensemble_pdb"
fname="ensemble_pdb_${label}.pdb"
cahra-metrics make_ensemble_pdb -d $ensemble_dir -o $output_ensemble_dir -e $EXT -f $fname -v
wait $!

