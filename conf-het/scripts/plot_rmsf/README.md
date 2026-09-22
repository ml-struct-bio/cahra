# RMSF plots

Run with the existing Python environment (Gemmi, NumPy, and Matplotlib):

```bash
python plot_rmsf.py
```

The default input directory is the parent of this script's directory. The
reference is `../../gt_confs/pdbs/0_6xou_0-200ns_115.pdb`, resolved relative to the
script. Outputs go to `rmsf_plots/` next to the script, regardless of the current
working directory:

- `rmsf_all_ensembles.{png,pdf}`: every prediction plus one black GT curve.
- `rmsf_absolute_difference_all_ensembles.{png,pdf}`: every prediction's absolute
  difference from GT.
- `individual/<label>_rmsf.{png,pdf}`: one prediction and GT.
- `individual/<label>_absolute_difference.{png,pdf}`: that prediction's absolute
  difference from GT.
- `rmsf_per_residue.csv`: residue identifiers, atom counts, mean RMSFs, and
  absolute differences, in angstroms.

For selected ensembles or a different output directory:

```bash
python plot_rmsf.py --ensembles poised-bobcat-1 poised-bobcat-2 --output-dir rmsf_selected
```

Atomic RMSFs are averaged arithmetically within each residue. Differences are
`abs(predicted residue mean - GT residue mean)`. This is neither an average of
absolute atomic differences nor an RMSF calculated from residue-centroid
trajectories. Arrays must use the same atom order as the first model of the
reference PDB; the script checks atom counts, finite/nonnegative values, and
identical GT arrays across submissions. Atom counts alone cannot verify order.
Folders with neither RMSF array are skipped during automatic discovery; a
selected submission missing RMSFs causes an error.

The default x-axis uses PDB residue numbers. Shading follows the supplied ranges:
D1 10–265, LNK 539–541, and D4 832–1014. Use `--no-domains` to remove it, or repeat
`--domain LABEL:START:END:COLOR` to replace the ranges. For a reference with multiple
chains or insertion codes, use `--x-axis index --no-domains`; indices are 1-based
in reference traversal order, and the CSV retains chain and residue identifiers.
RMSF plots share one y-axis range; difference plots share a separate range, so
individual ensemble plots are comparable. Rerunning overwrites matching outputs.

See `python plot_rmsf.py --help` for input paths, output formats, and DPI options.
