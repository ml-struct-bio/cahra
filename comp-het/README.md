# Compositional Heterogeneity Challenge Evaluation

Evaluate particle classifications for the CAHRA compositional heterogeneity
challenge. The tool matches submitted particles to ground-truth species and
exposure labels, reports the composition of each class, and calculates a purity
score out of 100.

The four targets are PKM2, GDH, ADH and Aldolase. Gamma globulin and Beta-Amylase
are decoys.

## Install and run

Install the Python dependencies and run the bundled control:

```bash
cd comp-het
<create and/or activate python3 environment>
python -m pip install -r requirements.txt
python -m comphet
```

By default, the tool reads submission folders directly inside `comp-het/`, uses
`source_labels.pkl.gz` as the reference, and writes to `comp-het/results/`.
The bundled `amaranth_rattlesnake/` control scores 72.4/100 in the printed table.

To evaluate your own results, put one particle file per class in a submission
folder, then run the same command:

```text
comp-het/
├── amaranth_rattlesnake/
│   ├── pkm2_class1.star
│   ├── gdh_class2.star
│   └── adh_class3.star
└── my_results/
    ├── class1.cs
    ├── class2.cs
    └── class3.cs
```

Each immediate subfolder containing `.cs` or `.star` files is treated as a
submission. Hidden folders and deeper subfolders are skipped. Plotting supports
up to nine submissions per run.

To use other input and output locations:

```bash
python -m comphet --results-dir /path/to/submissions --out /path/to/output
python -m comphet --help
```

`--results-dir` must contain the submission folders, rather than the class files
themselves. `--dpi` sets figure resolution (default: 300).

## Input files

Use CryoSPARC `.cs` or RELION `.star` particle metadata, with one file per class.
Particle images are not read.

| Format | Required for evaluation | Optional diagnostic fields |
|---|---|---|
| `.cs` | `blob/path`, `blob/idx` | `ctf/df1_A`, `ctf/df2_A` |
| `.star` | `rlnImageName` | `rlnDefocusU`, `rlnDefocusV` |

For a STAR file with multiple data blocks, particle metadata must be in the
`particles` block. Stack names must retain the challenge's numeric stack IDs,
such as `0022.mrcs`. CryoSPARC UID prefixes are stripped from stack names before
matching.

The class label is the full filename without its extension. For example,
`pkm2_class1.star` gets the label `pkm2_class1`. Each class must have one file,
and labels must be unique within a submission. Files such as `class1.cs` and
`class1.star` have the same label and cause evaluation to stop with an error.
Separate submissions may use the same labels. Species are assigned from the
matched particles, not from filenames.

## Scoring

Each class is assigned to the species with the largest particle count, summed
over both exposure subsets. This does not require a majority. For each target,
the highest-purity class assigned to that target is scored:

```text
target_score = 25 × matched particles of the target species in the selected class
                  / all matched particles in the selected class
```

A target with no class assigned to it scores zero. The four target scores sum
to a maximum of 100. Other species, including decoys, reduce the selected
class's purity. There is no separate decoy penalty.

Duplicate and unmatched particle identifiers produce warnings in the console
and `summary_qc.csv`; evaluation continues. Repeated identifiers count once
within each file, and only matched particles contribute to composition and
scores. Particles shared across class files within the same submission also
produce a warning, but each class is evaluated separately. Overlap between
different submissions is expected and is not flagged.

The score does not reward recovery of a larger fraction of the target
population. A pure class containing only 10% of a target's particles still
earns 25 points. Recovery is reported separately as `recall_pct`.

## Outputs

The score table is printed, and six files are written to the output folder:

| File | Contents |
|---|---|
| `summary_scores.csv` | Rank, selected class, purity and points for each target, and total score per submission |
| `summary_composition.csv` | Matched counts, `pct_of_class` and `recall_pct` per class, species and exposure |
| `summary_qc.csv` | Match counts, duplicate row counts and percentages, counts shared with other class files, warnings, metadata errors and defocus diagnostics per file |
| `composition_heatmap.png` | Species percentages for each class |
| `submission_small_multiples.png` | Class composition as stacked bars, with one panel per submission |
| `exposure_split.png` | Observed F30/L30 split and the split expected from each class's species composition |

The `team` column contains the submission folder name. In the composition
table, `recall_pct` is the row's count divided by the total number of particles
of that species in the reference, multiplied by 100. Sum it over F30 and L30
to get recall for a species within a class.

Files with no matches appear only in the QC table; submissions with no matches
have no score. If nothing matches anywhere in the run, it still completes:
the score and composition CSVs contain headers only, and the figures state
that no particles matched. These outputs replace results from earlier runs.

Missing, invalid or suspect identifiers remain a separate error condition:
evaluation writes the QC diagnostics and stops before scoring. In that case,
other CSVs and figures from earlier runs are left unchanged.

## Particle matching

The primary key is the particle's index within its stack plus the numeric stack
ID. Directory paths are ignored, and CryoSPARC's zero-based indices are
converted to one-based indices.

Defocus matching is diagnostic only. The tool checks defocus values against the
reference when a file is missing the required identifier columns, contains
invalid identifiers, or has at least 1,000 rows referencing at most three stacks.
The last condition is treated as a sign of damaged stack metadata, even if the
identifiers can be parsed and found in the reference.

`rlnDefocusU` and `rlnDefocusV` (or the corresponding CryoSPARC fields) may remain
intact when stack or image-index metadata has been lost or corrupted. Both
defocus fields are needed for this check; missing or unusable values are reported.

The diagnostic rounds defocus U and V to three decimal places and counts
matching unique keys. These values are not unique for every reference particle,
so a defocus match does not establish a particle's identity. Diagnostic matches
are never included in composition, recall, scores or figures, even if every
defocus key matches.

The tool prints the reason for the diagnostic and its results, saves
`summary_qc.csv`, and stops the entire run before scoring. Correct the particle
identifiers before rerunning. In the QC table, `match_method` is `path_idx` for
normal matching or `defocus_diagnostic` when a diagnostic is attempted.
`primary_error` records the identifier problem, and the `n_defocus_*` columns
record diagnostic counts separately from accepted particle matches.

Files that pass these metadata checks use the primary key. Unmatched identifiers
alone do not trigger defocus matching.

## Reference data

`source_labels.pkl.gz` is a compressed pandas DataFrame containing
`primary_key`, `defocus_key`, `species` and `exposure` for 1,613,664 particles.
It is sufficient to run the evaluation without the full reference STAR file.

```python
import pandas as pd
labels = pd.read_pickle("source_labels.pkl.gz")
```

The full `source.star` is available with
[EMPIAR-13973](https://empiar.pdbj.org/en/entry/13973/). It can be downloaded directly with `wget "https://ftp.ebi.ac.uk/empiar/world_availability/13973/data/Particles/particles.star"` or similar. Use it directly with
`--reference` (also accepts `source.star.gz`).

```bash
python -m comphet --reference /path/to/source.star
```
You can regenerate the compact file via:
```bash
python -m comphet.make_reference /path/to/source.star source_labels.pkl.gz
```

The bundled reference contains:

| Species | Particles | % of dataset | F30 / L30 (%) |
|---|---:|---:|---|
| PKM2 | 656,101 | 40.7% | 50 / 50 |
| ADH | 239,594 | 14.8% | 80 / 20 |
| Aldolase | 213,305 | 13.2% | 80 / 20 |
| GDH | 185,201 | 11.5% | 80 / 20 |
| Gamma-Globulin (decoy) | 162,065 | 10.0% | 100 / 0 |
| Beta-Amylase (decoy) | 157,398 | 9.8% | 100 / 0 |
| Total | 1,613,664 | | |

F30 and L30 are the two exposure subsets; percentages above are rounded.
In `exposure_split.png`, the expected split is calculated by weighting each
species' reference exposure ratio by its abundance in the class. A difference
between the observed and expected splits indicates an exposure imbalance
relative to that baseline.
