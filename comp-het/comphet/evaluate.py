"""Entry point: score every submission folder and write the figures.

    python -m comphet
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd

from .matching import build_summary, find_submissions, load_reference
from .scoring import class_level_table, compute_target_scores, scores_display_table
from .viz import (plot_composition_heatmap, plot_exposure_split,
                  plot_submission_small_multiples)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_REFERENCE = REPO_ROOT / "source_labels.pkl.gz"


def run(results_dir: Path, reference: Path, out_dir: Path, dpi: int = 300) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    found = find_submissions(results_dir)
    if not found:
        raise SystemExit(
            f"No submissions found in {results_dir}.\n"
            "Each submission should be its own subfolder of .cs/.star files, e.g.\n"
            "    my_results/class1.cs, my_results/class2.cs")
    names = sorted({t for t, _, _ in found})
    print(f"Found {len(found)} class file(s) across {len(names)} submission(s): {', '.join(names)}")

    print(f"Loading reference labels from {reference.name} ...")
    ref = load_reference(reference)

    qc, comp = build_summary(results_dir, ref)
    qc_path = out_dir / "summary_qc.csv"
    qc.to_csv(qc_path, index=False)

    warnings = qc[qc["warning"] != ""]
    failed = qc[qc["primary_error"] != ""]
    if len(warnings):
        print("\nQC warnings:")
        for _, r in warnings.iterrows():
            print(f"  WARNING: {r['team']}/{r['file']}: {r['warning']}")
    if len(failed):
        print("\nMetadata diagnostics:")
        for _, r in failed.iterrows():
            print(f"  {r['team']}/{r['file']}: {r['primary_error']}")
            if pd.notna(r["n_defocus_unique"]):
                print(f"    Defocus diagnostic: {int(r['n_defocus_matched']):,} of "
                      f"{int(r['n_defocus_unique']):,} unique keys matched; "
                      f"{int(r['n_defocus_invalid_rows']):,} rows had unusable defocus values.")
            print(f"    {r['defocus_note']}")

    if len(failed):
        raise SystemExit(
            f"Evaluation stopped: {len(failed)} file(s) have missing, invalid or suspect particle identifiers.\n"
            f"QC and defocus diagnostics were written to {qc_path}.\n"
            "Correct the identifiers and rerun. No scores or figures were produced for this run; "
            "any existing results are unchanged.")
    class_tbl = class_level_table(comp)

    scores = compute_target_scores(comp, class_tbl)
    scores.insert(0, "rank", range(1, len(scores) + 1))

    comp.to_csv(out_dir / "summary_composition.csv", index=False)
    scores.to_csv(out_dir / "summary_scores.csv", index=False)

    order = scores["team"].tolist()  # figures run best-scoring first
    for fig, name in [
        (plot_composition_heatmap(class_tbl, order=order), "composition_heatmap.png"),
        (plot_submission_small_multiples(class_tbl), "submission_small_multiples.png"),
        (plot_exposure_split(comp, class_tbl, ref["totals_by_exposure"], order=order),
         "exposure_split.png"),
    ]:
        fig.savefig(out_dir / name, dpi=dpi, facecolor="#fcfcfb", bbox_inches="tight")

    print()
    if scores.empty:
        print("No scores available: no particles matched the reference. "
              "Result tables are empty and figures are marked accordingly.")
    else:
        print(scores_display_table(scores).to_string(index=False))
    print(f"\nWrote 3 figures and 3 CSVs to {out_dir}/")


def main() -> None:
    p = argparse.ArgumentParser(
        prog="python -m comphet",
        description="Score cryo-EM class assignments against the challenge ground truth.")
    p.add_argument("--results-dir", type=Path, default=REPO_ROOT,
                   help="folder holding one subfolder per submission (default: the repo root)")
    p.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE,
                   help="source_labels.pkl.gz, or a full source.star / source.star.gz")
    p.add_argument("--out", type=Path, default=REPO_ROOT / "results",
                   help="output folder for figures and CSVs (default: results/)")
    p.add_argument("--dpi", type=int, default=300)
    a = p.parse_args()
    run(a.results_dir, a.reference, a.out, a.dpi)


if __name__ == "__main__":
    main()
