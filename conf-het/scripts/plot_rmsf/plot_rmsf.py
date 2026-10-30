#!/usr/bin/env python3
"""Plot saved atomic RMSFs averaged per residue; no RMSFs are recomputed.

The RMSF arrays must follow the atom order of the reference PDB (first model).
An atom-count check cannot establish atom identity/order; use the same topology
as was used to calculate the RMSFs. Every reference atom contributes equally to
its residue's arithmetic mean. Absolute differences are calculated AFTER this
averaging: abs(mean(pred atom RMSFs) - mean(GT atom RMSFs)).

Requires gemmi, numpy, and matplotlib. Run `python plot_rmsf.py --help` for options.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import gemmi
import matplotlib

matplotlib.use("Agg")  # Save figures on systems without a display.
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
DOMAINS = (
    (10, 265, "D1", "orange"),
    (539, 541, "LNK", "blue"),
    (832, 1014, "D4", "green"),
    (151, 202, "Switch A", "violet"),
    (398, 413, "Switch B", "cadetblue"),
    (483, 527, "Switch C", "grey"),
)


def read_reference(path, x_axis):
    """Return residue metadata, contiguous atom slices, and plotting coordinates."""
    structure = gemmi.read_structure(str(path))
    if len(structure) == 0:
        raise ValueError(f"No models in reference: {path}")
    residues, starts, counts = [], [], []
    offset = 0
    for chain in structure[0]:
        for residue in chain:
            if not len(residue):
                continue
            residues.append(
                (
                    chain.name,
                    str(residue.seqid),
                    residue.name,
                    residue.seqid.num,
                )
            )
            starts.append(offset)
            counts.append(len(residue))
            offset += len(residue)
    if not residues:
        raise ValueError(f"No atoms in reference: {path}")
    if x_axis == "residue-number":
        x = np.array([r[3] for r in residues])
        if len({r[0] for r in residues}) != 1 or np.any(np.diff(x) <= 0):
            raise ValueError(
                "Reference has multiple chains, insertion codes, or "
                "non-increasing residue numbers. Use --x-axis index "
                "and --no-domains (or supply index-based --domain ranges)."
            )
    else:
        x = np.arange(1, len(residues) + 1)
    return residues, np.asarray(starts), np.asarray(counts), x


def read_ensembles(root, labels, atom_count):
    """Load selected predictions and verify their common ground-truth array."""
    paths = (
        [root / label / "results" / "all_metrics.json" for label in labels]
        if labels
        else sorted(root.glob("*/results/all_metrics.json"))
    )
    predictions, ground_truth = {}, None
    for path in paths:
        with path.open() as handle:
            data = json.load(handle)
        stats = data.get("ensemble", {}).get("ensemble_stats", {})
        if "pred_rmsf" not in stats and "gt_rmsf" not in stats and not labels:
            print(
                f"Skipping {path.parent.parent.name}: no RMSF arrays.",
                file=sys.stderr,
            )
            continue
        arrays = {}
        for key in ("pred_rmsf", "gt_rmsf"):
            if key not in stats:
                raise ValueError(
                    f"{path}: missing ensemble.ensemble_stats.{key}"
                )
            values = np.asarray(stats[key], dtype=float)
            if values.ndim != 1 or len(values) != atom_count:
                raise ValueError(
                    f"{path}: {key} has shape {values.shape}; expected "
                    f"({atom_count},) to match reference atom order."
                )
            if not np.isfinite(values).all() or np.any(values < 0):
                raise ValueError(
                    f"{path}: {key} contains nonfinite or negative RMSFs."
                )
            arrays[key] = values
        if ground_truth is None:
            ground_truth = arrays["gt_rmsf"]
        elif not np.array_equal(ground_truth, arrays["gt_rmsf"]):
            raise ValueError(
                f"{path}: GT RMSFs differ from the other ensembles."
            )
        predictions[path.parent.parent.name] = arrays["pred_rmsf"]
    if not predictions:
        raise ValueError(f"No ensemble RMSF arrays found under {root}")
    return predictions, ground_truth


def make_plot(
    x, curves, output, title, ylabel, xlabel, domains, formats, dpi, ymax
):
    """Write one figure with consistent colors and optional domain shading."""
    fig, ax = plt.subplots(figsize=(14, 5))
    for start, end, label, color in domains:
        ax.axvspan(start, end, color=color, alpha=0.15, zorder=0, label=label)
    for label, values, color, is_gt in curves:
        ax.plot(
            x,
            values,
            label=label,
            color=color,
            linewidth=2.2 if is_gt else 1.1,
            alpha=1.0 if is_gt else 0.85,
            zorder=4 if is_gt else 2,
        )
    ax.set(
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        xlim=(x[0] - 0.5, x[-1] + 0.5),
        ylim=(0, ymax),
    )
    ax.grid(axis="y", alpha=0.2)
    # List curves before domain patches in the external legend.
    handles, labels = ax.get_legend_handles_labels()
    n = len(domains)
    ax.legend(
        handles[n:] + handles[:n],
        labels[n:] + labels[:n],
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        fontsize=8,
    )
    fig.tight_layout()
    try:
        for extension in formats:
            fig.savefig(
                output.with_suffix(f".{extension}"),
                dpi=dpi,
                bbox_inches="tight",
            )
    finally:
        plt.close(fig)


def parse_domain(value):
    try:
        label, start, end, color = value.split(":")
        start, end = float(start), float(end)
        if not label or not np.isfinite([start, end]).all() or start >= end:
            raise ValueError
        if not matplotlib.colors.is_color_like(color):
            raise ValueError
        return start, end, label, color
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Use LABEL:START:END:COLOR with finite START < END, e.g. D1:10:265:orange"
        ) from None


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=HERE.parent,
        help="Directory containing submitted ensembles (default: ../)",
    )
    parser.add_argument(
        "--reference-pdb",
        type=Path,
        default=HERE.parent.parent
        / "gt_confs/pdbs/0_6xou_0-200ns_115_reordered.pdb",
    )
    parser.add_argument("--output-dir", type=Path, default=HERE / "rmsf_plots")
    parser.add_argument(
        "--ensembles",
        nargs="+",
        metavar="LABEL",
        help="Only plot these submission directory names (default: all)",
    )
    parser.add_argument(
        "--x-axis",
        choices=("residue-number", "index"),
        default="residue-number",
        help="PDB residue number or 1-based residue index",
    )
    domains = parser.add_mutually_exclusive_group()
    domains.add_argument("--no-domains", action="store_true")
    domains.add_argument(
        "--domain",
        action="append",
        type=parse_domain,
        help="Replace default shading; repeat LABEL:START:END:COLOR in x-axis units",
    )
    parser.add_argument(
        "--formats", nargs="+", choices=("png", "pdf"), default=["png", "pdf"]
    )
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    if args.ensembles and any(
        Path(s).name != s or s in (".", "..") for s in args.ensembles
    ):
        parser.error(
            "--ensembles must contain submission directory names, not paths"
        )

    try:
        residues, starts, counts, x = read_reference(
            args.reference_pdb, args.x_axis
        )
        predictions, gt = read_ensembles(
            args.input_dir, args.ensembles, int(counts.sum())
        )
        gt = np.add.reduceat(gt, starts) / counts
        predictions = {
            label: np.add.reduceat(values, starts) / counts
            for label, values in predictions.items()
        }
        differences = {
            label: np.abs(values - gt) for label, values in predictions.items()
        }
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(1, f"Error: {error}\n")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    individual_dir = args.output_dir / "individual"
    individual_dir.mkdir(exist_ok=True)
    labels = list(predictions)
    colors = {
        label: plt.get_cmap("tab20")(i % 20) for i, label in enumerate(labels)
    }
    domain_ranges = [] if args.no_domains else (args.domain or DOMAINS)
    xlabel = (
        "Residue number (PDB)"
        if args.x_axis == "residue-number"
        else "Residue index (1-based)"
    )
    rmsf_max = max(
        float(gt.max()), *(float(v.max()) for v in predictions.values())
    )
    diff_max = max(float(v.max()) for v in differences.values())
    common = dict(
        x=x,
        xlabel=xlabel,
        domains=domain_ranges,
        formats=args.formats,
        dpi=args.dpi,
    )
    gt_curve = ("Ground truth", gt, "black", True)

    def plot_rmsf(selected, output, title):
        make_plot(
            curves=[(s, predictions[s], colors[s], False) for s in selected]
            + [gt_curve],
            output=output,
            title=title,
            ylabel=r"RMSF ($\AA$)",
            ymax=1.05 * rmsf_max if rmsf_max else 1,
            **common,
        )

    def plot_difference(selected, output, title):
        make_plot(
            curves=[(s, differences[s], colors[s], False) for s in selected],
            output=output,
            title=title,
            ylabel=r"Absolute RMSF difference ($\AA$)",
            ymax=1.05 * diff_max if diff_max else 1,
            **common,
        )

    plot_rmsf(
        labels,
        args.output_dir / "rmsf_all_ensembles",
        "Per-residue RMSF: all ensembles",
    )
    plot_difference(
        labels,
        args.output_dir / "rmsf_absolute_difference_all_ensembles",
        "Absolute per-residue RMSF difference from ground truth",
    )
    for label in labels:
        plot_rmsf(
            [label],
            individual_dir / f"{label}_rmsf",
            f"Per-residue RMSF: {label}",
        )
        plot_difference(
            [label],
            individual_dir / f"{label}_absolute_difference",
            f"Absolute per-residue RMSF difference: {label}",
        )

    with (args.output_dir / "rmsf_per_residue.csv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "residue_index",
                "chain",
                "residue_id",
                "residue_name",
                "atom_count",
                "gt_rmsf_A",
            ]
            + [f"{s}_rmsf_A" for s in labels]
            + [f"{s}_absolute_difference_A" for s in labels]
        )
        for i, (chain, seqid, name, _) in enumerate(residues):
            writer.writerow(
                [i + 1, chain, seqid, name, counts[i], gt[i]]
                + [predictions[s][i] for s in labels]
                + [differences[s][i] for s in labels]
            )
    print(
        f"Plotted {len(labels)} pred ensembles and GT across {len(residues)} residues "
        f"({counts.sum()} atoms).\nOutput: {args.output_dir.resolve()}"
    )


if __name__ == "__main__":
    main()
