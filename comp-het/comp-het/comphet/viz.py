"""The three evaluation figures: composition heatmap, per-submission small
multiples, and the F30/L30 exposure split."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

from .constants import JUNK_SPECIES, SPECIES_ORDER, TARGET_SPECIES

SPECIES_COLOR = {
    "PKM2": "#2a78d6",
    "GDH": "#eb6834",
    "ADH": "#1baf7a",
    "Aldolase": "#eda100",
    "Gamma-Globulin": "#52514e",
    "Beta-Amylase": "#898781",
}
# The two decoy species get neutral greys plus a hatch, so a junk-dominated
# class is recognisable at a glance and not just by reading the legend.
SPECIES_HATCH = {"Gamma-Globulin": "////", "Beta-Amylase": "\\\\\\\\"}

# Pale header-band tints behind the (bold, printed) submission name, so these
# only need to be mutually distinct. Assigned in fixed order and never cycled:
# plot_composition_heatmap raises rather than quietly reusing a tint.
SUBMISSION_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                     "#008300", "#4a3aa7", "#e34948", "#6f6f6f"]

# Exposure split: one hue, L30 as the lightened version of F30. Violet is not a
# species colour above, so it cannot be misread as PKM2/GDH/etc.
F30_COLOR = "#4a3aa7"

INK, INK_SECONDARY, INK_MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, SURFACE = "#e1e0d9", "#fcfcfb"


def _pale_tint(hex_color, alpha=0.16, base="#ffffff"):
    """Blend a hue toward white -- for a fill sitting behind dark text."""
    c = tuple(int(hex_color.lstrip("#")[i:i + 2], 16) for i in (0, 2, 4))
    b = tuple(int(base.lstrip("#")[i:i + 2], 16) for i in (0, 2, 4))
    mixed = tuple(round(alpha * ci + (1 - alpha) * bi) for ci, bi in zip(c, b))
    return "#{:02x}{:02x}{:02x}".format(*mixed)


def _grouped_rows(class_tbl: pd.DataFrame, order=None):
    """Shared row layout: classes grouped under a header row per submission,
    submissions in `order` (default alphabetical; any missing go last), classes
    within a submission by dominant species then purity. Returns
    (rows, colors) where rows are ("header", name) / ("data", class row)."""
    all_names = sorted(class_tbl["team"].unique())
    order = all_names if order is None else list(order) + [t for t in all_names if t not in order]
    pos = {t: i for i, t in enumerate(order)}

    dom_rank = {sp: i for i, sp in enumerate(SPECIES_ORDER)}
    tbl = class_tbl.copy()
    tbl["_pos"] = tbl["team"].map(pos)
    tbl["_dom_rank"] = tbl["dominant_species"].map(dom_rank)
    tbl = tbl.sort_values(["_pos", "_dom_rank", "purity_pct"],
                          ascending=[True, True, False]).reset_index(drop=True)

    # Tints are assigned alphabetically, independent of display order, so a
    # submission keeps its tint whether rows are sorted by rank or by name.
    if len(all_names) > len(SUBMISSION_COLORS):
        raise ValueError(
            f"{len(all_names)} submissions but only {len(SUBMISSION_COLORS)} SUBMISSION_COLORS -- "
            "add a distinct colour in viz.py rather than letting tints repeat")
    colors = dict(zip(all_names, SUBMISSION_COLORS))

    rows, cur = [], None
    for _, r in tbl.iterrows():
        if r["team"] != cur:
            rows.append(("header", r["team"]))
            cur = r["team"]
        rows.append(("data", r))
    return rows, colors


def _row_label(r) -> str:
    return f"{r['class_label']}  (n={int(r['n_particles']):,})"


def _no_matches_figure(title):
    """Explicitly mark an output figure when no composition can be calculated."""
    fig, ax = plt.subplots(figsize=(9.2, 3))
    ax.set_title(title, loc="left", color=INK)
    ax.text(0.5, 0.5, "No particles matched the reference.\nSee summary_qc.csv for details.",
            ha="center", va="center", transform=ax.transAxes, color=INK_SECONDARY)
    ax.set_axis_off()
    fig.tight_layout()
    return fig


def plot_composition_heatmap(class_tbl: pd.DataFrame, order=None):
    """Every submitted class (rows, grouped by submission) x every species
    (columns); cell = % of that class matching that species."""
    if class_tbl.empty:
        return _no_matches_figure("Submitted-class composition vs. ground-truth species")
    rows, colors = _grouped_rows(class_tbl, order)
    n = len(rows)

    mat = np.full((n, len(SPECIES_ORDER)), np.nan)
    row_labels, header_rows = [], []
    for i, (kind, payload) in enumerate(rows):
        if kind == "header":
            row_labels.append("")
            header_rows.append((i, payload))
        else:
            mat[i] = payload[SPECIES_ORDER].to_numpy()
            row_labels.append(_row_label(payload))

    fig, ax = plt.subplots(figsize=(9.2, max(6, 0.34 * n + 1.5)))
    cmap = LinearSegmentedColormap.from_list(
        "blue_seq", ["#fcfcfb", "#cde2fb", "#6da7ec", "#2a78d6", "#0d366b"])
    cmap.set_bad(SURFACE)
    im = ax.imshow(np.ma.masked_invalid(mat), cmap=cmap, vmin=0, vmax=100, aspect="auto")

    ncols = len(SPECIES_ORDER)
    for i, name in header_rows:
        ax.axhspan(i - 0.5, i + 0.5, color=_pale_tint(colors[name]), zorder=0)
        ax.text((ncols - 1) / 2.0, i, name, ha="center", va="center",
                fontsize=10.5, fontweight="bold", color=INK)

    ax.set_xticks(range(ncols))
    ax.set_xticklabels(SPECIES_ORDER, rotation=30, ha="right", fontsize=10)
    for tick, sp in zip(ax.get_xticklabels(), SPECIES_ORDER):
        if sp in JUNK_SPECIES:
            tick.set_color(INK_SECONDARY)
            tick.set_fontstyle("italic")

    ax.set_yticks(range(n))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.tick_params(axis="y", left=False)

    for i, (kind, _) in enumerate(rows):
        if kind == "header":
            continue
        for j in range(ncols):
            v = mat[i, j]
            if v >= 1.0:
                ax.text(j, i, f"{v:.0f}%", ha="center", va="center", fontsize=7.5,
                        color="white" if v >= 55 else INK)

    ax.axvline(len(TARGET_SPECIES) - 0.5, color=INK, linewidth=1.2, linestyle="--", alpha=0.6)
    ax.set_xticks(np.arange(-0.5, ncols, 1), minor=True)
    ax.grid(which="minor", axis="x", color=GRID, linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("% of class particles", fontsize=9, color=INK_SECONDARY)
    ax.set_title("Submitted-class composition vs. ground-truth species\n"
                 "(dashed line separates the 4 target complexes from the 2 junk/decoy species)",
                 fontsize=11, loc="left", color=INK)
    fig.tight_layout()
    return fig


def plot_submission_small_multiples(class_tbl: pd.DataFrame):
    """One panel per submission, one 100%-stacked bar per class."""
    if class_tbl.empty:
        return _no_matches_figure("Per-submission class composition")
    names = sorted(class_tbl["team"].unique())
    ncols = min(4, len(names))
    nrows = int(np.ceil(len(names) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 2.6 * nrows), squeeze=False)
    rank = {sp: i for i, sp in enumerate(SPECIES_ORDER)}

    for k, name in enumerate(names):
        ax = axes[k // ncols][k % ncols]
        sub = class_tbl[class_tbl["team"] == name].sort_values(
            "dominant_species", key=lambda s: s.map(rank))
        y = np.arange(len(sub))
        left = np.zeros(len(sub))
        for sp in SPECIES_ORDER:
            vals = sub[sp].to_numpy()
            ax.barh(y, vals, left=left, color=SPECIES_COLOR[sp], hatch=SPECIES_HATCH.get(sp),
                    edgecolor=SURFACE, linewidth=0.6, height=0.7)
            left += vals
        ax.set_yticks(y)
        ax.set_yticklabels([f"{cl} (n={int(c):,})" for cl, c in zip(sub["class_label"], sub["n_particles"])],
                           fontsize=7.5)
        ax.set_xlim(0, 100)
        ax.set_title(name, fontsize=10, fontweight="bold", color=INK, loc="left")
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(left=False)
        ax.set_xlabel("% of class", fontsize=8, color=INK_SECONDARY)
        ax.grid(axis="x", color=GRID, linewidth=0.5)

    for j in range(len(names), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    handles = [mpatches.Patch(facecolor=SPECIES_COLOR[sp], hatch=SPECIES_HATCH.get(sp),
                              edgecolor=SURFACE, label=sp) for sp in SPECIES_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Per-submission class composition (100% stacked)", fontsize=13, y=1.02)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    return fig


def plot_exposure_split(comp: pd.DataFrame, class_tbl: pd.DataFrame,
                        totals_by_exposure: pd.DataFrame, order=None):
    """F30 / L30 split of every submitted class, laid out like the heatmap.

    Exposure is largely fixed by species in this dataset (Gamma-Globulin and Beta-Amylase
    are F30-only, GDH/ADH/Aldolase ~80/20, PKM2 50/50), so a class's raw split
    mostly reflects its species mix. Each bar therefore also carries a marker at
    the split that mix *predicts* -- a boundary left of the marker means the
    class holds more L30 than its species content accounts for, right means
    less. That gap is the actual exposure-bias signal.
    """
    if class_tbl.empty:
        return _no_matches_figure("F30 / L30 split per submitted class")
    rows, colors = _grouped_rows(class_tbl, order)
    n = len(rows)
    l30_color = _pale_tint(F30_COLOR, alpha=0.5)

    by_exposure = comp.pivot_table(index=["team", "class_label"], columns="exposure",
                                   values="count", aggfunc="sum", fill_value=0)
    by_species = comp.pivot_table(index=["team", "class_label"], columns="species",
                                  values="count", aggfunc="sum", fill_value=0)
    baseline_l30 = totals_by_exposure.get("L30", 0) / totals_by_exposure.sum(axis=1)

    fig, ax = plt.subplots(figsize=(9.2, max(6, 0.34 * n + 1.5)))
    row_labels = []
    for i, (kind, payload) in enumerate(rows):
        if kind == "header":
            row_labels.append("")
            ax.axhspan(i - 0.5, i + 0.5, color=_pale_tint(colors[payload]), zorder=0)
            ax.text(50, i, payload, ha="center", va="center",
                    fontsize=10.5, fontweight="bold", color=INK)
            continue

        r = payload
        key = (r["team"], r["class_label"])
        row_labels.append(_row_label(r))
        counts = by_exposure.loc[key]
        total = counts.sum()
        f30_pct = 100 * counts.get("F30", 0) / total
        l30_pct = 100 - f30_pct

        sp_counts = by_species.loc[key]
        expected_l30 = 100 * float((sp_counts * baseline_l30.reindex(sp_counts.index).fillna(0)).sum() / total)
        expected_f30 = 100 - expected_l30

        ax.barh(i, f30_pct, color=F30_COLOR, height=0.7, zorder=2)
        ax.barh(i, l30_pct, left=f30_pct, color=l30_color, height=0.7, zorder=2)
        ax.plot([expected_f30, expected_f30], [i - 0.42, i + 0.42],
                color=INK, linewidth=1.6, solid_capstyle="butt", zorder=3)
        ax.text(f30_pct / 2, i, f"{f30_pct:.0f}%", ha="center", va="center",
                fontsize=7.5, color=SURFACE, zorder=4)
        ax.text(101.5, i, f"L30 {l30_pct:.0f}%", ha="left", va="center",
                fontsize=7.5, color=INK_SECONDARY, zorder=4)

    ax.set_yticks(range(n))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.tick_params(axis="y", left=False)
    ax.set_ylim(n - 0.5, -0.5)  # first row at the top, like the heatmap
    ax.set_xlim(0, 115)
    ax.set_xticks(range(0, 101, 20))
    ax.set_xlabel("% of class particles", fontsize=9, color=INK_SECONDARY)
    ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=1)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    handles = [
        mpatches.Patch(facecolor=F30_COLOR, label="F30"),
        mpatches.Patch(facecolor=l30_color, label="L30"),
        Line2D([], [], color=INK, linewidth=1.6, label="split predicted by the class's species mix"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=8.5, bbox_to_anchor=(0.5, 0.0))
    ax.set_title("F30 / L30 split per submitted class\n"
                 "(marker = split expected from species content alone; left of it means excess L30, right means excess F30)",
                 fontsize=11, loc="left", color=INK, pad=12)
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    return fig
