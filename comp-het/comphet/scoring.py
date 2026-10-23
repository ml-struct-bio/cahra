"""Purity-weighted scoring: 25 points per target species, scaled by precision.

    target_score = 25 * (particles of that species in the class)
                      / (total particles in that class)

"The class assigned to that target" means whichever of a submission's classes
have that species as their dominant (plurality) identity -- no manual "this
class is my PKM2 class" declaration is needed. If a target was split across
several classes, the single highest-purity one is scored: pooling them would
let a large, dirty class drag down a genuinely clean one. A target that is
never any class's dominant species scores 0.

This scores purity, not completeness: a class holding 10% of a target's
particles perfectly cleanly earns the full 25 points. Junk contamination is
penalised automatically -- decoy particles inside a target class lower that
class's purity -- so no separate junk term is needed.
"""
import pandas as pd

from .constants import SPECIES_ORDER, TARGET_SPECIES


def class_level_table(comp: pd.DataFrame) -> pd.DataFrame:
    """One row per (team, class_label): % composition per species, particle count,
    and the class's dominant species with that species' share (its purity).

    The dominant species is the plurality across the WHOLE class. `comp` holds
    one row per (species, exposure), so the species shares must be summed over
    exposure first -- taking the largest single row instead would let a species
    that sits in one exposure beat a larger species split across both. PKM2 is
    split 50/50, so that mistake both mislabels PKM2-dominant classes and halves
    their apparent purity.
    """
    if comp.empty:
        return pd.DataFrame(columns=["team", "class_label", *SPECIES_ORDER,
                                     "n_particles", "dominant_species", "purity_pct"])

    pivot = comp.pivot_table(index=["team", "class_label"], columns="species",
                             values="pct_of_class", aggfunc="sum", fill_value=0.0)
    for sp in SPECIES_ORDER:
        if sp not in pivot.columns:
            pivot[sp] = 0.0
    pivot = pivot[SPECIES_ORDER]

    counts = comp.groupby(["team", "class_label"])["count"].sum().rename("n_particles")
    dominant = pivot.idxmax(axis=1).rename("dominant_species")
    purity = pivot.max(axis=1).rename("purity_pct")

    return pivot.join(counts).join(dominant).join(purity).reset_index()


def compute_target_scores(comp: pd.DataFrame, class_tbl: pd.DataFrame,
                          target_species=TARGET_SPECIES) -> pd.DataFrame:
    """Per submission: purity, points (out of 25) per target, and the total /100.
    Sorted best-first."""
    rows = []
    for team in sorted(class_tbl["team"].unique()):
        team_classes = class_tbl[class_tbl["team"] == team]
        row = {"team": team}
        total = 0.0
        for sp in target_species:
            candidates = team_classes[team_classes["dominant_species"] == sp]
            if len(candidates) == 0:
                purity, chosen = 0.0, "-"
            else:
                # If a target was split across several classes, score the single
                # cleanest one. Pooling them would let a large, dirty class drag
                # down a genuinely clean one and understate what was achieved.
                best = candidates.loc[candidates["purity_pct"].idxmax()]
                purity = best["purity_pct"] / 100.0
                chosen = best["class_label"]
            score = 25 * purity
            total += score
            row[f"{sp}_class"] = chosen
            row[f"{sp}_purity_pct"] = round(purity * 100, 1)
            row[f"{sp}_score"] = round(score, 2)
        row["total_score"] = round(total, 2)
        rows.append(row)
    columns = ["team"]
    for sp in target_species:
        columns.extend([f"{sp}_class", f"{sp}_purity_pct", f"{sp}_score"])
    columns.append("total_score")
    return (pd.DataFrame(rows, columns=columns).sort_values("total_score", ascending=False)
            .reset_index(drop=True))


def scores_display_table(scores_detail: pd.DataFrame, target_species=TARGET_SPECIES) -> pd.DataFrame:
    """Compact view: rank, submission, purity% + points per target, total /100."""
    display_cols = ["rank", "team"]
    for sp in target_species:
        display_cols += [f"{sp}_purity_pct", f"{sp}_score"]
    display_cols += ["total_score"]

    out = scores_detail[display_cols].rename(columns={
        **{f"{sp}_purity_pct": f"{sp} %" for sp in target_species},
        **{f"{sp}_score": f"{sp} pts" for sp in target_species},
        "team": "submission", "total_score": "Total /100",
    })
    round_to = {f"{sp} %": 1 for sp in target_species}
    round_to.update({f"{sp} pts": 1 for sp in target_species})
    round_to["Total /100"] = 1
    return out.round(round_to)
