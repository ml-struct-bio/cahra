"""Derive the compact reference file (source_labels.pkl.gz) from the full source.star.

The repo ships the compact file because the full source.star is ~467 MB -- far
past what GitHub accepts. The pickle keeps only what scoring needs (the two
join keys plus species and exposure labels) and is ~13 MB. Nothing is lost for
evaluation: the dropped columns are coordinates, CTF extras and provenance
paths that the matching never reads.

It is an ordinary pandas DataFrame, so you can look at it directly:

    import pandas as pd
    pd.read_pickle("source_labels.pkl.gz")   # pandas handles the .gz itself

You only need this script to regenerate it, e.g. after downloading the full
source.star from EMPIAR-13973:

    python -m comphet.make_reference /path/to/source.star source_labels.pkl.gz
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import starfile

from .matching import normalize_star_key

NEEDED = ["rlnImageName", "rlnDefocusU", "rlnDefocusV", "source_species", "source_exposure"]


def build(source_star: Path, out_path: Path) -> None:
    ref = starfile.read(str(source_star))
    ref = ref["particles"] if isinstance(ref, dict) else ref
    missing = [c for c in NEEDED if c not in ref.columns]
    if missing:
        raise SystemExit(f"{source_star} is missing required column(s): {missing}")

    idx_stack = ref["rlnImageName"].map(normalize_star_key)
    out = pd.DataFrame({
        "primary_key": np.array([i * 1_000_000 + s for i, s in idx_stack], dtype=np.int64),
        "defocus_key": (np.round(ref["rlnDefocusU"].to_numpy(), 3) * 1_000_000.0
                        + np.round(ref["rlnDefocusV"].to_numpy(), 3)),
        "species": pd.Categorical(ref["source_species"]),
        "exposure": pd.Categorical(ref["source_exposure"]),
    }).sort_values("primary_key").reset_index(drop=True)

    out.to_pickle(out_path)  # .gz extension -> pandas compresses automatically
    print(f"wrote {out_path} ({Path(out_path).stat().st_size / 1e6:.1f} MB) "
          f"from {len(out):,} particles")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("source_star", type=Path, help="full source.star (e.g. from EMPIAR-13973)")
    p.add_argument("out_path", type=Path, nargs="?", default=Path("source_labels.pkl.gz"))
    a = p.parse_args()
    build(a.source_star, a.out_path)


if __name__ == "__main__":
    main()
