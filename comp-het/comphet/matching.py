"""Match submitted particle sets (.cs / .star) back to the ground-truth labels.

Every submission references particle images by a per-stack index plus a stack
filename. Directory names and cryoSPARC UID prefixes differ between pipelines,
but the (index-within-stack, stack number) pair survives re-import, so that
pair is the primary join key back to source.star.

If identifiers are missing, invalid or concentrated in suspiciously few stacks,
CTF defocus can provide diagnostic match counts. Those matches never contribute
to composition or scores; the caller must report the metadata errors and stop
evaluation.
"""
import gzip
import re
import shutil
import tempfile
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import starfile

from .constants import SPECIES_RENAME

_NUM_RE = re.compile(r"0*(\d+)")

# Large files referencing very few stacks require a diagnostic and review of
# the stack metadata. Defocus matches from this check must never be scored.
DEGENERATE_ROW_THRESHOLD = 1000
DEGENERATE_STACK_THRESHOLD = 3

PARTICLE_SUFFIXES = (".cs", ".star")


def _stack_num(basename: str) -> int:
    """Numeric part of a stack filename, e.g. '0022.mrcs' -> 22."""
    stem = basename.split(".")[0]
    match = _NUM_RE.match(stem)
    if match is None:
        raise ValueError(f"Stack filename has no numeric ID: {basename!r}")
    return int(match.group(1))


def normalize_star_key(rln_image_name: str) -> Tuple[int, int]:
    """'0045@Data/0022.mrcs' -> (index, stack_number)."""
    if not isinstance(rln_image_name, str):
        raise ValueError(f"Invalid rlnImageName: {rln_image_name!r}")
    idx, path = rln_image_name.split("@")
    idx = int(idx)
    if idx < 1:
        raise ValueError("rlnImageName must use a positive, one-based image index")
    basename = path.rsplit("/", 1)[-1]
    return idx, _stack_num(basename)


def normalize_cs_keys(paths: np.ndarray, idxs: np.ndarray):
    """cryoSPARC blob/path + blob/idx (0-based) -> (combined key, stack number)."""
    keys = np.empty(len(paths), dtype=np.int64)
    stacks = np.empty(len(paths), dtype=np.int64)
    for i, (p, ix) in enumerate(zip(paths, idxs)):
        if isinstance(p, bytes):
            p = p.decode()
        if not isinstance(p, str):
            raise ValueError(f"Invalid blob/path: {p!r}")
        idx = int(ix)
        if idx < 0 or idx != ix:
            raise ValueError(f"blob/idx must be a nonnegative integer: {ix!r}")
        basename = p.rsplit("/", 1)[-1]
        # Strip cryoSPARC import UID prefix(es). A stack re-imported through a
        # second project picks up a second UID (<uid>_<uid>_0023.mrcs), so take
        # the stack number from the LAST token rather than splitting off one.
        basename = basename.rsplit("_", 1)[-1]
        s = _stack_num(basename)
        keys[i] = (idx + 1) * 1_000_000 + s
        stacks[i] = s
    return keys, stacks


def read_submission(path: Path) -> dict:
    """Read identifiers and optional defocus, retaining errors for diagnostics.

    Missing or invalid identifiers leave keys/stacks unset. Defocus is kept
    unconverted here so missing or malformed CTF values cannot break otherwise
    valid identifier matching.
    """
    path = Path(path)
    if path.suffix == ".cs":
        data = np.load(str(path), allow_pickle=False)
        columns = data.dtype.names or ()
        required = ("blob/path", "blob/idx")
        defocus_columns = ("ctf/df1_A", "ctf/df2_A")
    elif path.suffix == ".star":
        d = starfile.read(str(path))
        data = d["particles"] if isinstance(d, dict) else d
        columns = data.columns
        required = ("rlnImageName",)
        defocus_columns = ("rlnDefocusU", "rlnDefocusV")
    else:
        raise ValueError(f"Unsupported file type: {path} (expected .cs or .star)")

    keys = stacks = None
    primary_error = ""
    missing = [name for name in required if name not in columns]
    if missing:
        primary_error = f"Missing required column(s): {', '.join(missing)}"
    else:
        try:
            if path.suffix == ".cs":
                keys, stacks = normalize_cs_keys(data["blob/path"], data["blob/idx"])
            else:
                idx_stack = data["rlnImageName"].map(normalize_star_key)
                keys = np.array([i * 1_000_000 + s for i, s in idx_stack], dtype=np.int64)
                stacks = np.array([s for i, s in idx_stack], dtype=np.int64)
        except (ValueError, TypeError, OverflowError) as exc:
            keys = stacks = None
            primary_error = f"Invalid particle identifier: {exc}"

    defocus_u, defocus_v = (
        np.asarray(data[name]) if name in columns else None
        for name in defocus_columns
    )
    return dict(n_rows_raw=len(data), keys=keys, stacks=stacks,
                primary_error=primary_error, defocus_u=defocus_u, defocus_v=defocus_v)


def class_label_from_path(path: Path) -> str:
    """Use the full filename without its extension as the class label."""
    return path.stem


def find_submissions(results_dir: Path) -> List[Tuple[str, str, Path]]:
    """Every (submission, class_label, path) under results_dir.

    Each immediate subdirectory holding .cs/.star files is one submission, named
    after the directory. Directories with no particle files (the package itself,
    an output folder, .git) are skipped. Class labels must be unique within each
    submission; collisions are rejected before any particle files are read.
    """
    results_dir = Path(results_dir)
    out = []
    for sub_dir in sorted(p for p in results_dir.iterdir() if p.is_dir() and not p.name.startswith(".")):
        files = sorted(f for f in sub_dir.iterdir()
                       if f.is_file() and f.suffix in PARTICLE_SUFFIXES)
        labels = {}
        for f in files:
            label = class_label_from_path(f)
            if label in labels:
                raise SystemExit(
                    f"Duplicate class label '{label}' in submission '{sub_dir.name}':\n"
                    f"  {labels[label].name}\n"
                    f"  {f.name}\n"
                    "Use a distinct class label for each class and keep only one file per class. "
                    "Files are not combined.")
            labels[label] = f
            out.append((sub_dir.name, label, f))
    return sorted(out, key=lambda t: (t[0], t[1]))


def _reference_from_arrays(primary, defocus, species, exposure) -> dict:
    """Package the ground truth into sorted arrays for vectorised lookup."""
    for old, new in SPECIES_RENAME.items():
        species = np.where(species == old, new, species)
    p_order = np.argsort(primary)
    d_order = np.argsort(defocus)
    counts = pd.Series(species).value_counts()
    totals_by_exposure = (pd.DataFrame({"species": species, "exposure": exposure})
                          .groupby(["species", "exposure"]).size().unstack(fill_value=0))
    return dict(
        primary_sorted=primary[p_order], primary_species=species[p_order],
        primary_exposure=exposure[p_order],
        defocus_sorted=defocus[d_order], defocus_species=species[d_order],
        defocus_exposure=exposure[d_order],
        totals=counts, totals_by_exposure=totals_by_exposure,
    )


def load_reference(ref_path: Path) -> dict:
    """Load the ground truth as sorted lookup arrays.

    Accepts the shipped compact source_labels.pkl.gz, or a full source.star /
    source.star.gz if you would rather work from the original (see
    comphet.make_reference and EMPIAR-13973).
    """
    ref_path = Path(ref_path)
    if not ref_path.exists():
        raise SystemExit(f"Reference file not found: {ref_path}")

    if ref_path.name.endswith((".pkl", ".pkl.gz", ".pickle")):
        df = pd.read_pickle(ref_path)  # pandas handles the .gz itself
        return _reference_from_arrays(
            df["primary_key"].to_numpy(), df["defocus_key"].to_numpy(),
            df["species"].to_numpy().astype("U"), df["exposure"].to_numpy().astype("U"))

    if ref_path.suffix == ".gz":
        with tempfile.NamedTemporaryFile(suffix=".star", delete=False) as tmp:
            with gzip.open(ref_path, "rb") as fi:
                shutil.copyfileobj(fi, tmp)
            tmp_path = tmp.name
        try:
            ref = starfile.read(tmp_path)["particles"].reset_index(drop=True)
        finally:
            Path(tmp_path).unlink(missing_ok=True)
    else:
        ref = starfile.read(str(ref_path))["particles"].reset_index(drop=True)

    idx_stack = ref["rlnImageName"].map(normalize_star_key)
    primary = np.array([i * 1_000_000 + s for i, s in idx_stack], dtype=np.int64)
    defocus = (np.round(ref["rlnDefocusU"].to_numpy(), 3) * 1_000_000.0
               + np.round(ref["rlnDefocusV"].to_numpy(), 3))
    return _reference_from_arrays(primary, defocus,
                                  ref["source_species"].to_numpy().astype("U"),
                                  ref["source_exposure"].to_numpy().astype("U"))


def _lookup(query, sorted_keys, species, exposure, tol=None):
    """Vectorised exact lookup of `query` against `sorted_keys`.

    Returns (species, exposure) for the entries that matched. `tol` is for the
    float defocus key; the integer primary key compares exactly.
    """
    if len(query) == 0:
        return np.array([], dtype="U"), np.array([], dtype="U")
    pos = np.clip(np.searchsorted(sorted_keys, query), 0, len(sorted_keys) - 1)
    found = sorted_keys[pos]
    hit = np.isclose(found, query, rtol=0, atol=tol) if tol else found == query
    return species[pos][hit], exposure[pos][hit]


def _defocus_diagnostic(sub: dict, ref: dict) -> dict:
    """Count possible defocus matches without returning particle assignments."""
    result = dict(n_defocus_unique=None, n_defocus_matched=None,
                  n_defocus_unmatched=None, n_defocus_invalid_rows=None,
                  defocus_note="")
    if sub["defocus_u"] is None or sub["defocus_v"] is None:
        result["defocus_note"] = "Diagnostic unavailable: both defocus columns are required"
        return result

    u = np.asarray(pd.to_numeric(sub["defocus_u"], errors="coerce"), dtype=np.float64)
    v = np.asarray(pd.to_numeric(sub["defocus_v"], errors="coerce"), dtype=np.float64)
    valid = np.isfinite(u) & np.isfinite(v)
    dkeys = np.unique(np.round(u[valid], 3) * 1_000_000.0 + np.round(v[valid], 3))
    species, _ = _lookup(dkeys, ref["defocus_sorted"], ref["defocus_species"],
                         ref["defocus_exposure"], tol=1e-6)
    result.update(n_defocus_unique=len(dkeys), n_defocus_matched=len(species),
                  n_defocus_unmatched=len(dkeys) - len(species),
                  n_defocus_invalid_rows=int((~valid).sum()),
                  defocus_note="Diagnostic only; defocus matches are not used for evaluation")
    return result


def build_summary(results_dir: Path, ref: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """One pass over every submitted class file under results_dir.

    Returns (qc, comp):
      qc   -- one row per file: primary matches, metadata errors, diagnostic counts.
      comp -- composition from primary identifiers only; never from defocus.

    The caller must stop evaluation if any primary_error is nonempty. Duplicate
    and unmatched identifiers produce QC warnings and do not stop evaluation.

    The `team` column names the submission folder; it is called `team` because
    the challenge was run between teams, but it is just the folder name.
    """
    qc_rows, comp_rows = [], []
    keys_by_team = {}

    for team, class_label, f in find_submissions(results_dir):
        sub = read_submission(f)
        n_raw = sub["n_rows_raw"]
        row = dict(team=team, class_label=class_label, file=f.name,
                   n_rows_raw=n_raw, n_unique_particles=None,
                   n_stacks_referenced=None, match_method="path_idx",
                   n_matched=0, n_unmatched=None, pct_duplicated_rows=None,
                   n_duplicated_rows=None, n_shared_particles=None, warning="",
                   primary_error=sub["primary_error"],
                   n_defocus_unique=None, n_defocus_matched=None,
                   n_defocus_unmatched=None, n_defocus_invalid_rows=None,
                   defocus_note="")

        if not row["primary_error"]:
            uniq_keys = np.unique(sub["keys"])
            n_unique = len(uniq_keys)
            n_stacks = len(np.unique(sub["stacks"]))
            row.update(n_unique_particles=n_unique, n_stacks_referenced=n_stacks,
                       n_duplicated_rows=n_raw - n_unique,
                       pct_duplicated_rows=100 * (1 - n_unique / n_raw) if n_raw else 0.0)
            if n_raw >= DEGENERATE_ROW_THRESHOLD and n_stacks <= DEGENERATE_STACK_THRESHOLD:
                row["primary_error"] = (
                    f"Suspected damaged stack metadata: {n_raw:,} rows reference "
                    f"only {n_stacks} stack(s)")

        if row["primary_error"]:
            row["match_method"] = "defocus_diagnostic"
            row.update(_defocus_diagnostic(sub, ref))
            qc_rows.append(row)
            continue

        species_arr, exposure_arr = _lookup(
            uniq_keys, ref["primary_sorted"], ref["primary_species"], ref["primary_exposure"])
        n_matched = len(species_arr)
        row.update(n_matched=n_matched, n_unmatched=n_unique - n_matched,
                   n_shared_particles=0)
        keys_by_team.setdefault(team, []).append((row, uniq_keys))
        qc_rows.append(row)

        if n_matched:
            sub_df = pd.DataFrame({"species": species_arr, "exposure": exposure_arr})
            for (sp, ex), cnt in sub_df.groupby(["species", "exposure"]).size().items():
                comp_rows.append(dict(team=team, class_label=class_label, file=f.name,
                                      species=sp, exposure=ex, count=cnt))

    # Count overlap between class files within a submission. Each array is
    # already unique, so repeated keys here mean a particle occurs in more
    # than one file. Classes remain independent for scoring.
    for files in keys_by_team.values():
        if len(files) < 2:
            continue
        keys, counts = np.unique(np.concatenate([keys for _, keys in files]), return_counts=True)
        shared = keys[counts > 1]
        for row, keys in files:
            row["n_shared_particles"] = len(np.intersect1d(keys, shared, assume_unique=True))

    for row in qc_rows:
        notes = []
        if row["n_duplicated_rows"]:
            notes.append(f"{row['n_duplicated_rows']:,} duplicate row(s); repeated identifiers "
                         "are counted once within this file")
        if row["n_shared_particles"]:
            notes.append(f"{row['n_shared_particles']:,} particle identifier(s) also occur in "
                         "other class files in this submission; each class is evaluated separately")
        if row["n_unmatched"]:
            notes.append(f"{row['n_unmatched']:,} of {row['n_unique_particles']:,} unique "
                         "particle identifiers did not match the reference; only matched "
                         "particles contribute to composition and scores")
        if not row["primary_error"] and row["n_matched"] == 0:
            notes.append("No matched particles in this file; no score can be calculated from it")
        row["warning"] = "; ".join(notes)

    qc = pd.DataFrame(qc_rows)
    comp = pd.DataFrame(comp_rows, columns=["team", "class_label", "file", "species", "exposure",
                                          "count", "pct_of_class", "recall_pct"])
    if comp.empty:
        return qc, comp

    class_totals = comp.groupby(["team", "class_label"])["count"].transform("sum")
    comp["pct_of_class"] = 100 * comp["count"] / class_totals
    comp["recall_pct"] = 100 * comp["count"] / comp["species"].map(ref["totals"])
    return qc, comp
