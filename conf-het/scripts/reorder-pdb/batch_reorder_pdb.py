#!/usr/bin/env python3
"""Mirror trajectory directories under an output root, with per-PDB logs."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import gemmi

from reorder_pdb import ATOM_RECORDS, reorder, verify_preservation


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument(
        "--mirror-base",
        required=True,
        type=Path,
        help="Paths relative to this directory are retained under output-root",
    )
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("directories", nargs="+", type=Path)
    args = parser.parse_args()
    base = args.mirror_base.resolve()
    jobs = []
    for directory in args.directories:
        files = sorted(
            p
            for p in directory.rglob("*")
            if p.is_file() and p.suffix.lower() == ".pdb"
        )
        if not files:
            parser.error(f"No PDBs found in {directory}")
        for source in files:
            relative = source.resolve().relative_to(base)
            destination = (
                args.output_root
                / relative.parent
                / (source.stem + "_reordered.pdb")
            )
            log = destination.with_suffix(".reorder.log")
            jobs.append((source, destination, log))
    all_outputs = [
        p.resolve() for _, destination, log in jobs for p in (destination, log)
    ]
    manifest = args.output_root / "reorder_batch_manifest.json"
    all_outputs.append(manifest.resolve())
    if len(set(all_outputs)) != len(all_outputs):
        parser.error("Duplicate output paths")
    inputs = {p.resolve() for p, _, _ in jobs} | {args.reference.resolve()}
    for path in all_outputs:
        if path in inputs or path.exists():
            parser.error(f"Refusing to overwrite: {path}")
    original_hashes = {p: digest(p) for p in inputs}
    results = []
    for index, (source, destination, log) in enumerate(jobs, 1):
        destination.parent.mkdir(parents=True, exist_ok=True)
        messages = []
        with log.open("x") as handle:

            def report(message):
                messages.append(message)
                print(message, file=handle, flush=True)

            reorder(args.reference, source, destination, report=report)
            before = source.read_bytes().splitlines(keepends=True)
            after = destination.read_bytes().splitlines(keepends=True)
            verify_preservation(before, after)
            atom_lines = [line for line in after if line[:6] in ATOM_RECORDS]
            # These trajectories have one model; verify serials from saved bytes.
            structure = gemmi.read_structure(str(destination))
            if len(structure) != 1:
                raise ValueError(
                    f"Batch validator expects one model: {source}"
                )
            if [int(line[6:11]) for line in atom_lines] != list(
                range(1, len(atom_lines) + 1)
            ):
                raise AssertionError(f"Nonconsecutive serials: {destination}")
            if sum(
                len(residue) for chain in structure[0] for residue in chain
            ) != len(atom_lines):
                raise AssertionError(
                    f"Gemmi atom count differs: {destination}"
                )
            if digest(source) != original_hashes[source.resolve()]:
                raise AssertionError(f"Input changed: {source}")
            print(
                "VALIDATED saved output: all non-serial atom bytes and"
                " non-atom records preserved; "
                "serials consecutive; Gemmi atom count matches; input"
                " SHA-256 unchanged",
                file=handle,
            )
        results.append(
            dict(
                source=str(source.resolve()),
                output=str(destination.resolve()),
                log=str(log.resolve()),
                atoms=len(atom_lines),
                warnings=sum(m.startswith("WARNING ") for m in messages),
                input_sha256=original_hashes[source.resolve()],
                output_sha256=digest(destination),
                validated=True,
            )
        )
        if index % 20 == 0 or index == len(jobs):
            print(
                f"Processed and validated {index}/{len(jobs)} PDBs", flush=True
            )
    for path, original in original_hashes.items():
        if digest(path) != original:
            raise AssertionError(f"Input changed during batch: {path}")
    with manifest.open("x") as handle:
        json.dump(
            dict(
                reference=str(args.reference.resolve()),
                reference_sha256=original_hashes[args.reference.resolve()],
                original_files_unchanged=True,
                files=results,
            ),
            handle,
            indent=2,
        )
        handle.write("\n")
    print(
        "Counts by output directory:",
        dict(Counter(str(Path(r["output"]).parent) for r in results)),
    )
    print(
        "Warning counts per file:",
        dict(Counter(r["warnings"] for r in results)),
    )
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
