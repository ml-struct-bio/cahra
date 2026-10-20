#!/usr/bin/env python3
"""Reorder PDB atoms using a deposited PDB, preserving original record bytes."""

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
import re

import gemmi

ATOM_RECORDS = (b"ATOM  ", b"HETATM")


@dataclass
class Residue:
    key: tuple
    slots: list = field(default_factory=list)
    lines: list = field(default_factory=list)

    @property
    def name(self):
        return self.key[3].decode("ascii").strip()

    @property
    def label(self):
        return f"{
            self.name
        } chain {
            self.key[1].decode().strip() or '(blank)'
        } residue {self.key[2].decode().strip()}"

    @property
    def names(self):
        return [line[12:16].strip().decode("ascii") for line in self.lines]


def read_pdb(path):
    lines = path.read_bytes().splitlines(keepends=True)
    models = []
    chains = None
    current = None
    for index, line in enumerate(lines):
        record = line[:6]
        if record in (b"MODEL ", b"ENDMDL"):
            chains = None
            current = None
        elif record == b"TER   ":
            current = None
        elif record in ATOM_RECORDS:
            if len(line.rstrip(b"\r\n")) < 54:
                raise ValueError(f"{path}:{index + 1}: truncated atom record")
            if chains is None:
                chains = {}
                models.append(chains)
            key = (record, line[21:22], line[22:27], line[17:20])
            if current is None or current.key != key:
                current = Residue(key)
                chains.setdefault(key[1], []).append(current)
            current.slots.append(index)
            current.lines.append(line)
    if not models:
        raise ValueError(f"{path}: no atoms found")
    return lines, models


def align(query, target):
    scoring = gemmi.AlignmentScoring()
    scoring.match, scoring.mismatch = 3, -1
    scoring.gapo, scoring.gape = -5, -1
    result = gemmi.align_string_sequences(
        [r.name for r in query], [r.name for r in target], [], scoring
    )
    pairs, missing = {}, []
    q = t = 0
    for count, op in re.findall(r"(\d+)([MID])", result.cigar_str()):
        count = int(count)
        if op == "M":
            for offset in range(count):
                pairs[q + offset] = target[t + offset]
            q += count
            t += count
        elif op == "I":
            q += count
        else:
            missing.extend(target[t : t + count])
            t += count
    if q != len(query) or t != len(target):
        raise ValueError("Unexpected Gemmi alignment CIGAR")
    return pairs, missing


def templates_from(chains):
    variants = defaultdict(Counter)
    for residues in chains.values():
        for residue in residues:
            # Alternate conformations do not change the atom-name order.
            variants[residue.name][tuple(dict.fromkeys(residue.names))] += 1
    return variants


def choose_template(variants, residue):
    candidates = variants.get(residue.name)
    if not candidates:
        return None
    names = set(residue.names)
    # Prefer matching atom content, then the most frequent deposited variant.
    return min(
        candidates,
        key=lambda order: (
            len(names.symmetric_difference(order)),
            -candidates[order],
            order,
        ),
    )


def reorder(reference, source, destination, chain_map=None, report=print):
    _, ref_models = read_pdb(reference)
    if len(ref_models) != 1:
        raise ValueError("Reference must contain exactly one model")
    ref_chains = ref_models[0]
    variants = templates_from(ref_chains)
    lines, models = read_pdb(source)
    # These records contain atom serial references or accompany individual atoms.
    # They require a separate policy to satisfy strict non-atom byte preservation.
    unsupported = {
        line[:6].decode().strip()
        for line in lines
        if line[:6] in (b"ANISOU", b"SIGATM", b"SIGUIJ", b"CONECT")
    }
    if unsupported:
        raise ValueError(
            f"Unsupported serial-linked records: {', '.join(sorted(unsupported))}"
        )
    if destination.resolve() in (source.resolve(), reference.resolve()):
        raise ValueError("Output must not overwrite an input")
    output = list(lines)
    warnings = 0

    def warn(message):
        nonlocal warnings
        warnings += 1
        report(f"WARNING {source.name}: {message}")

    for model_number, chains in enumerate(models, 1):
        total = sum(
            len(r.lines) for residues in chains.values() for r in residues
        )
        if total > 99999:
            raise ValueError(
                "More than 99,999 atoms per model; decimal PDB serials would overflow"
            )
        assigned = {}
        for chain in chains:
            if chain_map and chain in chain_map:
                target = chain_map[chain]
                if target not in ref_chains:
                    raise ValueError(
                        f"Mapped reference chain {target!r} does not exist"
                    )
            elif len(chains) == len(ref_chains) == 1:
                target = next(iter(ref_chains))
            elif chain in ref_chains:
                target = chain
            else:
                raise ValueError(
                    f"Ambiguous chain correspondence for {chain!r}; specify --chain-map"
                )
            if target in assigned.values():
                raise ValueError(
                    "Multiple trajectory chains map to one reference chain; use"
                    " separate references"
                )
            assigned[chain] = target
        for chain in ref_chains.keys() - set(assigned.values()):
            for residue in ref_chains[chain]:
                warn(
                    f"model {model_number}: missing from trajectory:"
                    f" reference {residue.label}"
                )
        for chain, residues in chains.items():
            pairs, missing = align(residues, ref_chains[assigned[chain]])
            for residue in missing:
                warn(
                    f"model {model_number}: missing from trajectory:"
                    f" reference {residue.label}"
                )
            for index, residue in enumerate(residues):
                matched = pairs.get(index)
                if matched is None:
                    warn(
                        f"model {model_number}: {residue.label}"
                        f" has no aligned deposited residue"
                    )
                elif matched.name != residue.name:
                    warn(
                        f"model {model_number}: residue name differs:"
                        f" {residue.label} / reference {matched.label}"
                    )
                # Alignment is diagnostic only; the trajectory residue type
                # determines atom ordering, regardless of reference identity.
                order = choose_template(variants, residue)
                if order is None:
                    warn(
                        f"model {model_number}: {residue.label}: no deposited"
                        f" atom-order template; retaining atom order"
                    )
                    continue
                names = set(residue.names)
                extra = sorted(names - set(order))
                absent = sorted(set(order) - names)
                if extra or absent:
                    warn(
                        f"model {model_number}: {residue.label}: atom-name mismatch; "
                        f"trajectory-only={extra}, template-only={absent};"
                        f" unmatched atoms retained last"
                    )
                rank = {name: i for i, name in enumerate(order)}
                reordered = sorted(
                    residue.lines,
                    key=lambda line: rank.get(
                        line[12:16].strip().decode("ascii"), len(rank)
                    ),
                )
                for slot, line in zip(residue.slots, reordered):
                    output[slot] = line
        # Serial numbers follow final text order and restart for each MODEL.
        slots = sorted(
            slot
            for residues in chains.values()
            for r in residues
            for slot in r.slots
        )
        for serial, slot in enumerate(slots, 1):
            output[slot] = (
                output[slot][:6] + f"{serial:5d}".encode() + output[slot][11:]
            )

    verify_preservation(lines, output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("xb") as handle:
        handle.write(b"".join(output))
    report(
        f"WROTE {destination}: "
        f"{sum(line[:6] in ATOM_RECORDS for line in output)} atoms, "
        f"{len(models)} model(s), {warnings} warnings; "
        f"byte-preservation check passed"
    )


def verify_preservation(before, after):
    if len(before) != len(after):
        raise AssertionError("Line count changed")
    for old, new in zip(before, after):
        if old[:6] not in ATOM_RECORDS and old != new:
            raise AssertionError("Non-atom record changed")

    def payload(line):
        return line[:6] + line[11:]

    # Check separately per residue and model, so even identical residues cannot
    # exchange coordinates unnoticed.
    old_group, new_group, key = [], [], None

    def check():
        if Counter(old_group) != Counter(new_group):
            raise AssertionError("Atom data changed outside the serial field")

    for old, new in zip(before, after):
        next_key = (old[:6], old[17:27]) if old[:6] in ATOM_RECORDS else None
        if next_key != key or next_key is None:
            check()
            old_group, new_group, key = [], [], next_key
        if next_key is not None:
            old_group.append(payload(old))
            new_group.append(payload(new))
    check()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--chain-map",
        action="append",
        default=[],
        metavar="TRAJ:REF",
        help="Explicit chain correspondence; use '_' for a blank chain",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="Trajectory PDB files (shell globs supported)",
    )
    args = parser.parse_args()
    mapping = {}
    for item in args.chain_map:
        if not re.fullmatch(r".:.", item):
            parser.error(
                "--chain-map must be TRAJ:REF, each a single character"
            )
        mapping[item[0].replace("_", " ").encode()] = (
            item[2].replace("_", " ").encode()
        )
    destinations = [args.output_dir / p.name for p in args.inputs]
    if len(set(p.resolve() for p in destinations)) != len(destinations):
        parser.error("Input basenames must be unique")
    for destination in destinations:
        if destination.exists():
            parser.error(
                f"Refusing to overwrite existing output: {destination}"
            )
    for source, destination in zip(args.inputs, destinations):
        try:
            reorder(args.reference, source, destination, chain_map=mapping)
        except (ValueError, OSError) as exc:
            parser.exit(1, f"ERROR {source}: {exc}\n")


if __name__ == "__main__":
    main()
