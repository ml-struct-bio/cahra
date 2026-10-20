#!/usr/bin/python3

import tempfile
from pathlib import Path
import unittest

from reorder_pdb import read_pdb, reorder


def atom(serial, name, residue="ALA", number=1, chain="A", alt=" "):
    return (
        f"ATOM  {serial:5d} {name:>4s}{alt}{residue:3s} {chain}{number:4d}    "
        f"{serial + .123:8.3f}{2.0:8.3f}{3.0:8.3f}  0.75 12.34      SYST C  \r\n"
    ).encode()


class ReorderTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.reference = self.root / "reference.pdb"
        self.source = self.root / "trajectory.pdb"
        self.output = self.root / "out.pdb"
        self.messages = []

    def run_reorder(self, **kwargs):
        reorder(
            self.reference,
            self.source,
            self.output,
            report=self.messages.append,
            **kwargs,
        )
        return self.output.read_bytes().splitlines(keepends=True)

    def test_columns_alternates_and_different_identifiers(self):
        self.reference.write_bytes(
            b"".join(
                atom(i, n, number=31)
                for i, n in enumerate(["N", "CA", "C", "O", "CB"], 1)
            )
        )
        original = [b"TITLE     preserve this exactly\r\n"]
        for i, (n, alt) in enumerate(
            [
                ("CB", "B"),
                ("CB", "A"),
                ("O", " "),
                ("N", " "),
                ("C", " "),
                ("CA", " "),
            ],
            11,
        ):
            original.append(atom(i, n, number=8, chain="Z", alt=alt))
        original.extend([b"TER      17      ALA Z   8\r\n", b"END"])
        self.source.write_bytes(b"".join(original))
        result = self.run_reorder()
        self.assertEqual(
            [line[12:16].strip() for line in result[1:-2]],
            [b"N", b"CA", b"C", b"O", b"CB", b"CB"],
        )
        self.assertEqual(
            [int(line[6:11]) for line in result[1:-2]], list(range(1, 7))
        )
        for old, new in zip(
            [original[i] for i in [4, 6, 5, 3, 1, 2]], result[1:-2]
        ):
            self.assertEqual(old[:6] + old[11:], new[:6] + new[11:])
        self.assertEqual(result[0], original[0])
        self.assertEqual(result[-2:], original[-2:])
        self.assertEqual(self.source.read_bytes(), b"".join(original))

    def test_insertions_deletions_renaming_and_unknown_atoms(self):
        names = [
            "ALA",
            "CYS",
            "GLU",
            "ARG",
            "TRP",
            "LEU",
            "GLY",
            "SER",
            "THR",
            "TYR",
        ]
        self.reference.write_bytes(
            b"".join(
                atom(i * 10 + j, n, residue, i + 31)
                for i, residue in enumerate(names)
                for j, n in enumerate(["N", "CA", "C", "O"])
            )
        )
        trajectory_names = ["MET"] + names[:4] + names[5:-1] + ["GLN"]
        self.source.write_bytes(
            b"".join(
                atom(i * 10 + j, n, residue, i + 1, "Z")
                for i, residue in enumerate(trajectory_names)
                for j, n in enumerate(["XX", "O", "C", "CA", "N"])
            )
        )
        self.run_reorder()
        log = "\n".join(self.messages)
        self.assertIn("no aligned deposited residue", log)
        self.assertIn("missing from trajectory: reference TRP", log)
        self.assertIn("residue name differs: GLN", log)
        self.assertIn("trajectory-only=['XX']", log)
        self.assertIn("no deposited atom-order template", log)
        _, models = read_pdb(self.output)
        self.assertEqual(models[0][b"Z"][1].names, ["N", "CA", "C", "O", "XX"])

    def test_multiple_models_reset_serials(self):
        self.reference.write_bytes(atom(1, "N") + atom(2, "CA"))
        self.source.write_bytes(
            b"".join(
                b"MODEL        "
                + str(i).encode()
                + b"\n"
                + atom(8, "CA")
                + atom(9, "N")
                + b"ENDMDL\n"
                for i in (1, 2)
            )
        )
        result = self.run_reorder()
        self.assertEqual(
            [int(line[6:11]) for line in result if line[:6] == b"ATOM  "],
            [1, 2, 1, 2],
        )

    def test_refuses_existing_output_and_linked_records(self):
        self.reference.write_bytes(atom(1, "N"))
        self.source.write_bytes(atom(5, "N"))
        self.output.write_bytes(b"keep me")
        with self.assertRaises(FileExistsError):
            self.run_reorder()
        self.assertEqual(self.output.read_bytes(), b"keep me")
        self.source.write_bytes(atom(5, "N") + b"CONECT    5    6\n")
        with self.assertRaisesRegex(ValueError, "serial-linked"):
            self.run_reorder()


if __name__ == "__main__":
    unittest.main()
