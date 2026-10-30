"""Regression tests for primary matching and diagnostic-only defocus lookup."""
import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import starfile

from comphet.evaluate import run
from comphet.matching import build_summary, find_submissions, load_reference
import matplotlib.pyplot as plt


class MatchingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.addCleanup(plt.close, "all")
        self.root = Path(self.temp.name)
        self.inputs = self.root / "inputs"
        self.submission = self.inputs / "example"
        self.submission.mkdir(parents=True)
        self.output = self.root / "output"

        self.indices = np.arange(1000) // 4 + 1
        self.stacks = np.arange(1000) % 4 + 1
        self.u = np.arange(1000, dtype=float) + 1000
        self.v = self.u + 0.25
        labels = pd.DataFrame({
            "primary_key": self.indices * 1_000_000 + self.stacks,
            "defocus_key": self.u * 1_000_000.0 + self.v,
            "species": "PKM2",
            "exposure": "F30",
        })
        self.reference_path = self.root / "reference.pkl"
        labels.to_pickle(self.reference_path)
        self.ref = load_reference(self.reference_path)

    def star_data(self):
        return pd.DataFrame({
            "rlnImageName": [f"{i}@particles/{s:04d}.mrcs"
                             for i, s in zip(self.indices, self.stacks)],
            "rlnDefocusU": self.u,
            "rlnDefocusV": self.v,
        })

    def cs_data(self):
        return np.rec.fromarrays([
            np.array([f"particles/uid_{s:04d}.mrcs" for s in self.stacks], dtype="S32"),
            self.indices - 1, self.u, self.v,
        ], names=["blob/path", "blob/idx", "ctf/df1_A", "ctf/df2_A"])

    def write_star(self, data, name="class1.star"):
        starfile.write({"particles": data}, self.submission / name, overwrite=True)

    def write_cs(self, data, name="class1.cs"):
        with (self.submission / name).open("wb") as stream:
            np.save(stream, data, allow_pickle=False)

    def test_duplicate_labels_are_rejected_for_different_filenames(self):
        cases = [("runA_class3.star", "runA_class3.cs", "runA_class3"),
                 ("class3.star", "class3.cs", "class3"),
                 ("particles.star", "particles.cs", "particles")]
        for first, second, label in cases:
            with self.subTest(first=first, second=second), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                submission = root / "example"
                submission.mkdir()
                (submission / first).touch()
                (submission / second).touch()
                with self.assertRaisesRegex(SystemExit, f"Duplicate class label '{label}'") as error:
                    find_submissions(root)
                self.assertIn(first, str(error.exception))
                self.assertIn(second, str(error.exception))

    def test_duplicate_labels_stop_evaluation_before_loading_data(self):
        (self.submission / "runA_class3.star").touch()
        (self.submission / "runA_class3.cs").touch()
        with patch("comphet.evaluate.load_reference") as load:
            with self.assertRaisesRegex(SystemExit, "Duplicate class label"):
                run(self.inputs, self.reference_path, self.output)
            load.assert_not_called()
        self.assertEqual(list(self.output.iterdir()), [])

    def test_full_filename_labels_keep_classes_separate(self):
        data = self.star_data()
        self.write_star(data.iloc[:100], name="runA_class3.star")
        self.write_star(data.iloc[100:300], name="runB_class3.v2.star")
        qc, comp = build_summary(self.inputs, self.ref)
        self.assertEqual(qc["class_label"].tolist(), ["runA_class3", "runB_class3.v2"])
        self.assertEqual(comp.groupby("class_label")["count"].sum().to_dict(),
                         {"runA_class3": 100, "runB_class3.v2": 200})

    def test_same_filename_in_separate_submissions_is_allowed(self):
        self.write_star(self.star_data())
        other = self.inputs / "other_submission"
        other.mkdir()
        starfile.write({"particles": self.star_data()}, other / "class1.star")
        qc, comp = build_summary(self.inputs, self.ref)
        self.assertEqual(qc["class_label"].tolist(), ["class1", "class1"])
        self.assertEqual(qc["team"].tolist(), ["example", "other_submission"])
        self.assertEqual(comp.groupby("team")["count"].sum().to_dict(),
                         {"example": 1000, "other_submission": 1000})

    def test_valid_identifiers_above_stack_threshold_need_no_defocus(self):
        self.write_star(self.star_data().drop(columns=["rlnDefocusU", "rlnDefocusV"]))
        data = self.cs_data()
        self.write_cs(data[["blob/path", "blob/idx"]], name="class2.cs")
        qc, comp = build_summary(self.inputs, self.ref)
        self.assertEqual(qc["match_method"].tolist(), ["path_idx", "path_idx"])
        self.assertEqual(qc["n_matched"].tolist(), [1000, 1000])
        self.assertTrue(qc["n_defocus_matched"].isna().all())
        self.assertEqual(comp["count"].sum(), 2000)
        self.assertEqual(comp.groupby("class_label")["count"].sum().to_dict(),
                         {"class1": 1000, "class2": 1000})

    def test_bad_defocus_does_not_affect_valid_identifiers(self):
        data = self.star_data()
        data["rlnDefocusU"] = "invalid"
        data["rlnDefocusV"] = np.nan
        self.write_star(data)
        qc, comp = build_summary(self.inputs, self.ref)
        self.assertEqual(qc.iloc[0]["n_matched"], 1000)
        self.assertEqual(qc.iloc[0]["match_method"], "path_idx")
        self.assertEqual(comp["count"].sum(), 1000)

    def test_stack_diagnostic_thresholds_even_with_valid_unique_identifiers(self):
        cases = [(999, 1, False), (999, 3, False), (1000, 1, True),
                 (1000, 2, True), (1000, 3, True), (1000, 4, False)]
        for n_rows, n_stacks, diagnostic in cases:
            with self.subTest(rows=n_rows, stacks=n_stacks):
                indices = np.arange(1000) // n_stacks + 1
                stacks = np.arange(1000) % n_stacks + 1
                labels = pd.read_pickle(self.reference_path)
                labels["primary_key"] = indices * 1_000_000 + stacks
                labels.to_pickle(self.reference_path)
                ref = load_reference(self.reference_path)
                data = self.star_data()
                data["rlnImageName"] = [f"{i}@{s:04d}.mrcs" for i, s in zip(indices, stacks)]
                self.write_star(data.iloc[:n_rows])
                qc, comp = build_summary(self.inputs, ref)
                row = qc.iloc[0]
                self.assertEqual(row["n_unique_particles"], n_rows)
                self.assertEqual(row["n_stacks_referenced"], n_stacks)
                if diagnostic:
                    self.assertEqual(row["match_method"], "defocus_diagnostic")
                    self.assertIn("Suspected damaged stack metadata", row["primary_error"])
                    self.assertEqual(row["n_defocus_matched"], n_rows)
                    self.assertEqual(row["n_matched"], 0)
                    self.assertTrue(comp.empty)
                else:
                    self.assertEqual(row["match_method"], "path_idx")
                    self.assertEqual(row["n_matched"], n_rows)
                    self.assertEqual(comp["count"].sum(), n_rows)

    def test_cs_stack_diagnostic_logs_and_stops_with_or_without_defocus(self):
        for include_defocus in [True, False]:
            with self.subTest(include_defocus=include_defocus):
                data = self.cs_data()
                data["blob/path"] = b"particles/uid_0001.mrcs"
                if not include_defocus:
                    data = data[["blob/path", "blob/idx"]]
                self.write_cs(data)
                stdout = io.StringIO()
                with contextlib.redirect_stdout(stdout), self.assertRaisesRegex(SystemExit, "Evaluation stopped"):
                    run(self.inputs, self.reference_path, self.output)
                self.assertIn("1,000 rows reference only 1 stack(s)", stdout.getvalue())
                if include_defocus:
                    self.assertIn("1,000 of 1,000 unique keys matched", stdout.getvalue())
                else:
                    self.assertIn("Diagnostic unavailable", stdout.getvalue())
                self.assertEqual([p.name for p in self.output.iterdir()], ["summary_qc.csv"])
                qc = pd.read_csv(self.output / "summary_qc.csv")
                self.assertEqual(qc.iloc[0]["match_method"], "defocus_diagnostic")
                self.assertEqual(qc.iloc[0]["n_matched"], 0)

    def test_missing_star_identifier_produces_only_diagnostics(self):
        self.write_star(self.star_data().drop(columns="rlnImageName"))
        qc, comp = build_summary(self.inputs, self.ref)
        row = qc.iloc[0]
        self.assertIn("rlnImageName", row["primary_error"])
        self.assertEqual(row["match_method"], "defocus_diagnostic")
        self.assertEqual(row["n_defocus_matched"], 1000)
        self.assertEqual(row["n_matched"], 0)
        self.assertTrue(pd.isna(row["n_unique_particles"]))
        self.assertTrue(comp.empty)

    def test_either_missing_cs_identifier_produces_only_diagnostics(self):
        for missing in ["blob/path", "blob/idx"]:
            with self.subTest(missing=missing):
                data = self.cs_data()
                self.write_cs(data[[name for name in data.dtype.names if name != missing]])
                qc, comp = build_summary(self.inputs, self.ref)
                self.assertIn(missing, qc.iloc[0]["primary_error"])
                self.assertEqual(qc.iloc[0]["n_defocus_matched"], 1000)
                self.assertEqual(qc.iloc[0]["n_matched"], 0)
                self.assertTrue(comp.empty)

    def test_invalid_identifiers_produce_only_diagnostics(self):
        for identifier in ["broken", "1@particles/renamed.mrcs", "0@0001.mrcs", np.nan]:
            with self.subTest(identifier=identifier):
                data = self.star_data()
                data.loc[0, "rlnImageName"] = identifier
                self.write_star(data)
                qc, comp = build_summary(self.inputs, self.ref)
                self.assertIn("Invalid particle identifier", qc.iloc[0]["primary_error"])
                self.assertEqual(qc.iloc[0]["n_defocus_matched"], 1000)
                self.assertTrue(comp.empty)

    def test_diagnostic_requires_both_defocus_columns(self):
        for remaining in ["rlnDefocusU", "rlnDefocusV", None]:
            with self.subTest(remaining=remaining):
                data = self.star_data().drop(columns="rlnImageName")
                data["other_column"] = 1
                drop = [name for name in ["rlnDefocusU", "rlnDefocusV"] if name != remaining]
                self.write_star(data.drop(columns=drop))
                qc, comp = build_summary(self.inputs, self.ref)
                self.assertIn("both defocus columns", qc.iloc[0]["defocus_note"])
                self.assertTrue(pd.isna(qc.iloc[0]["n_defocus_matched"]))
                self.assertTrue(comp.empty)

    def test_diagnostic_counts_invalid_and_unmatched_defocus(self):
        data = self.star_data().drop(columns="rlnImageName")
        data["rlnDefocusU"] = data["rlnDefocusU"].astype(object)
        data.loc[0, "rlnDefocusU"] = "invalid"
        data.loc[1, "rlnDefocusV"] = np.inf
        data.loc[2, "rlnDefocusU"] = 99999
        self.write_star(data)
        qc, comp = build_summary(self.inputs, self.ref)
        self.assertEqual(qc.iloc[0]["n_defocus_invalid_rows"], 2)
        self.assertEqual(qc.iloc[0]["n_defocus_unique"], 998)
        self.assertEqual(qc.iloc[0]["n_defocus_matched"], 997)
        self.assertEqual(qc.iloc[0]["n_defocus_unmatched"], 1)
        self.assertTrue(comp.empty)

    def test_unmatched_identifiers_are_not_replaced_by_defocus(self):
        for n_unmatched in [1, 1000]:
            with self.subTest(n_unmatched=n_unmatched):
                data = self.star_data()
                data.loc[:n_unmatched - 1, "rlnImageName"] = [
                    f"{i + 1}@particles/{9000 + i % 4}.mrcs" for i in range(n_unmatched)]
                self.write_star(data)
                qc, comp = build_summary(self.inputs, self.ref)
                self.assertEqual(qc.iloc[0]["match_method"], "path_idx")
                self.assertEqual(qc.iloc[0]["n_unmatched"], n_unmatched)
                self.assertTrue(pd.isna(qc.iloc[0]["n_defocus_matched"]))
                self.assertEqual(comp["count"].sum(), 1000 - n_unmatched)

    def test_evaluation_logs_diagnostics_then_stops_even_with_valid_files(self):
        self.write_star(self.star_data(), name="valid_class1.star")
        self.write_star(self.star_data().drop(columns="rlnImageName"), name="invalid_class2.star")
        stdout = io.StringIO()
        with patch("comphet.evaluate.class_level_table") as score_input:
            with contextlib.redirect_stdout(stdout), self.assertRaisesRegex(SystemExit, "Evaluation stopped"):
                run(self.inputs, self.reference_path, self.output)
            score_input.assert_not_called()
        self.assertIn("1,000 of 1,000 unique keys matched", stdout.getvalue())
        self.assertIn("not used for evaluation", stdout.getvalue())
        self.assertEqual([p.name for p in self.output.iterdir()], ["summary_qc.csv"])
        qc = pd.read_csv(self.output / "summary_qc.csv")
        diagnostic = qc[qc["match_method"] == "defocus_diagnostic"].iloc[0]
        self.assertEqual(diagnostic["n_matched"], 0)
        self.assertEqual(diagnostic["n_defocus_matched"], 1000)

    def test_evaluation_logs_unavailable_diagnostic_before_stopping(self):
        self.write_star(self.star_data()[["rlnDefocusU"]])
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout), self.assertRaisesRegex(SystemExit, "Evaluation stopped"):
            run(self.inputs, self.reference_path, self.output)
        self.assertIn("Diagnostic unavailable", stdout.getvalue())
        self.assertEqual([p.name for p in self.output.iterdir()], ["summary_qc.csv"])

    def test_duplicate_and_unmatched_particles_warn_and_evaluation_continues(self):
        data = self.star_data()
        unmatched = data.iloc[:4].copy()
        unmatched["rlnImageName"] = [f"1@{9000 + i}.mrcs" for i in range(4)]
        self.write_star(pd.concat([data, data.iloc[[0, 0, 1]], unmatched], ignore_index=True))
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            run(self.inputs, self.reference_path, self.output, dpi=40)
        self.assertIn("WARNING: example/class1.star", stdout.getvalue())
        self.assertIn("3 duplicate row(s)", stdout.getvalue())
        self.assertIn("4 of 1,004 unique particle identifiers did not match", stdout.getvalue())
        qc = pd.read_csv(self.output / "summary_qc.csv").iloc[0]
        self.assertEqual(qc["n_rows_raw"], 1007)
        self.assertEqual(qc["n_duplicated_rows"], 3)
        self.assertEqual(qc["n_unmatched"], 4)
        self.assertEqual(qc["n_matched"], 1000)
        self.assertIn("3 duplicate row(s)", qc["warning"])
        self.assertIn("did not match", qc["warning"])
        self.assertEqual(pd.read_csv(self.output / "summary_composition.csv")["count"].sum(), 1000)
        self.assertEqual(pd.read_csv(self.output / "summary_scores.csv").iloc[0]["total_score"], 25)
        self.assertEqual(len(list(self.output.iterdir())), 6)

    def test_shared_particles_warn_only_within_the_same_submission(self):
        data = self.star_data()
        self.write_star(data.iloc[:100], name="class1.star")
        self.write_star(data.iloc[50:150], name="class2.star")
        other = self.inputs / "other_submission"
        other.mkdir()
        starfile.write({"particles": data.iloc[:100]}, other / "class1.star")
        qc, comp = build_summary(self.inputs, self.ref)
        shared = qc.set_index(["team", "class_label"])["n_shared_particles"].to_dict()
        self.assertEqual(shared, {("example", "class1"): 50, ("example", "class2"): 50,
                                  ("other_submission", "class1"): 0})
        for message in qc.loc[qc["team"] == "example", "warning"]:
            self.assertIn("50 particle identifier(s) also occur in other class files", message)
        self.assertEqual(qc.loc[qc["team"] == "other_submission", "warning"].iloc[0], "")
        self.assertEqual(comp.groupby(["team", "class_label"])["count"].sum().to_dict(),
                         {("example", "class1"): 100, ("example", "class2"): 100,
                          ("other_submission", "class1"): 100})

    def test_fully_unmatched_class_warns_without_blocking_valid_class(self):
        data = self.star_data()
        self.write_star(data, name="valid.star")
        data["rlnImageName"] = [f"{i + 1}@{9000 + i % 4}.mrcs" for i in range(len(data))]
        self.write_star(data, name="unmatched.star")
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            run(self.inputs, self.reference_path, self.output, dpi=40)
        self.assertIn("WARNING: example/unmatched.star", stdout.getvalue())
        self.assertIn("No matched particles in this file", stdout.getvalue())
        scores = pd.read_csv(self.output / "summary_scores.csv")
        self.assertEqual(scores.iloc[0]["PKM2_class"], "valid")
        self.assertEqual(scores.iloc[0]["total_score"], 25)
        qc = pd.read_csv(self.output / "summary_qc.csv").set_index("class_label")
        self.assertEqual(qc.loc["unmatched", "n_unmatched"], 1000)
        self.assertEqual(len(list(self.output.iterdir())), 6)

    def test_no_primary_matches_completes_with_explicit_empty_results(self):
        data = self.star_data()
        data["rlnImageName"] = [f"{i + 1}@{9000 + i % 4}.mrcs" for i in range(len(data))]
        self.write_star(data)
        self.output.mkdir()
        # Successful completion must replace stale scores and plots too.
        (self.output / "summary_scores.csv").write_text("old scores\n")
        (self.output / "composition_heatmap.png").write_text("old figure\n")
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            run(self.inputs, self.reference_path, self.output, dpi=40)
        self.assertIn("WARNING: example/class1.star", stdout.getvalue())
        self.assertIn("No scores available", stdout.getvalue())
        self.assertEqual(len(list(self.output.iterdir())), 6)
        self.assertTrue(pd.read_csv(self.output / "summary_scores.csv").empty)
        self.assertTrue(pd.read_csv(self.output / "summary_composition.csv").empty)
        qc = pd.read_csv(self.output / "summary_qc.csv").iloc[0]
        self.assertEqual(qc["n_unmatched"], 1000)
        for filename in ["composition_heatmap.png", "submission_small_multiples.png", "exposure_split.png"]:
            self.assertTrue((self.output / filename).read_bytes().startswith(b"\x89PNG"))
        for number in plt.get_fignums():
            text = [item.get_text() for ax in plt.figure(number).axes for item in ax.texts]
            self.assertTrue(any("No particles matched the reference" in item for item in text))


if __name__ == "__main__":
    unittest.main()
