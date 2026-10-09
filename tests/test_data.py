"""Automated analytical acceptance tests using the included source CSV."""
from __future__ import annotations

import csv
from contextlib import closing
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from database import (
    baseline_breakdown, baseline_project_counts, baseline_samples, connect,
    filtered_frequencies, frequencies, male_responder_baseline_b_mean, read_frame,
)
from load_data import CSV_PATH, HEADERS, load
from stats import (
    cohort_frequencies, compare, independent_values, standardized_mean_difference,
    written_conclusion,
)


def _benjamini_hochberg(p_values):
    """BH adjusted p-values, implemented here so the test does not call statsmodels."""
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adjusted = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n)
    out[order] = adjusted
    return out


class TrialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workdir = tempfile.TemporaryDirectory()
        cls.db = Path(cls.workdir.name) / "test_clinical_trial.db"
        cls.stats = load(CSV_PATH, cls.db)

    @classmethod
    def tearDownClass(cls):
        cls.workdir.cleanup()

    def test_loading_counts(self):
        self.assertEqual(self.stats, {"projects": 3, "subjects": 3500,
                                      "samples": 10500, "cell_counts": 52500})
        with closing(connect(self.db)) as con:
            self.assertFalse(con.execute("PRAGMA foreign_key_check").fetchall())

    def test_frequency_shape_and_totals(self):
        data = frequencies(self.db)
        self.assertEqual(list(data.columns), ["sample", "total_count", "population", "count", "percentage"])
        self.assertEqual(data.shape, (52500, 5))
        grouped = data.groupby("sample").agg(rows=("population", "size"),
                                               percentage_total=("percentage", "sum"))
        self.assertTrue((grouped["rows"] == 5).all())
        self.assertTrue(np.allclose(grouped["percentage_total"], 100.0))
        self.assertTrue((data["total_count"] > 0).all())
        one = data[data["sample"] == data.iloc[0]["sample"]]
        self.assertTrue(np.allclose(one["percentage"], 100.0 * one["count"] / one["total_count"]))

    def test_baseline_filters(self):
        base = baseline_samples(self.db)
        self.assertEqual(len(base), 656)
        self.assertEqual(base["subject"].nunique(), 656)
        self.assertTrue((base["condition"] == "melanoma").all())
        self.assertTrue((base["treatment"] == "miraclib").all())
        self.assertTrue((base["sample_type"] == "PBMC").all())
        self.assertTrue((base["time_from_treatment_start"] == 0).all())

    def test_baseline_breakdowns(self):
        projects = dict(baseline_breakdown("project", self.db).values)
        responses = dict(baseline_breakdown("response", self.db).values)
        sexes = dict(baseline_breakdown("sex", self.db).values)
        self.assertEqual(projects, {"prj1": 384, "prj3": 272})
        self.assertEqual(responses, {"no": 325, "yes": 331})
        self.assertEqual(sexes, {"F": 312, "M": 344})

    def test_b_cell_mean_all_treatments_and_types(self):
        mean, count = male_responder_baseline_b_mean(self.db)
        self.assertEqual(count, 485)
        self.assertEqual(f"{mean:.2f}", "10206.15")

    def test_cohort_statistical_unit_is_subject(self):
        raw = cohort_frequencies(self.db)
        self.assertEqual(raw["sample"].nunique(), 1968)
        averages = independent_values(raw)
        self.assertEqual(averages["subject_id"].nunique(), 656)
        self.assertEqual(len(averages), 656 * 5)
        results = compare(raw)
        self.assertEqual(len(results), 5)
        self.assertTrue((results["responder_n"] == 331).all())
        self.assertTrue((results["non_responder_n"] == 325).all())
        self.assertTrue(results["adjusted_p_value"].between(0, 1).all())
        self.assertTrue((results["adjusted_p_value"] + 1e-8 >= results["p_value"]).all())
        self.assertTrue(np.isfinite(results["cohens_d"]).all())
        self.assertTrue(np.isfinite(results["hedges_g"]).all())
        self.assertTrue((np.sign(results["hedges_g"]) == np.sign(results["mean_difference"])).all())
        self.assertTrue((results["ci_lower"] <= results["mean_difference"]).all())
        self.assertTrue((results["mean_difference"] <= results["ci_upper"]).all())
        self.assertEqual(set(results.loc[results["significant"], "population"]), {"cd4_t_cell"})
        conclusion = written_conclusion(results, "All visits (patient average)")
        self.assertIn("331 responders", conclusion)
        self.assertIn("325 non-responders", conclusion)
        self.assertIn("CD4 T cells", conclusion)
        self.assertIn("does not show that miraclib caused", conclusion)
        cd4 = results.set_index("population").loc["cd4_t_cell"]
        self.assertAlmostEqual(float(cd4["mean_difference"]), 0.6355, places=3)
        self.assertAlmostEqual(float(cd4["ci_lower"]), 0.1976, places=3)
        self.assertAlmostEqual(float(cd4["ci_upper"]), 1.0735, places=3)
        self.assertAlmostEqual(float(cd4["hedges_g"]), 0.2222, places=3)
        self.assertAlmostEqual(float(cd4["adjusted_p_value"]), 0.02257, places=4)
        self.assertTrue(np.allclose(_benjamini_hochberg(results["p_value"].to_numpy()), results["adjusted_p_value"]))
        day7 = compare(raw, day=7)
        self.assertEqual(set(day7.loc[day7["significant"], "population"]), {"cd4_t_cell"})
        day7_p = float(day7.set_index("population").loc["cd4_t_cell", "adjusted_p_value"])
        self.assertGreater(day7_p, 0.04)
        self.assertLess(day7_p, 0.05)
        self.assertFalse(compare(raw, day=14)["significant"].any())
        day0 = compare(raw, day=0)
        self.assertFalse(day0["significant"].any())
        self.assertIn("no population differs", written_conclusion(day0, "Day 0"))
        for day in [0, 7, 14]:
            day_result = compare(raw, day=day)
            self.assertEqual(len(day_result), 5)
            self.assertTrue(np.isfinite(day_result["hedges_g"]).all())

    def test_reloading_does_not_duplicate(self):
        load(CSV_PATH, self.db)
        with closing(connect(self.db)) as con:
            self.assertEqual(con.execute("SELECT COUNT(*) FROM samples").fetchone()[0], 10500)
            self.assertEqual(con.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0], 52500)

    def test_project_counts_include_explicit_zero(self):
        projects = dict(baseline_project_counts(self.db).values)
        self.assertEqual(projects, {"prj1": 384, "prj2": 0, "prj3": 272})
        self.assertEqual(int(baseline_project_counts(self.db)["samples"].sum()), 656)

    def test_b_cell_query_spans_treatments_and_sample_types(self):
        frame = read_frame(
            """
            SELECT COUNT(DISTINCT su.treatment) AS treatments,
                   COUNT(DISTINCT s.sample_type) AS sample_types
            FROM cell_counts cc
            JOIN samples s ON s.sample_id = cc.sample_id
            JOIN subjects su ON su.subject_id = s.subject_id
            WHERE cc.population = 'b_cell' AND su.condition = 'melanoma' AND su.sex = 'M'
              AND su.response = 'yes' AND s.time_from_treatment_start = 0
            """,
            db_path=self.db,
        )
        self.assertGreater(int(frame.loc[0, "treatments"]), 1)
        self.assertGreater(int(frame.loc[0, "sample_types"]), 1)

    def test_filtered_frequencies_and_search(self):
        everything = filtered_frequencies({}, db_path=self.db)
        self.assertEqual(list(everything.columns), ["sample", "total_count", "population", "count", "percentage"])
        self.assertEqual(len(everything), 52500)
        b_cells = filtered_frequencies({"population": "b_cell"}, db_path=self.db)
        self.assertEqual(len(b_cells), 10500)
        self.assertTrue((b_cells["population"] == "b_cell").all())
        sample = everything.iloc[0]["sample"]
        found = filtered_frequencies({}, search=sample, db_path=self.db)
        self.assertEqual(found["sample"].nunique(), 1)
        self.assertEqual(len(found), 5)
        injected = filtered_frequencies({}, search="' OR 1=1 --", db_path=self.db)
        self.assertTrue(injected.empty)
        unrecorded = filtered_frequencies({"response": "unrecorded"}, db_path=self.db)
        self.assertEqual(len(unrecorded), 7110)

    def test_source_metadata_domains(self):
        with closing(connect(self.db)) as con:
            sexes = {row[0] for row in con.execute("SELECT DISTINCT sex FROM subjects")}
            conditions = {row[0] for row in con.execute("SELECT DISTINCT condition FROM subjects")}
            treatments = {row[0] for row in con.execute("SELECT DISTINCT treatment FROM subjects")}
            sample_types = {row[0] for row in con.execute("SELECT DISTINCT sample_type FROM samples")}
            days = {row[0] for row in con.execute("SELECT DISTINCT time_from_treatment_start FROM samples")}
            populations = {row[0] for row in con.execute("SELECT population FROM cell_populations")}
            incomplete_samples = con.execute(
                "SELECT COUNT(*) FROM (SELECT sample_id FROM cell_counts GROUP BY sample_id HAVING COUNT(*) != 5)"
            ).fetchone()[0]
            incomplete_subjects = con.execute(
                "SELECT COUNT(*) FROM (SELECT subject_id FROM samples GROUP BY subject_id HAVING COUNT(*) != 3)"
            ).fetchone()[0]
            null_response = con.execute("SELECT COUNT(*) FROM subjects WHERE response IS NULL").fetchone()[0]
            response_no = con.execute("SELECT COUNT(*) FROM subjects WHERE response = 'no'").fetchone()[0]
        self.assertEqual(sexes, {"M", "F"})
        self.assertEqual(conditions, {"melanoma", "carcinoma", "healthy"})
        self.assertEqual(treatments, {"miraclib", "phauximab", "none"})
        self.assertEqual(sample_types, {"PBMC", "WB"})
        self.assertEqual(days, {0, 7, 14})
        self.assertEqual(populations, {"b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"})
        self.assertEqual(incomplete_samples, 0)
        self.assertEqual(incomplete_subjects, 0)
        self.assertEqual(null_response, 474)
        self.assertEqual(response_no, 1489)

    def test_breakdown_rejects_unknown_dimension(self):
        with self.assertRaises(ValueError):
            baseline_breakdown("gender", self.db)

    def test_effect_size_known_values(self):
        cohens_d, hedges_g = standardized_mean_difference(
            np.array([1.0, 2.0, 3.0]), np.array([1.0, 1.0, 1.0])
        )
        self.assertAlmostEqual(cohens_d, 1.41421356237, places=6)
        self.assertAlmostEqual(hedges_g, 1.13137084990, places=6)
        zero_d, zero_g = standardized_mean_difference(np.array([1.0, 1.0]), np.array([1.0, 1.0]))
        self.assertEqual((zero_d, zero_g), (0.0, 0.0))

    def test_conclusion_does_not_invent_significance(self):
        quiet = pd.DataFrame([{
            "population": "b_cell", "responder_n": 10, "non_responder_n": 12,
            "responder_mean": 1.0, "non_responder_mean": 1.2, "mean_difference": -0.2,
            "ci_lower": -0.5, "ci_upper": 0.1, "cohens_d": -0.1, "hedges_g": -0.09,
            "p_value": 0.2, "adjusted_p_value": 0.4, "significant": False,
        }])
        text = written_conclusion(quiet, "Day 0")
        self.assertIn("no population differs", text)
        self.assertIn("10 responders", text)
        self.assertNotIn("B cells are lower", text)

    def test_loader_rejects_invalid_sources(self):
        folder = Path(self.workdir.name)

        def write(name, rows):
            path = folder / name
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(HEADERS))
                writer.writeheader()
                writer.writerows(rows)
            return path

        def row(**overrides):
            base = {key: "0" for key in HEADERS}
            base.update({
                "project": "prj1", "subject": "s1", "condition": "melanoma", "age": "40",
                "sex": "F", "treatment": "miraclib", "response": "yes", "sample": "a",
                "sample_type": "PBMC", "time_from_treatment_start": "0",
                "b_cell": "10", "cd8_t_cell": "10", "cd4_t_cell": "10",
                "nk_cell": "10", "monocyte": "10",
            })
            base.update(overrides)
            return base

        with self.assertRaises(ValueError):
            load(write("dup.csv", [row(), row(subject="s2")]), folder / "dup.db")
        with self.assertRaises(ValueError):
            load(write("zero.csv", [row(
                b_cell="0", cd8_t_cell="0", cd4_t_cell="0", nk_cell="0", monocyte="0"
            )]), folder / "zero.db")
        with self.assertRaises(ValueError):
            load(write("neg.csv", [row(b_cell="-1")]), folder / "neg.db")
        with self.assertRaises(ValueError):
            load(write("meta.csv", [
                row(sample="a1"),
                row(sample="a2", age="41", time_from_treatment_start="7"),
            ]), folder / "meta.db")
        bad = folder / "bad.csv"
        bad.write_text("sample_id,indication\n1,melanoma\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            load(bad, folder / "bad.db")


if __name__ == "__main__":
    unittest.main()
