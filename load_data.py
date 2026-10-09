"""Rebuild a validated SQLite database from cell-count.csv.

Run directly from any working directory: python load_data.py
"""
from __future__ import annotations

import csv
from contextlib import closing
import os
import sqlite3
import tempfile
from pathlib import Path

from database import DATABASE, ROOT

CSV_PATH = ROOT / "cell-count.csv"
SCHEMA_PATH = ROOT / "schema.sql"
POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")
HEADERS = (
    "project", "subject", "condition", "age", "sex", "treatment", "response",
    "sample", "sample_type", "time_from_treatment_start", *POPULATIONS
)


def load(csv_path: Path = CSV_PATH, database_path: Path = DATABASE) -> dict:
    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    projects: set[str] = set()
    subjects: dict[str, tuple] = {}
    samples: list[tuple] = []
    counts: list[tuple] = []
    used_samples: set[str] = set()

    with csv_path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or set(reader.fieldnames) != set(HEADERS):
            raise ValueError(f"CSV columns differ from expected fields: {reader.fieldnames}")
        for index, row in enumerate(reader, start=2):
            project = row["project"].strip()
            subject = row["subject"].strip()
            sample = row["sample"].strip()
            if not all((project, subject, sample)):
                raise ValueError(f"Missing project, subject or sample on CSV line {index}")
            if sample in used_samples:
                raise ValueError(f"Duplicate sample id {sample!r} on CSV line {index}")
            used_samples.add(sample)
            age = int(row["age"])
            day = int(row["time_from_treatment_start"])
            response = row["response"].strip().lower() or None
            subject_data = (
                subject, project, row["condition"].strip().lower(), age,
                row["sex"].strip().upper(), row["treatment"].strip().lower(), response
            )
            if subject in subjects and subjects[subject] != subject_data:
                raise ValueError(f"Inconsistent metadata for subject {subject!r}")
            subjects[subject] = subject_data
            projects.add(project)
            samples.append((sample, subject, row["sample_type"].strip().upper(), day))
            cell_values = [int(row[population]) for population in POPULATIONS]
            if any(n < 0 for n in cell_values) or sum(cell_values) == 0:
                raise ValueError(f"Invalid cell counts for sample {sample!r}")
            counts.extend((sample, population, n) for population, n in zip(POPULATIONS, cell_values))

    if not samples:
        raise ValueError("CSV is empty")

    # Build a complete database in a temporary file and only replace the
    # previous version after successful validation and commit.
    database_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=".clinical_trial_", suffix=".db", dir=database_path.parent)
    os.close(fd)
    temp_path = Path(temp_name)
    try:
        with closing(sqlite3.connect(temp_path)) as con:
            con.execute("PRAGMA foreign_keys = ON")
            con.executescript(SCHEMA_PATH.read_text(encoding="utf-8"))
            with con:
                con.executemany("INSERT INTO projects(project_id) VALUES (?)", [(p,) for p in sorted(projects)])
                con.executemany(
                    "INSERT INTO cell_populations(population) VALUES (?)", [(p,) for p in POPULATIONS]
                )
                con.executemany("INSERT INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)", subjects.values())
                con.executemany("INSERT INTO samples VALUES (?, ?, ?, ?)", samples)
                con.executemany("INSERT INTO cell_counts VALUES (?, ?, ?)", counts)
            errors = con.execute("PRAGMA foreign_key_check").fetchall()
            if errors:
                raise ValueError(f"Foreign key errors: {errors[:5]}")
            count_rows = con.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0]
            if count_rows != len(samples) * len(POPULATIONS):
                raise ValueError("Loaded cell count record total does not match source rows")
        os.replace(temp_path, database_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()

    return {"projects": len(projects), "subjects": len(subjects),
            "samples": len(samples), "cell_counts": len(counts)}


if __name__ == "__main__":
    totals = load()
    print(f"Database: {DATABASE}")
    print("Loaded: " + ", ".join(f"{key}={val:,}" for key, val in totals.items()))
