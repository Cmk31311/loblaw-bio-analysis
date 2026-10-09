"""Repository-relative SQLite connections and reusable parameterized queries."""
from __future__ import annotations

import sqlite3
from contextlib import closing
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DATABASE = ROOT / "clinical_trial.db"


def connect(db_path: Path = DATABASE) -> sqlite3.Connection:
    if not db_path.is_file():
        raise FileNotFoundError(
            f"Database missing: {db_path}. Run `make pipeline` (or `python load_data.py`) first."
        )
    con = sqlite3.connect(db_path)
    # Not stored in the file; each connection has to turn this on.
    con.execute("PRAGMA foreign_keys = ON")
    return con


def read_frame(query: str, params: tuple = (), db_path: Path = DATABASE) -> pd.DataFrame:
    with closing(connect(db_path)) as con:
        return pd.read_sql_query(query, con, params=params)


def frequencies(db_path: Path = DATABASE) -> pd.DataFrame:
    return read_frame(
        "SELECT sample, total_count, population, count, percentage "
        "FROM sample_frequencies ORDER BY sample, population", db_path=db_path
    )


_FREQUENCY_COLUMNS = {
    "project": "project",
    "condition": "condition",
    "treatment": "treatment",
    "response": "response",
    "sex": "sex",
    "sample_type": "sample_type",
    "time_from_treatment_start": "time_from_treatment_start",
    "population": "population",
}


def filtered_frequencies(
    filters: dict | None = None, search: str = "", db_path: Path = DATABASE
) -> pd.DataFrame:
    """Five required frequency columns, with optional exact filters and an id search."""
    clauses: list[str] = []
    params: list = []
    for key, value in (filters or {}).items():
        column = _FREQUENCY_COLUMNS.get(key)
        if column is None or value is None or value == "All":
            continue
        if key == "response" and value == "unrecorded":
            clauses.append("response IS NULL")
            continue
        clauses.append(f"{column} = ?")
        params.append(value)
    if search.strip():
        clauses.append("(sample LIKE ? OR subject_id LIKE ?)")
        token = f"%{search.strip()}%"
        params.extend((token, token))
    where = (" WHERE " + " AND ".join(clauses)) if clauses else ""
    return read_frame(
        "SELECT sample, total_count, population, count, percentage "
        "FROM frequency_details" + where + " ORDER BY sample, population",
        tuple(params), db_path=db_path,
    )


def baseline_samples(db_path: Path = DATABASE) -> pd.DataFrame:
    """One row for each baseline, melanoma, miraclib, PBMC sample."""
    return read_frame(
        """
        SELECT su.project_id AS project, su.subject_id AS subject,
               s.sample_id AS sample, s.sample_type, s.time_from_treatment_start,
               su.condition, su.treatment, su.response, su.sex
        FROM samples s
        JOIN subjects su ON su.subject_id = s.subject_id
        WHERE su.condition = ? AND su.treatment = ?
          AND s.sample_type = ? AND s.time_from_treatment_start = ?
        ORDER BY su.project_id, su.subject_id, s.sample_id
        """,
        ("melanoma", "miraclib", "PBMC", 0), db_path=db_path
    )


def baseline_breakdown(dimension: str, db_path: Path = DATABASE) -> pd.DataFrame:
    """SQL-backed counts with explicit distinction between samples and subjects."""
    dimensions = {
        "project": ("su.project_id", "COUNT(*)", "samples"),
        "response": ("su.response", "COUNT(DISTINCT su.subject_id)", "subjects"),
        "sex": ("su.sex", "COUNT(DISTINCT su.subject_id)", "subjects"),
    }
    if dimension not in dimensions:
        raise ValueError("dimension must be project, response or sex")
    column, aggregate, label = dimensions[dimension]
    return read_frame(
        f"""
        SELECT {column} AS {dimension}, {aggregate} AS {label}
        FROM samples s JOIN subjects su ON su.subject_id = s.subject_id
        WHERE su.condition = ? AND su.treatment = ?
          AND s.sample_type = ? AND s.time_from_treatment_start = ?
        GROUP BY {column} ORDER BY {column}
        """,
        ("melanoma", "miraclib", "PBMC", 0), db_path=db_path
    )


def baseline_project_counts(db_path: Path = DATABASE) -> pd.DataFrame:
    """Qualifying baseline samples in every project, including explicit zeros."""
    return read_frame(
        """
        SELECT p.project_id AS project, COUNT(s.sample_id) AS samples
        FROM projects p
        LEFT JOIN subjects su ON su.project_id = p.project_id
          AND su.condition = ? AND su.treatment = ?
        LEFT JOIN samples s ON s.subject_id = su.subject_id
          AND s.sample_type = ? AND s.time_from_treatment_start = ?
        GROUP BY p.project_id
        ORDER BY p.project_id
        """,
        ("melanoma", "miraclib", "PBMC", 0), db_path=db_path
    )


def male_responder_baseline_b_mean(db_path: Path = DATABASE) -> tuple[float, int]:
    """DO NOT restrict sample_type or treatment: the assignment asks for all."""
    with closing(connect(db_path)) as con:
        result = con.execute(
            """
            SELECT AVG(cc.count), COUNT(*)
            FROM cell_counts cc
            JOIN samples s ON s.sample_id = cc.sample_id
            JOIN subjects su ON su.subject_id = s.subject_id
            WHERE cc.population = ? AND su.condition = ? AND su.sex = ?
              AND su.response = ? AND s.time_from_treatment_start = ?
            """,
            ("b_cell", "melanoma", "M", "yes", 0),
        ).fetchone()
    if result[0] is None:
        raise ValueError("No samples match the B-cell mean query")
    return float(result[0]), int(result[1])
