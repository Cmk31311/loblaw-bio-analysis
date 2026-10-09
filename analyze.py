"""Run Parts 2-4 from SQLite, export reproducible outputs, print key answers."""
from __future__ import annotations

import json

from database import (
    ROOT, baseline_breakdown, baseline_project_counts, baseline_samples,
    frequencies, male_responder_baseline_b_mean,
)
from stats import TIMEPOINTS, cohort_frequencies, compare, written_conclusion

OUTPUTS = ROOT / "outputs"


def main() -> None:
    OUTPUTS.mkdir(exist_ok=True)
    sample_frequencies = frequencies()
    sample_frequencies.to_csv(OUTPUTS / "sample_frequencies.csv", index=False)
    base = baseline_samples()
    base.to_csv(OUTPUTS / "baseline_miraclib_pbmc.csv", index=False)
    projects = baseline_project_counts()
    projects.to_csv(OUTPUTS / "baseline_counts_project.csv", index=False)
    for dimension in ("response", "sex"):
        baseline_breakdown(dimension).to_csv(OUTPUTS / f"baseline_counts_{dimension}.csv", index=False)

    cohort = cohort_frequencies()
    statistical_results = {}
    comparisons = {}
    for name, day in TIMEPOINTS.items():
        tag = "all_visits" if day is None else f"day_{day}"
        result = compare(cohort, day=day)
        result.to_csv(OUTPUTS / f"responder_comparison_{tag}.csv", index=False)
        statistical_results[tag] = result.loc[result["significant"], "population"].tolist()
        comparisons[name] = result

    avg, n = male_responder_baseline_b_mean()
    summary = {
        "sample_count": int(sample_frequencies["sample"].nunique()),
        "frequency_rows": int(len(sample_frequencies)),
        "eligible_baseline_samples": int(len(base)),
        "baseline_subjects": int(base["subject"].nunique()),
        "b_cell_mean_male_melanoma_responder_day_0": round(avg, 2),
        "b_cell_mean_qualifying_samples": n,
        "baseline_samples_by_project": {
            str(row.project): int(row.samples) for row in projects.itertuples(index=False)
        },
        "significant_populations_by_timepoint": statistical_results,
        "conclusions": {
            name: written_conclusion(result, name) for name, result in comparisons.items()
        },
    }
    (OUTPUTS / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (OUTPUTS / "answer.txt").write_text(f"{avg:.2f}\n")
    print(f"Part 2: {len(sample_frequencies):,} frequency rows from {summary['sample_count']:,} samples")
    print(f"Part 3: completed 5 tests at each of 4 timepoint settings")
    print(f"Part 4: {len(base)} qualifying baseline samples")
    print(f"Part 4: average B-cell count = {avg:.2f} (n={n})")
    print("Part 4 projects: " + ", ".join(f"{row.project}={int(row.samples)}" for row in projects.itertuples(index=False)))
    for name, result in comparisons.items():
        print()
        print(written_conclusion(result, name))
    print(f"Exports: {OUTPUTS}")


if __name__ == "__main__":
    main()
