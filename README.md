# Loblaw Bio Analysis

Loads `cell-count.csv` into SQLite, exports the frequency and response tables, and opens a Streamlit app that queries that database.

Python 3.10 or newer. No API keys.

## Run

```bash
make setup
make pipeline
make test
make dashboard
```

Then open http://127.0.0.1:8501.

`make setup` installs `requirements.txt`. `make pipeline` runs `load_data.py` and `analyze.py`. `make dashboard` starts Streamlit on `0.0.0.0:8501`.

If `python` is not on `PATH`, the Makefile uses `python3.11`, then `python3`. `python load_data.py` takes no arguments and finds the CSV next to the script, so the working directory does not matter.

## GitHub Codespaces

1. Push the repo, including `cell-count.csv`. Leave `clinical_trial.db` and `outputs/` untracked.
2. Code → Codespaces → Create codespace.
3. In the codespace terminal, from the repo root, run the four `make` commands above.
4. Open port 8501 from the Ports tab. If the page is blank, set the port to public.
5. Stop the app with Ctrl+C.

`.streamlit/config.toml` sets headless mode and turns off CORS and XSRF checks so the forwarded port can connect.

## Layout

```text
cell-count.csv → load_data.py + schema.sql → clinical_trial.db
                                              ├─ analyze.py → outputs/
                                              └─ app.py
```

| File | What it does |
|---|---|
| `load_data.py` | Checks the CSV and replaces `clinical_trial.db`. Running it again does not duplicate rows. |
| `schema.sql` | Tables, foreign keys, and the frequency views. |
| `database.py` | Queries used by the analysis and the app. |
| `stats.py` | Patient-level Welch tests, intervals, effect sizes, and FDR. |
| `analyze.py` | Writes the CSV outputs. |
| `app.py` | Dashboard. Numbers come from the database at view time. |
| `tests/test_data.py` | Counts, filters, and the statistical checks. |

`projects` → `subjects` → `samples` → `cell_counts`, plus `cell_populations`. The view `sample_frequencies` has `sample`, `total_count`, `population`, `count`, `percentage`.

The CSV columns are `project`, `subject`, `condition`, `age`, `sex`, `treatment`, `response`, `sample`, `sample_type`, `time_from_treatment_start`, `b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, `monocyte`.

## Methods

Frequency for a population is `100 * count / total_count`, where `total_count` is the sum of the five counts in that sample.

The response comparison keeps melanoma, miraclib, PBMC, and response `yes` or `no`. Blank responses stay NULL and are left out of that comparison. They are not coded as non-responders.

Each patient is one observation. The combined view averages that patient's visit percentages. Day 0, 7, and 14 each use the sample from that day. The test is a two-sided Welch t-test. The difference is responder mean minus non-responder mean, in percentage points, with a 95% Welch interval. Cohen's d uses the pooled SD. Hedges' g applies `1 - 3 / (4*(n1+n2) - 9)`. Benjamini-Hochberg is applied to the five populations in the selected comparison, not across timepoints. Adjusted p < 0.05 is the cutoff used in the app.

The day-0 subset is melanoma, miraclib, PBMC, day 0. Project counts are samples, including a project with zero rows. Response and sex counts are distinct patients.

The B-cell mean is a different question: melanoma, male, response yes, day 0, every treatment and sample type. It averages the absolute B-cell count and is printed to two decimals.

## Assumptions

Subject fields are constant across visits. The loader rejects conflicts. Condition, treatment, and response are stored lowercase. Sex is `M` or `F`. Sample type is uppercase.

The five measured populations are the denominator. A higher share for one population is a lower share for the others, so the five tests are not independent. FDR does not fix that. Age, sex, and project are not covariates.

A difference is an association with the recorded response. Day 0 is before treatment. Day 7 and day 14 are on treatment, so they are not baseline predictors. There is no fitted prediction model.

Treatments in the file are `miraclib`, `phauximab`, and `none`.

## Outputs

`make pipeline` writes `clinical_trial.db` and:

- `outputs/sample_frequencies.csv` (52,500 rows)
- `outputs/responder_comparison_all_visits.csv`
- `outputs/responder_comparison_day_0.csv`
- `outputs/responder_comparison_day_7.csv`
- `outputs/responder_comparison_day_14.csv`
- `outputs/baseline_miraclib_pbmc.csv`
- `outputs/baseline_counts_project.csv`
- `outputs/baseline_counts_response.csv`
- `outputs/baseline_counts_sex.csv`
- `outputs/answer.txt`
- `outputs/analysis_summary.json`

The database and `outputs/` are gitignored. `make pipeline` recreates them.

For this CSV the checks are 10,500 samples, 3,500 subjects, 52,500 count rows, 656 melanoma/miraclib/PBMC patients (1,968 measurements), day-0 samples 656 (prj1 384, prj2 0, prj3 272), response 331/325, sex 344 male and 312 female, and a male melanoma responder day-0 B-cell mean of 10206.15 on 485 rows.

In the combined visit comparison, CD4 T cells are higher in responders by 0.64 percentage points (95% CI 0.20 to 1.07, Hedges' g 0.22, adjusted p 0.0226). The other four populations are not significant after correction. Day 0 and day 14 have no population under the cutoff. Day 7 CD4 T cells have adjusted p 0.0498, which is just under 0.05, with Hedges' g 0.20.

## Deploy

The app is not deployed. To host it, push the repo and create a Streamlit Community Cloud app with main file `app.py`. If `clinical_trial.db` is absent, the app builds it from the CSV. Put the real URL here after that works.

## Limits

The tests do not adjust for age, sex, or project, and they do not use a mixed model for the three visits. Multiplicity control is only within the five populations of one comparison.
