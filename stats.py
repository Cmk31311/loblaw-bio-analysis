"""Patient-independent responder comparisons for melanoma/miraclib/PBMC."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

from database import DATABASE, read_frame

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
POPULATION_LABELS = {
    "b_cell": "B cells",
    "cd8_t_cell": "CD8 T cells",
    "cd4_t_cell": "CD4 T cells",
    "nk_cell": "NK cells",
    "monocyte": "monocytes",
}
TIMEPOINTS = {"All visits (patient average)": None, "Day 0": 0, "Day 7": 7, "Day 14": 14}


def cohort_frequencies(db_path: Path = DATABASE) -> pd.DataFrame:
    return read_frame(
        """
        SELECT subject_id, sample, population, percentage, response,
               time_from_treatment_start AS day
        FROM frequency_details
        WHERE condition = ? AND treatment = ? AND sample_type = ?
          AND response IN ('yes', 'no')
        ORDER BY subject_id, day, population
        """,
        ("melanoma", "miraclib", "PBMC"), db_path=db_path
    )


def independent_values(data: pd.DataFrame, day: int | None = None) -> pd.DataFrame:
    """One row per patient/population for each comparison, avoiding pseudoreplication."""
    if day is not None:
        data = data.loc[data["day"] == day]
    if data.empty:
        return pd.DataFrame(columns=["subject_id", "response", "population", "percentage"])
    return data.groupby(["subject_id", "response", "population"], as_index=False)["percentage"].mean()


def standardized_mean_difference(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Cohen's d and Hedges' g for responder mean minus non-responder mean.

    The standardizer is the pooled standard deviation. Hedges' g multiplies d by
    the small-sample correction 1 - 3 / (4*(n1+n2) - 9).
    """
    responder = np.asarray(a, dtype=float)
    non_responder = np.asarray(b, dtype=float)
    n1, n2 = len(responder), len(non_responder)
    difference = float(np.mean(responder) - np.mean(non_responder))
    variance_1 = float(np.var(responder, ddof=1))
    variance_2 = float(np.var(non_responder, ddof=1))
    pooled = ((n1 - 1) * variance_1 + (n2 - 1) * variance_2) / (n1 + n2 - 2)
    if pooled == 0:
        if difference == 0:
            return 0.0, 0.0
        signed = math.copysign(math.inf, difference)
        return signed, signed
    cohens_d = difference / math.sqrt(pooled)
    hedges_g = cohens_d * (1 - (3 / (4 * (n1 + n2) - 9)))
    return float(cohens_d), float(hedges_g)


def effect_size_band(g: float) -> str:
    """Conventional absolute-value bands. Descriptive only, not clinical cutoffs."""
    magnitude = abs(g)
    if magnitude < 0.2:
        return "negligible"
    if magnitude < 0.5:
        return "small"
    if magnitude < 0.8:
        return "medium"
    return "large"


def _welch_ci(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """95% confidence interval for (responder mean - non-responder mean)."""
    na, nb = len(a), len(b)
    va = float(np.var(a, ddof=1))
    vb = float(np.var(b, ddof=1))
    diff = float(np.mean(a) - np.mean(b))
    term_a, term_b = va / na, vb / nb
    se2 = term_a + term_b
    if se2 == 0:
        return diff, diff
    df = se2**2 / (term_a**2 / (na - 1) + term_b**2 / (nb - 1))
    margin = float(scipy_stats.t.ppf(0.975, df) * math.sqrt(se2))
    return diff - margin, diff + margin


def compare(data: pd.DataFrame, day: int | None = None) -> pd.DataFrame:
    values = independent_values(data, day)
    results = []
    for population in POPULATIONS:
        current = values[values["population"] == population]
        a = current[current["response"] == "yes"]["percentage"].to_numpy(dtype=float)
        b = current[current["response"] == "no"]["percentage"].to_numpy(dtype=float)
        if len(a) < 2 or len(b) < 2:
            raise ValueError(f"Too few patient observations in {population} at day={day}")
        result = scipy_stats.ttest_ind(a, b, equal_var=False, alternative="two-sided")
        lower, upper = _welch_ci(a, b)
        cohens_d, hedges_g = standardized_mean_difference(a, b)
        results.append({
            "population": population,
            "responder_n": len(a), "non_responder_n": len(b),
            "responder_mean": float(np.mean(a)), "non_responder_mean": float(np.mean(b)),
            "mean_difference": float(np.mean(a) - np.mean(b)),
            "ci_lower": lower, "ci_upper": upper,
            "cohens_d": cohens_d, "hedges_g": hedges_g,
            "p_value": float(result.pvalue),
        })
    frame = pd.DataFrame(results)
    if frame["p_value"].isna().any():
        raise ValueError("Undefined p-values: check within-group variance")
    frame["adjusted_p_value"] = multipletests(frame["p_value"], method="fdr_bh")[1]
    frame["significant"] = frame["adjusted_p_value"] < 0.05
    return frame


def written_conclusion(results: pd.DataFrame, timepoint_label: str) -> str:
    """Describe the comparison that was actually run. Significance comes only from results."""
    if results.empty:
        return "No patients met the cohort filters, so no comparison was run."
    responders = int(results["responder_n"].iloc[0])
    non_responders = int(results["non_responder_n"].iloc[0])
    if timepoint_label.startswith("All visits"):
        unit = (
            "Each patient contributes one value per population: the mean of that "
            "patient's visit-level relative frequencies. Repeated visits are not "
            "treated as independent patients."
        )
    else:
        unit = (
            f"Each patient contributes the single {timepoint_label} PBMC measurement."
        )
    opening = (
        f"{timepoint_label}: {responders} responders and {non_responders} non-responders "
        f"(melanoma, miraclib, PBMC). {unit} "
        f"Two-sided Welch t-tests, Benjamini-Hochberg adjusted across these five populations "
        f"(adjusted p < 0.05)."
    )
    order = {name: index for index, name in enumerate(POPULATIONS)}
    ranked = results.assign(_order=results["population"].map(order)).sort_values("_order")
    significant = ranked.loc[ranked["significant"]].drop(columns="_order")
    if significant.empty:
        finding = " After correction, no population differs between responders and non-responders."
    else:
        details = []
        for row in significant.itertuples(index=False):
            direction = "higher" if row.mean_difference > 0 else "lower"
            label = POPULATION_LABELS.get(row.population, row.population)
            details.append(
                f"{label} are {direction} in responders "
                f"({row.mean_difference:+.2f} percentage points, "
                f"95% CI {row.ci_lower:.2f} to {row.ci_upper:.2f}, "
                f"Hedges' g {row.hedges_g:+.2f} ({effect_size_band(row.hedges_g)}), "
                f"adjusted p {row.adjusted_p_value:.4f})"
            )
        verb = "differs" if len(significant) == 1 else "differ"
        names = ", ".join(POPULATION_LABELS.get(name, name) for name in significant["population"])
        finding = f" After correction, {len(significant)} of 5 populations {verb}: {names}. " + " ".join(
            item if item.endswith(".") else item + "." for item in details
        )
        if len(significant) < len(ranked):
            others = ", ".join(
                POPULATION_LABELS.get(name, name)
                for name in ranked.loc[~ranked["significant"], "population"]
            )
            finding += f" The other populations ({others}) do not differ after correction."
    caveat = (
        " This does not show that miraclib caused the difference, and it is not a prediction model. "
        "Day 0 is pretreatment. Day 7 and day 14 are on treatment."
    )
    return opening + finding + caveat
