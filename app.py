"""SQLite-backed clinical analytics dashboard. Presentation only; queries live in database.py and stats.py."""
from __future__ import annotations

import html
from datetime import datetime

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from database import (
    DATABASE, baseline_breakdown, baseline_project_counts, baseline_samples,
    filtered_frequencies, male_responder_baseline_b_mean, read_frame,
)
from load_data import load
from stats import (
    POPULATION_LABELS, POPULATIONS, TIMEPOINTS, cohort_frequencies, compare,
    independent_values, written_conclusion,
)

INK = "#101010"
MUTED = "#6F6F6A"
LINE = "#E6E6E3"
GRID = "#EFEFEC"
PAPER = "#F8F8F6"
RED = "#FF3638"
RESPONDER = "#101010"
NONRESPONDER = "#5E738C"
FONT = "Inter, ui-sans-serif, system-ui, sans-serif"
PAGE_SIZE = 100
PLOTLY = {"displaylogo": False, "modeBarButtonsToRemove": ["lasso2d", "select2d"]}
SEX_LABELS = {"M": "Male", "F": "Female"}
TIMEPOINT_SHORT = {"All visits (patient average)": "All visits"}
POP_COLORS = {
    "b_cell": "#161616",
    "cd8_t_cell": "#4E6E8E",
    "cd4_t_cell": "#8C5E4A",
    "nk_cell": "#8A8A84",
    "monocyte": "#C8C1B4",
}
NAV = [
    ("", "Overview"),
    ("frequencies", "Cell Frequencies"),
    ("response", "Response Analysis"),
    ("subset", "Subset Analysis"),
]

st.set_page_config(page_title="Loblaw Bio Analysis", layout="wide", initial_sidebar_state="collapsed")

st.markdown(
    """
    <style>
    [data-testid="stSidebar"], [data-testid="stSidebarCollapsedControl"] {display: none;}
    .block-container {padding-top: 4.75rem; padding-bottom: 3rem; max-width: 1140px;}
    [data-testid="stDataFrame"] {font-variant-numeric: tabular-nums lining-nums;}
    a:focus-visible, button:focus-visible {outline: 2px solid #101010; outline-offset: 2px;}
    .navpill {
        display: flex; align-items: center; justify-content: space-between; gap: 12px;
        background: #FFFFFF; border: 1px solid #E6E6E3; border-radius: 999px;
        padding: 6px 14px 6px 16px; margin: 0 0 1.35rem;
    }
    .brandlink {
        display: inline-flex; align-items: center; gap: 8px; color: #101010;
        text-decoration: none; font-weight: 600; letter-spacing: -0.03em; font-size: 0.92rem;
        white-space: nowrap;
    }
    .dot {width: 7px; height: 7px; border-radius: 50%; background: #FF3638; display: inline-block; flex: 0 0 7px;}
    .navlinks {display: flex; align-items: center; gap: 2px; flex-wrap: wrap; justify-content: center;}
    .navlinks a {
        color: #101010; text-decoration: none; font-size: 0.84rem; font-weight: 500;
        padding: 6px 12px; border-radius: 999px;
    }
    .navlinks a.on {background: #101010; color: #FFFFFF;}
    .navstatus {color: #6F6F6A; font-size: 0.75rem; white-space: nowrap;}
    .page-title {
        font-size: 1.85rem; font-weight: 600; letter-spacing: -0.045em; line-height: 1.05;
        color: #101010; margin: 0 0 0.3rem;
    }
    .lede {color: #6F6F6A; font-size: 0.92rem; line-height: 1.45; margin: 0 0 0.9rem; max-width: 44rem;}
    .metrics {
        display: grid; border-top: 1px solid #E6E6E3; border-bottom: 1px solid #E6E6E3;
        margin: 0.2rem 0 1.15rem;
    }
    .metrics-4 {grid-template-columns: repeat(4, 1fr);}
    .metrics-3 {grid-template-columns: repeat(3, 1fr);}
    .metrics-2 {grid-template-columns: repeat(2, 1fr);}
    .metrics div {padding: 12px 0 14px;}
    .metrics div + div {border-left: 1px solid #E6E6E3; padding-left: 18px;}
    .metrics span {display: block; color: #6F6F6A; font-size: 0.75rem; font-weight: 500;}
    .metrics strong {
        display: block; margin-top: 2px; color: #101010; font-size: 1.7rem; font-weight: 600;
        letter-spacing: -0.04em; font-variant-numeric: tabular-nums lining-nums; line-height: 1.1;
    }
    .filter-note {color: #6F6F6A; font-size: 0.8rem; margin: 0.15rem 0 0.85rem;}
    .filter-note.on {color: #101010;}
    .key {color: #6F6F6A; font-size: 0.8rem; margin: 0 0 0.35rem;}
    .swatch {display: inline-block; width: 8px; height: 8px; border-radius: 2px; margin: 0 6px 1px 12px;}
    .key .swatch:first-child {margin-left: 0;}
    .note {
        border: 1px solid #E6E6E3; background: #FFFFFF; border-radius: 12px;
        padding: 14px 16px; color: #101010; font-size: 0.92rem; line-height: 1.55; margin: 0.15rem 0 0.4rem;
    }
    .empty {
        color: #6F6F6A; border: 1px solid #E6E6E3; border-radius: 12px; background: #FFFFFF;
        padding: 1.25rem 1rem; text-align: center;
    }
    .criteria {display: flex; flex-wrap: wrap; gap: 8px 18px; margin: 0 0 0.8rem; color: #6F6F6A; font-size: 0.8rem;}
    .criteria b {color: #101010; font-weight: 600;}
    table.detail {width: 100%; border-collapse: collapse; font-size: 0.86rem; font-variant-numeric: tabular-nums;}
    table.detail th {text-align: left; font-weight: 500; color: #6F6F6A; font-size: 0.75rem; padding: 0 0 6px;}
    table.detail td {border-top: 1px solid #EFEFEC; padding: 7px 0; color: #101010;}
    table.detail td.num, table.detail th.num {text-align: right;}
    @media (max-width: 800px) {
        .navpill {border-radius: 18px; align-items: flex-start; flex-direction: column;}
        .metrics-4, .metrics-3 {grid-template-columns: 1fr 1fr;}
        .metrics div:nth-child(odd) {border-left: 0; padding-left: 0;}
    }
    </style>
    """,
    unsafe_allow_html=True,
)

if not DATABASE.exists():
    with st.spinner("Building the SQLite database from cell-count.csv"):
        try:
            load()
        except Exception as exc:
            st.error(f"Database unavailable. Run `make pipeline` first. Details: {exc}")
            st.stop()

revision = DATABASE.stat().st_mtime_ns


@st.cache_data(show_spinner="Loading clinical data")
def sql_cached(query: str, params: tuple, db_revision: int) -> pd.DataFrame:
    _ = db_revision
    return read_frame(query, params)


@st.cache_data(show_spinner="Preparing responder statistics")
def responder_data(db_revision: int) -> pd.DataFrame:
    _ = db_revision
    return cohort_frequencies()


@st.cache_data(show_spinner="Filtering frequencies")
def cached_frequencies(db_revision: int, filter_items: tuple, search: str) -> pd.DataFrame:
    _ = db_revision
    return filtered_frequencies(dict(filter_items), search)


def population_name(key: str) -> str:
    label = POPULATION_LABELS.get(key, str(key))
    return label[:1].upper() + label[1:]


def download_csv(label: str, data: pd.DataFrame, filename: str) -> None:
    st.download_button(label, data.to_csv(index=False).encode("utf-8"), filename, "text/csv")


def page_header(title: str, note: str) -> None:
    st.markdown(
        f"<h1 class='page-title'>{html.escape(title)}</h1><p class='lede'>{html.escape(note)}</p>",
        unsafe_allow_html=True,
    )


def metric_strip(items: list[tuple[str, str]]) -> None:
    cells = "".join(
        f"<div><span>{html.escape(label)}</span><strong>{html.escape(value)}</strong></div>"
        for label, value in items
    )
    st.markdown(f"<div class='metrics metrics-{len(items)}'>{cells}</div>", unsafe_allow_html=True)


def _rgba(hex_color: str, alpha: float) -> str:
    value = hex_color.lstrip("#")
    red, green, blue = int(value[0:2], 16), int(value[2:4], 16), int(value[4:6], 16)
    return f"rgba({red},{green},{blue},{alpha})"


def apply_chart(fig, height: int | None = None, margin_top: int = 36, grid: str = "y"):
    fig.update_layout(
        template="simple_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family=FONT, color=INK, size=12),
        margin=dict(l=4, r=8, t=margin_top, b=4),
        hoverlabel=dict(bgcolor="#FFFFFF", bordercolor=LINE, font=dict(family=FONT, size=12, color=INK)),
        title=dict(font=dict(family=FONT, size=13, color=INK), x=0, xanchor="left"),
        legend=dict(orientation="h", yanchor="top", y=-0.18, x=0, font=dict(size=11, color=MUTED)),
    )
    axis = dict(tickfont=dict(family=FONT, size=11, color=MUTED), title_font=dict(family=FONT, size=12, color=MUTED))
    fig.update_xaxes(showgrid=grid == "x", gridcolor=GRID, gridwidth=1, zeroline=False, showline=False, **axis)
    fig.update_yaxes(showgrid=grid == "y", gridcolor=GRID, gridwidth=1, zeroline=False, showline=False, **axis)
    if height is not None:
        fig.update_layout(height=height)
    return fig


def show_chart(fig, key: str | None = None, on_select="ignore", alt: str | None = None):
    kwargs = {"width": "stretch", "config": PLOTLY, "alt": alt}
    if key is not None:
        kwargs["key"] = key
    if on_select != "ignore":
        kwargs["on_select"] = on_select
        kwargs["selection_mode"] = "points"
    return st.plotly_chart(fig, **kwargs)


def _point_value(point) -> str | None:
    getter = point.get if isinstance(point, dict) else lambda key, default=None: getattr(point, key, default)
    data = getter("customdata")
    if isinstance(data, (list, tuple)):
        data = data[0] if data else None
    if data not in (None, ""):
        return str(data)
    for key in ("y", "x"):
        value = getter(key)
        if value not in (None, ""):
            return str(value)
    return None


def choose_filter(widget_key: str, chart_key: str, allowed: set[str]):
    def _callback() -> None:
        raw = st.session_state.get(chart_key)
        selection = raw.get("selection", {}) if isinstance(raw, dict) else getattr(raw, "selection", {})
        points = selection.get("points", []) if isinstance(selection, dict) else getattr(selection, "points", [])
        if not points:
            return
        value = _point_value(points[0])
        if value not in allowed:
            return
        current = st.session_state.get(widget_key, "All")
        st.session_state[widget_key] = "All" if current == value else value
    return _callback


def hbar(frame: pd.DataFrame, label: str, value: str, title: str, selected: str, height: int = 210):
    labels = frame[label].astype(str)
    colors = [
        INK if selected in ("All", None) or item == selected else "#D4D4D0"
        for item in labels
    ]
    fig = go.Figure(go.Bar(
        y=labels,
        x=frame[value],
        orientation="h",
        marker=dict(color=colors, line=dict(width=0)),
        text=frame[value].map(lambda item: f"{item:,}"),
        textposition="outside",
        textfont=dict(family=FONT, size=11, color=MUTED),
        customdata=labels,
        cliponaxis=False,
        hovertemplate="%{y}<br>%{x:,}<extra></extra>",
    ))
    peak = float(frame[value].max()) if len(frame) else 1.0
    fig.update_layout(title=dict(text=title), bargap=0.45, showlegend=False)
    apply_chart(fig, height=height, margin_top=40, grid="x")
    fig.update_xaxes(title_text=None, range=[0, peak * 1.45 if peak else 1])
    fig.update_yaxes(
        title_text=None,
        categoryorder="array",
        categoryarray=list(reversed(labels.tolist())),
    )
    return fig


def visit_figure(frame: pd.DataFrame, title: str):
    fig = go.Figure(go.Scatter(
        x=frame["visit"],
        y=frame["n"],
        mode="lines+markers",
        line=dict(color=INK, width=1.6),
        marker=dict(size=8, color=INK),
        hovertemplate="%{x}<br>%{y:,}<extra></extra>",
    ))
    peak = float(frame["n"].max()) if len(frame) else 1.0
    fig.update_layout(title=dict(text=title))
    apply_chart(fig, height=280, margin_top=44)
    fig.update_yaxes(title_text=None, range=[0, peak * 1.28 if peak else 1])
    fig.update_xaxes(title_text=None)
    return fig


def composition_figure(frame: pd.DataFrame, title: str, height: int = 168):
    fig = go.Figure()
    for name in POPULATIONS:
        match = frame.loc[frame["population"] == name]
        share = float(match["percentage"].iloc[0]) if not match.empty else 0.0
        fig.add_trace(go.Bar(
            y=[title],
            x=[share],
            orientation="h",
            name=population_name(name),
            marker=dict(color=POP_COLORS[name], line=dict(width=0)),
            hovertemplate=f"{population_name(name)}<br>%{{x:.2f}}%<extra></extra>",
        ))
    fig.update_layout(barmode="stack", title=dict(text=title), bargap=0.35)
    apply_chart(fig, height=height, margin_top=36, grid="x")
    fig.update_layout(
        margin=dict(l=8, r=130, t=36, b=16),
        legend=dict(orientation="v", y=0.5, yanchor="middle", x=1.02, xanchor="left", font=dict(size=11, color=MUTED)),
    )
    fig.update_xaxes(title_text=None, range=[0, 100])
    fig.update_yaxes(title_text=None, showticklabels=False)
    return fig


def box_figure(frame: pd.DataFrame, population: str, height: int = 300):
    subset = frame.loc[frame["population"] == population]
    fig = go.Figure()
    for label, color in (("Responder", RESPONDER), ("Non-responder", NONRESPONDER)):
        part = subset.loc[subset["response"] == label]
        fig.add_trace(go.Box(
            x=[label] * len(part),
            y=part["percentage"],
            name=label,
            marker=dict(color=color, size=3.5, opacity=0.45),
            line=dict(color=color, width=1.4),
            fillcolor=_rgba(color, 0.12),
            boxpoints="outliers",
            customdata=part["patient"],
            hovertemplate="%{customdata}<br>%{y:.2f}%<extra>" + label + "</extra>",
            showlegend=False,
        ))
    fig.update_layout(title=dict(text=population_name(population)), boxmode="group")
    apply_chart(fig, height=height, margin_top=40)
    fig.update_yaxes(title_text="Percent", rangemode="tozero")
    fig.update_xaxes(title_text=None, categoryorder="array", categoryarray=["Responder", "Non-responder"])
    return fig


def forest_figure(stats_table: pd.DataFrame):
    fig = go.Figure()
    fig.add_vline(x=0, line_width=1, line_color="#C8C8C4")
    for row in stats_table.itertuples(index=False):
        color = RED if row.significant else INK
        label = population_name(row.population)
        fig.add_trace(go.Scatter(
            x=[row.mean_difference],
            y=[label],
            mode="markers",
            marker=dict(size=8, color=color, line=dict(width=0)),
            error_x=dict(
                type="data",
                symmetric=False,
                array=[row.ci_upper - row.mean_difference],
                arrayminus=[row.mean_difference - row.ci_lower],
                color=color,
                thickness=1.2,
                width=5,
            ),
            hovertemplate=(
                f"{label}<br>Difference %{{x:.3f}} percentage points"
                f"<br>95% CI {row.ci_lower:.3f} to {row.ci_upper:.3f}"
                f"<br>Hedges' g {row.hedges_g:.3f}"
                f"<br>Adjusted p {row.adjusted_p_value:.4f}<extra></extra>"
            ),
            showlegend=False,
        ))
    apply_chart(fig, height=320, margin_top=12, grid="x")
    lower = float(stats_table["ci_lower"].min())
    upper = float(stats_table["ci_upper"].max())
    pad = (upper - lower) * 0.12
    fig.update_xaxes(range=[lower - pad, upper + pad], title_text="Responder − non-responder (percentage points)")
    fig.update_yaxes(
        title_text=None,
        categoryorder="array",
        categoryarray=[population_name(name) for name in reversed(POPULATIONS)],
    )
    return fig


def distinct_values(query: str) -> list:
    return sql_cached(query, (), revision)["value"].dropna().tolist()


def response_options() -> tuple[list[str], dict]:
    values = sql_cached(
        "SELECT DISTINCT response AS value FROM subjects ORDER BY value", (), revision
    )["value"].tolist()
    labels = ["All"]
    mapping: dict = {}
    for value in values:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            label, stored = "Unrecorded", "unrecorded"
        else:
            label, stored = str(value), value
        if label not in mapping and label != "All":
            labels.append(label)
            mapping[label] = stored
    return labels, mapping


def _subject_on(project: str, treatment: str, condition: str) -> tuple[str, list]:
    parts, params = [], []
    if project != "All":
        parts.append("su.project_id = ?")
        params.append(project)
    if treatment != "All":
        parts.append("su.treatment = ?")
        params.append(treatment)
    if condition != "All":
        parts.append("su.condition = ?")
        params.append(condition)
    suffix = (" AND " + " AND ".join(parts)) if parts else ""
    return suffix, params


def overview_totals(project: str, treatment: str, condition: str) -> dict:
    if project == "All" and treatment == "All" and condition == "All":
        row = sql_cached(
            """
            SELECT (SELECT COUNT(*) FROM samples) AS samples,
                   (SELECT COUNT(*) FROM subjects) AS patients,
                   (SELECT COUNT(*) FROM projects) AS projects,
                   (SELECT COUNT(*) FROM cell_counts) AS cell_counts
            """,
            (), revision,
        ).iloc[0]
        return {key: int(row[key]) for key in ("samples", "patients", "projects", "cell_counts")}
    suffix, params = _subject_on(project, treatment, condition)
    where = f"WHERE 1=1{suffix}" if suffix else ""
    # suffix is " AND col = ?" so attach it to a real predicate.
    where = ("WHERE " + suffix.removeprefix(" AND ")) if suffix else ""
    counts = sql_cached(
        f"""
        SELECT COUNT(DISTINCT s.sample_id) AS samples,
               COUNT(DISTINCT su.subject_id) AS patients,
               COUNT(DISTINCT su.project_id) AS projects
        FROM samples s JOIN subjects su ON su.subject_id = s.subject_id
        {where}
        """,
        tuple(params), revision,
    ).iloc[0]
    cells = sql_cached(
        f"""
        SELECT COUNT(*) AS cell_counts
        FROM cell_counts cc
        JOIN samples s ON s.sample_id = cc.sample_id
        JOIN subjects su ON su.subject_id = s.subject_id
        {where}
        """,
        tuple(params), revision,
    ).iloc[0]
    return {
        "samples": int(counts["samples"]),
        "patients": int(counts["patients"]),
        "projects": int(counts["projects"]),
        "cell_counts": int(cells["cell_counts"]),
    }


def visit_counts(project: str, treatment: str, condition: str, unit: str) -> pd.DataFrame:
    suffix, params = _subject_on(project, treatment, condition)
    aggregate = "COUNT(DISTINCT su.subject_id)" if unit == "Patients" else "COUNT(s.sample_id)"
    frame = sql_cached(
        f"""
        SELECT d.day AS day, {aggregate} AS n
        FROM (SELECT DISTINCT time_from_treatment_start AS day FROM samples) d
        LEFT JOIN samples s ON s.time_from_treatment_start = d.day
        LEFT JOIN subjects su ON su.subject_id = s.subject_id{suffix}
        GROUP BY d.day
        ORDER BY d.day
        """,
        tuple(params), revision,
    )
    frame["visit"] = frame["day"].map(lambda day: f"Day {int(day)}")
    return frame


def breakdown(dimension: str, project: str, treatment: str, condition: str, unit: str) -> pd.DataFrame:
    """One chart's categories stay visible; the other two filters apply."""
    filters = {"project": project, "treatment": treatment, "condition": condition}
    filters[dimension] = "All"
    suffix, params = _subject_on(filters["project"], filters["treatment"], filters["condition"])
    aggregate = "COUNT(DISTINCT su.subject_id)" if unit == "Patients" else "COUNT(s.sample_id)"
    if dimension == "project":
        query = f"""
            SELECT p.project_id AS label, {aggregate} AS n
            FROM projects p
            LEFT JOIN subjects su ON su.project_id = p.project_id{suffix}
            LEFT JOIN samples s ON s.subject_id = su.subject_id
            GROUP BY p.project_id
            ORDER BY p.project_id
        """
    elif dimension == "treatment":
        query = f"""
            SELECT t.treatment AS label, {aggregate} AS n
            FROM (SELECT DISTINCT treatment FROM subjects) t
            LEFT JOIN subjects su ON su.treatment = t.treatment{suffix}
            LEFT JOIN samples s ON s.subject_id = su.subject_id
            GROUP BY t.treatment
            ORDER BY n DESC, t.treatment
        """
    else:
        query = f"""
            SELECT c.condition AS label, {aggregate} AS n
            FROM (SELECT DISTINCT condition FROM subjects) c
            LEFT JOIN subjects su ON su.condition = c.condition{suffix}
            LEFT JOIN samples s ON s.subject_id = su.subject_id
            GROUP BY c.condition
            ORDER BY n DESC, c.condition
        """
    return sql_cached(query, tuple(params), revision)


def pooled_composition(project: str, treatment: str, condition: str) -> pd.DataFrame:
    clauses, params = [], []
    if project != "All":
        clauses.append("project = ?")
        params.append(project)
    if treatment != "All":
        clauses.append("treatment = ?")
        params.append(treatment)
    if condition != "All":
        clauses.append("condition = ?")
        params.append(condition)
    where = ("WHERE " + " AND ".join(clauses)) if clauses else ""
    frame = sql_cached(
        f"SELECT population, SUM(count) AS count FROM frequency_details {where} GROUP BY population",
        tuple(params), revision,
    )
    total = float(frame["count"].sum()) if not frame.empty else 0.0
    frame["percentage"] = frame["count"] / total * 100 if total else 0.0
    return frame


def sample_rows(sample_id: str) -> pd.DataFrame:
    frame = sql_cached(
        """
        SELECT sample, total_count, population, count, percentage
        FROM frequency_details WHERE sample = ?
        """,
        (sample_id,), revision,
    )
    frame["population"] = pd.Categorical(frame["population"], categories=POPULATIONS, ordered=True)
    return frame.sort_values("population")


def legend() -> None:
    st.markdown(
        f"<p class='key'><span class='swatch' style='background:{RESPONDER}'></span>Responder"
        f"<span class='swatch' style='background:{NONRESPONDER}'></span>Non-responder</p>",
        unsafe_allow_html=True,
    )


def reset_overview() -> None:
    st.session_state["ov_project"] = "All"
    st.session_state["ov_treatment"] = "All"
    st.session_state["ov_condition"] = "All"
    st.session_state["ov_unit"] = "Samples"


def show_overview() -> None:
    page_header(
        "Trial Overview",
        "Structure of the loaded trial. Filters here are exploratory and do not change the response or subset cohorts.",
    )
    projects = distinct_values("SELECT project_id AS value FROM projects ORDER BY 1")
    treatments = distinct_values("SELECT DISTINCT treatment AS value FROM subjects ORDER BY 1")
    conditions = distinct_values("SELECT DISTINCT condition AS value FROM subjects ORDER BY 1")
    filters = st.columns([1.1, 1.1, 1.1, 1.2, 0.7], gap="small")
    with filters[0]:
        project = st.selectbox("Project", ["All", *projects], key="ov_project")
    with filters[1]:
        treatment = st.selectbox("Treatment", ["All", *treatments], key="ov_treatment")
    with filters[2]:
        condition = st.selectbox("Condition", ["All", *conditions], key="ov_condition")
    with filters[3]:
        unit = st.segmented_control(
            "Count",
            ["Samples", "Patients"],
            default="Samples",
            key="ov_unit",
            required=True,
            width="content",
        )
    with filters[4]:
        st.markdown("<div style='height:1.7rem'></div>", unsafe_allow_html=True)
        st.button("Reset", on_click=reset_overview, key="ov_reset")
    if unit not in ("Samples", "Patients"):
        unit = "Samples"
    active = [f"{label} {value}" for label, value in (
        ("project", project), ("treatment", treatment), ("condition", condition),
    ) if value != "All"]
    if active:
        st.markdown(
            "<p class='filter-note on'><span class='dot'></span> "
            + html.escape(" · ".join(active))
            + ". Response analysis and subset analysis are unchanged.</p>",
            unsafe_allow_html=True,
        )
    totals = overview_totals(project, treatment, condition)
    metric_strip([
        ("Samples", f"{totals['samples']:,}"),
        ("Patients", f"{totals['patients']:,}"),
        ("Projects", f"{totals['projects']:,}"),
        ("Cell measurements", f"{totals['cell_counts']:,}"),
    ])

    visits = visit_counts(project, treatment, condition, unit)
    composition = pooled_composition(project, treatment, condition)
    left, right = st.columns([1.35, 1], gap="small")
    with left:
        show_chart(
            visit_figure(visits, f"{unit} by visit"),
            alt=f"{unit} on day 0, day 7, and day 14",
        )
    with right:
        if composition.empty or float(composition["count"].sum()) == 0:
            st.markdown("<p class='empty'>No cell counts in this filter.</p>", unsafe_allow_html=True)
        else:
            show_chart(
                composition_figure(composition, "Pooled share"),
                alt="Pooled share of the five populations in the current filter",
            )
            st.caption("Pooled from counts in the current filter, not a mean of percentages.")

    project_frame = breakdown("project", project, treatment, condition, unit)
    treatment_frame = breakdown("treatment", project, treatment, condition, unit)
    condition_frame = breakdown("condition", project, treatment, condition, unit)
    panels = st.columns(3, gap="small")
    charts = (
        (panels[0], project_frame, "By project", "ov_project", "ov_project_chart", project),
        (panels[1], treatment_frame, "By treatment", "ov_treatment", "ov_treatment_chart", treatment),
        (panels[2], condition_frame, "By condition", "ov_condition", "ov_condition_chart", condition),
    )
    for panel, frame, title, widget_key, chart_key, selected in charts:
        with panel:
            show_chart(
                hbar(frame, "label", "n", title, selected),
                key=chart_key,
                on_select=choose_filter(widget_key, chart_key, set(frame["label"].astype(str))),
                alt=title,
            )
    st.caption("Select a bar to filter the other views. Select it again to clear that filter.")


def reset_frequencies() -> None:
    for key, value in (
        ("freq_search", ""),
        ("freq_project", "All"),
        ("freq_condition", "All"),
        ("freq_treatment", "All"),
        ("freq_response", "All"),
        ("freq_sex", "All"),
        ("freq_type", "All"),
        ("freq_day", "All"),
        ("freq_population", "All"),
        ("freq_page", 0),
    ):
        st.session_state[key] = value


def show_frequencies() -> None:
    page_header(
        "Cell Frequencies",
        "One row is one population in one sample. Percentage uses that sample's five populations as the denominator.",
    )
    projects = distinct_values("SELECT project_id AS value FROM projects ORDER BY 1")
    conditions = distinct_values("SELECT DISTINCT condition AS value FROM subjects ORDER BY 1")
    treatments = distinct_values("SELECT DISTINCT treatment AS value FROM subjects ORDER BY 1")
    sexes = distinct_values("SELECT DISTINCT sex AS value FROM subjects ORDER BY 1")
    sample_types = distinct_values("SELECT DISTINCT sample_type AS value FROM samples ORDER BY 1")
    days = distinct_values("SELECT DISTINCT time_from_treatment_start AS value FROM samples ORDER BY 1")
    response_labels, response_map = response_options()

    search = st.text_input("Sample or patient", key="freq_search", placeholder="Search a sample id or subject id")
    row1 = st.columns(4, gap="small")
    row2 = st.columns(5, gap="small")
    with row1[0]:
        project = st.selectbox("Project", ["All", *projects], key="freq_project")
    with row1[1]:
        condition = st.selectbox("Condition", ["All", *conditions], key="freq_condition")
    with row1[2]:
        treatment = st.selectbox("Treatment", ["All", *treatments], key="freq_treatment")
    with row1[3]:
        response_label = st.selectbox("Response", response_labels, key="freq_response")
    with row2[0]:
        sex = st.selectbox("Sex", ["All", *sexes], key="freq_sex", format_func=lambda value: SEX_LABELS.get(value, value))
    with row2[1]:
        sample_type = st.selectbox("Sample type", ["All", *sample_types], key="freq_type")
    with row2[2]:
        day = st.selectbox(
            "Day", ["All", *days], key="freq_day",
            format_func=lambda value: "All" if value == "All" else f"Day {int(value)}",
        )
    with row2[3]:
        population = st.selectbox(
            "Population", ["All", *POPULATIONS], key="freq_population",
            format_func=lambda value: "All" if value == "All" else population_name(value),
        )
    with row2[4]:
        st.markdown("<div style='height:1.7rem'></div>", unsafe_allow_html=True)
        st.button("Reset", on_click=reset_frequencies, key="freq_reset")

    active: dict = {}
    active_labels = []
    for key, chosen, label in (
        ("project", project, "project"),
        ("condition", condition, "condition"),
        ("treatment", treatment, "treatment"),
        ("sex", sex, "sex"),
        ("sample_type", sample_type, "sample type"),
        ("population", population, "population"),
    ):
        if chosen != "All":
            active[key] = chosen
            shown = population_name(chosen) if key == "population" else SEX_LABELS.get(chosen, chosen)
            active_labels.append(f"{label} {shown}")
    if day != "All":
        active["time_from_treatment_start"] = int(day)
        active_labels.append(f"day {int(day)}")
    if response_label != "All":
        active["response"] = response_map[response_label]
        active_labels.append(f"response {response_label}")
    if search.strip():
        active_labels.append(f"search {search.strip()}")
    if active_labels:
        st.markdown(
            "<p class='filter-note on'><span class='dot'></span>"
            + html.escape(" · ".join(active_labels)) + "</p>",
            unsafe_allow_html=True,
        )

    data = cached_frequencies(revision, tuple(sorted(active.items(), key=lambda item: item[0])), search)
    signature = tuple(sorted(active.items())) + (search.strip(),)
    if st.session_state.get("freq_sig") != signature:
        st.session_state["freq_sig"] = signature
        st.session_state["freq_page"] = 0
    sample_count = int(data["sample"].nunique()) if not data.empty else 0
    median = f"{data['total_count'].median():,.0f}" if not data.empty else "—"
    st.markdown(
        f"<p class='filter-note'>{len(data):,} rows · {sample_count:,} samples · median total {median}</p>",
        unsafe_allow_html=True,
    )
    if data.empty:
        st.markdown("<p class='empty'>Nothing matches. Reset the filters or clear the search.</p>", unsafe_allow_html=True)
        return

    pages = max(1, (len(data) + PAGE_SIZE - 1) // PAGE_SIZE)
    page_index = min(int(st.session_state.get("freq_page", 0)), pages - 1)
    st.session_state["freq_page"] = page_index
    start = page_index * PAGE_SIZE
    view = data.iloc[start:start + PAGE_SIZE]
    event = st.dataframe(
        view,
        width="stretch",
        hide_index=True,
        height=460,
        on_select="rerun",
        selection_mode="single-row",
        key=f"freq_grid_{page_index}_{start}",
        column_config={
            "sample": "sample",
            "total_count": st.column_config.NumberColumn("total_count", format="%d"),
            "population": "population",
            "count": st.column_config.NumberColumn("count", format="%d"),
            "percentage": st.column_config.NumberColumn("percentage", format="%.2f"),
        },
    )
    pager = st.columns([1, 1, 3], gap="small")
    with pager[0]:
        st.button("Previous", disabled=page_index == 0, on_click=lambda: st.session_state.__setitem__("freq_page", page_index - 1))
    with pager[1]:
        st.button("Next", disabled=page_index >= pages - 1, on_click=lambda: st.session_state.__setitem__("freq_page", page_index + 1))
    with pager[2]:
        st.caption(f"{start + 1:,}–{start + len(view):,} of {len(data):,}")

    selected_rows = []
    if event is not None:
        selection = event["selection"] if isinstance(event, dict) else event.selection
        selected_rows = list(selection["rows"] if isinstance(selection, dict) else selection.rows)
    sample_id = None
    if selected_rows:
        sample_id = str(view.iloc[int(selected_rows[0])]["sample"])
    elif sample_count == 1:
        sample_id = str(data["sample"].iloc[0])

    detail = sample_rows(sample_id) if sample_id else None
    chart_frame = detail if detail is not None and not detail.empty else None
    if chart_frame is None:
        grouped = data.groupby("population", as_index=False)["count"].sum()
        total = float(grouped["count"].sum())
        grouped["percentage"] = grouped["count"] / total * 100 if total else 0.0
        chart_title = "Pooled share of these rows"
        chart_note = "Pooled from the filtered counts. Select a row to inspect one sample."
    else:
        grouped = chart_frame
        chart_title = str(sample_id)
        chart_note = "All five populations in this sample. Percentages sum to 100."
    chart_col, table_col = st.columns([1.3, 1], gap="medium")
    with chart_col:
        show_chart(composition_figure(grouped, chart_title, height=210), alt=chart_note)
        st.caption(chart_note)
    with table_col:
        if chart_frame is not None:
            body = "".join(
                "<tr>"
                f"<td>{html.escape(population_name(str(row.population)))}</td>"
                f"<td class='num'>{int(row.count):,}</td>"
                f"<td class='num'>{row.percentage:.2f}</td>"
                "</tr>"
                for row in chart_frame.itertuples(index=False)
            )
            total_count = int(chart_frame["total_count"].iloc[0])
            st.markdown(
                f"<p class='filter-note'>Sample {html.escape(sample_id)} · total {total_count:,}</p>"
                "<table class='detail'><thead><tr><th>Population</th><th class='num'>Count</th>"
                f"<th class='num'>Percent</th></tr></thead><tbody>{body}</tbody></table>",
                unsafe_allow_html=True,
            )
        else:
            st.caption("Select a table row for the sample's own counts and frequencies.")
    action, _rest = st.columns([1.3, 4])
    with action:
        download_csv("Download CSV", data, "sample_frequencies_filtered.csv")
    st.caption("The file keeps the five columns unrounded, for every matching row, not only this page.")


def results_table(stats_table: pd.DataFrame) -> None:
    shown = pd.DataFrame({
        "Population": stats_table["population"].map(population_name),
        "Responders": stats_table["responder_n"].astype(int),
        "Non-responders": stats_table["non_responder_n"].astype(int),
        "Difference": stats_table["mean_difference"],
        "CI lower": stats_table["ci_lower"],
        "CI upper": stats_table["ci_upper"],
        "Hedges' g": stats_table["hedges_g"],
        "Cohen's d": stats_table["cohens_d"],
        "p": stats_table["p_value"],
        "Adjusted p": stats_table["adjusted_p_value"],
        "Significant": stats_table["significant"].map(lambda flag: "Yes" if flag else "No"),
    })

    def paint(row: pd.Series) -> list[str]:
        styles = [""] * len(row)
        if row["Significant"] == "Yes":
            styles[-1] = "color: #FF3638; font-weight: 600"
        return styles

    styled = (
        shown.style
        .format({
            "Difference": "{:+.3f}",
            "CI lower": "{:.3f}",
            "CI upper": "{:.3f}",
            "Hedges' g": "{:+.3f}",
            "Cohen's d": "{:+.3f}",
            "p": "{:.4f}",
            "Adjusted p": "{:.4f}",
        })
        .apply(paint, axis=1)
        .hide(axis="index")
    )
    st.dataframe(styled, width="stretch", hide_index=True)


def show_response() -> None:
    page_header(
        "Response Analysis",
        "Melanoma, miraclib, PBMC, response yes or no. One value per patient. The timepoint changes the comparison, not who is included.",
    )
    current_label = st.segmented_control(
        "Timepoint",
        options=list(TIMEPOINTS),
        default=list(TIMEPOINTS)[0],
        required=True,
        key="resp_time",
        format_func=lambda key: TIMEPOINT_SHORT.get(key, key),
        help="All visits uses each patient's mean percentage. A single day uses that visit only.",
        width="content",
    )
    if current_label not in TIMEPOINTS:
        current_label = list(TIMEPOINTS)[0]
    focus = st.segmented_control(
        "Population",
        options=["All", *POPULATIONS],
        default="All",
        key="resp_focus",
        format_func=lambda value: "All five" if value == "All" else population_name(value),
        width="content",
    )
    if focus not in ("All", *POPULATIONS):
        focus = "All"
    selected_day = TIMEPOINTS[current_label]
    cohort = responder_data(revision)
    values = independent_values(cohort, selected_day)
    stats_table = compare(cohort, selected_day)
    patients = values[["subject_id", "response"]].drop_duplicates()
    responders = int((patients["response"] == "yes").sum())
    non_responders = int((patients["response"] == "no").sum())
    significant_count = int(stats_table["significant"].sum())
    metric_strip([
        ("Responders", f"{responders:,}"),
        ("Non-responders", f"{non_responders:,}"),
        ("Subjects", f"{responders + non_responders:,}"),
        ("Significant populations", f"{significant_count} of 5"),
    ])

    plot_values = values.rename(columns={"subject_id": "patient"}).copy()
    plot_values["response"] = plot_values["response"].map({"yes": "Responder", "no": "Non-responder"})
    legend()
    if focus == "All":
        top, bottom = st.columns(3, gap="small"), st.columns(3, gap="small")
        for slot, population in zip([*top, *bottom], POPULATIONS):
            with slot:
                show_chart(box_figure(plot_values, population, 280), alt=population_name(population))
    else:
        show_chart(box_figure(plot_values, focus, 420), alt=population_name(focus))
        others = [name for name in POPULATIONS if name != focus]
        slots = st.columns(4, gap="small")
        for slot, population in zip(slots, others):
            with slot:
                show_chart(box_figure(plot_values, population, 230), alt=population_name(population))
    st.caption("Axes start at zero and are scaled independently. Points are outliers beyond 1.5 IQR.")

    st.markdown("<h3>Tests</h3>", unsafe_allow_html=True)
    st.caption(
        "Difference is responder mean minus non-responder mean, in percentage points. "
        "Cohen's d and Hedges' g use the pooled standard deviation. "
        "Adjustment is only across these five populations."
    )
    results_table(stats_table)
    downloads = st.columns([1.25, 1.55, 3.2], gap="small")
    with downloads[0]:
        download_csv("Download statistics", stats_table, "responder_comparison.csv")
    with downloads[1]:
        download_csv("Download patient values", values, "patient_level_frequencies.csv")

    st.markdown("<h3>Confidence intervals</h3>", unsafe_allow_html=True)
    show_chart(forest_figure(stats_table), alt="Mean differences and 95 percent confidence intervals")
    st.caption("Positive values are higher in responders. Red marks adjusted p < 0.05.")

    st.markdown("<h3>Summary</h3>", unsafe_allow_html=True)
    st.markdown(
        f"<div class='note'>{html.escape(written_conclusion(stats_table, current_label))}</div>",
        unsafe_allow_html=True,
    )


def show_subset() -> None:
    page_header(
        "Subset Analysis",
        "Two fixed cohorts. The baseline block is not the cohort used for the B-cell mean.",
    )
    base = baseline_samples()
    projects = baseline_project_counts()
    responses = baseline_breakdown("response")
    sexes = baseline_breakdown("sex")

    st.markdown("<h3>Baseline cohort</h3>", unsafe_allow_html=True)
    st.markdown(
        "<p class='criteria'><span><b>Condition</b> melanoma</span><span><b>Treatment</b> miraclib</span>"
        "<span><b>Sample</b> PBMC</span><span><b>Day</b> 0</span></p>",
        unsafe_allow_html=True,
    )
    metric_strip([
        ("Eligible samples", f"{len(base):,}"),
        ("Patients", f"{base['subject'].nunique():,}"),
    ])
    response_plot = responses.copy()
    response_plot["response"] = response_plot["response"].map(
        {"yes": "Responder", "no": "Non-responder"}
    ).fillna(response_plot["response"])
    sex_plot = sexes.copy()
    sex_plot["sex"] = sex_plot["sex"].map(SEX_LABELS).fillna(sex_plot["sex"])
    panels = st.columns(3, gap="small")
    charts = (
        (panels[0], projects, "project", "samples", "Samples by project"),
        (panels[1], response_plot, "response", "subjects", "Patients by response"),
        (panels[2], sex_plot, "sex", "subjects", "Patients by sex"),
    )
    for panel, frame, x_col, y_col, title in charts:
        with panel:
            show_chart(hbar(frame, x_col, y_col, title, "All", height=220), alt=title)
    zero_projects = projects.loc[projects["samples"] == 0, "project"].astype(str).tolist()
    if zero_projects:
        st.caption("No qualifying baseline samples in " + ", ".join(zero_projects) + ".")
    response_total = int(responses["subjects"].sum())
    sex_total = int(sexes["subjects"].sum())
    patient_total = int(base["subject"].nunique())
    if response_total != patient_total or sex_total != patient_total:
        st.error(
            "Distinct patient counts do not match the baseline sample list "
            f"(response {response_total}, sex {sex_total}, samples {patient_total})."
        )
    if int(projects["samples"].sum()) != len(base):
        st.error("Project sample counts do not add up to the baseline sample list.")

    st.markdown("<h3>B-cell mean</h3>", unsafe_allow_html=True)
    st.markdown(
        "<p class='criteria'><span><b>Condition</b> melanoma</span><span><b>Sex</b> male</span>"
        "<span><b>Response</b> yes</span><span><b>Day</b> 0</span>"
        "<span><b>Treatment</b> all</span><span><b>Sample type</b> all</span></p>",
        unsafe_allow_html=True,
    )
    st.caption("Absolute B-cell count. This is not limited to miraclib or PBMC, and it is not a percentage.")
    average, count = male_responder_baseline_b_mean()
    metric_strip([
        ("Mean B-cell count", f"{average:.2f}"),
        ("Records", f"{count:,}"),
    ])

    st.markdown("<h3>Baseline samples</h3>", unsafe_allow_html=True)
    query = st.text_input("Find a baseline sample", placeholder="Sample id or subject id", key="subset_search")
    shown = base
    if query.strip():
        token = query.strip().lower()
        shown = base[
            base["sample"].astype(str).str.lower().str.contains(token, regex=False)
            | base["subject"].astype(str).str.lower().str.contains(token, regex=False)
        ]
    if shown.empty:
        st.markdown("<p class='empty'>No baseline sample matches that search.</p>", unsafe_allow_html=True)
    else:
        st.caption(f"{len(shown):,} of {len(base):,} baseline samples")
        st.dataframe(shown, width="stretch", hide_index=True, height=380)
    download_cols = st.columns(4, gap="small")
    with download_cols[0]:
        download_csv("Download samples", base, "baseline_miraclib_pbmc.csv")
    with download_cols[1]:
        download_csv("Project counts", projects, "baseline_counts_project.csv")
    with download_cols[2]:
        download_csv("Response counts", responses, "baseline_counts_response.csv")
    with download_cols[3]:
        download_csv("Sex counts", sexes, "baseline_counts_sex.csv")


def shell(current: str) -> None:
    built = datetime.fromtimestamp(DATABASE.stat().st_mtime).strftime("%Y-%m-%d")
    links = []
    for path, label in NAV:
        href = "/" if path == "" else f"/{path}"
        selected = " on" if current == path else ""
        current_attr = " aria-current='page'" if current == path else ""
        links.append(f"<a class='{selected.strip()}' href='{href}'{current_attr}>{html.escape(label)}</a>")
    st.markdown(
        "<div class='navpill'>"
        "<a class='brandlink' href='/'><span class='dot'></span>Loblaw Bio Analysis</a>"
        f"<nav class='navlinks'>{''.join(links)}</nav>"
        f"<span class='navstatus'>Database {html.escape(built)}</span>"
        "</div>",
        unsafe_allow_html=True,
    )


pages = st.navigation(
    [
        st.Page(show_overview, title="Overview", url_path="overview", default=True),
        st.Page(show_frequencies, title="Cell Frequencies", url_path="frequencies"),
        st.Page(show_response, title="Response Analysis", url_path="response"),
        st.Page(show_subset, title="Subset Analysis", url_path="subset"),
    ],
    position="hidden",
)
shell(pages.url_path)
pages.run()
with st.expander("Cohort definitions"):
    st.markdown(
        "Frequencies use the five populations in each sample as the denominator.\n\n"
        "Response analysis is melanoma, miraclib, and PBMC, with response yes or no. "
        "One observation per patient. Changing the timepoint does not change that cohort.\n\n"
        "The overview filters are exploratory. They do not change the response comparison "
        "or the subset queries.\n\n"
        "Subset analysis is the response cohort at day 0.\n\n"
        "The B-cell mean is separate: melanoma, male, response yes, day 0, "
        "every treatment and sample type."
    )
