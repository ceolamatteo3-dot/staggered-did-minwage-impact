from __future__ import annotations

import glob
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Paths and configuration
# ---------------------------------------------------------------------

PROJECT_DIR =  Path('C:/Users/Federico/Desktop/causal_inference_wage_employment')

QCEW_GLOB = str(PROJECT_DIR / "data" / "qcew" / "*Food services*.csv")
MW_PATH = PROJECT_DIR / "data" / "mw" / "mw_state_quarterly.xlsx"

OUTPUT_DIR = PROJECT_DIR / "output"
TABLE_DIR = OUTPUT_DIR / "tables"
OUT_PATH = OUTPUT_DIR / "merged_qcew_mw_panel.csv"

START_YEAR = 2009
END_YEAR = 2019

# Treatment searches begin after eight quarters of pre-treatment data.
TREATMENT_START_YEAR = 2011

# Main binary estimand: first state-initiated, binding increase of >= 5%.
MAJOR_HIKE_THRESHOLD = 0.05

# Used to distinguish genuine changes from floating-point differences.
WAGE_TOLERANCE = 0.005

# CPI-U annual averages, all items, U.S. city average, 1982-84 = 100.
CPI_U = {
    2009: 214.537,
    2010: 218.056,
    2011: 224.939,
    2012: 229.594,
    2013: 232.957,
    2014: 236.736,
    2015: 237.017,
    2016: 240.007,
    2017: 245.120,
    2018: 251.107,
    2019: 255.657,
}
CPI_2012 = CPI_U[2012]

VALID_STATE_FIPS = {
    1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
    35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50,
    51, 53, 54, 55, 56,
}

CENSUS_DIVISION = {
    # New England
    9: "New England", 23: "New England", 25: "New England",
    33: "New England", 44: "New England", 50: "New England",

    # Middle Atlantic
    34: "Middle Atlantic", 36: "Middle Atlantic", 42: "Middle Atlantic",

    # East North Central
    17: "East North Central", 18: "East North Central",
    26: "East North Central", 39: "East North Central",
    55: "East North Central",

    # West North Central
    19: "West North Central", 20: "West North Central",
    27: "West North Central", 29: "West North Central",
    31: "West North Central", 38: "West North Central",
    46: "West North Central",

    # South Atlantic
    10: "South Atlantic", 11: "South Atlantic", 12: "South Atlantic",
    13: "South Atlantic", 24: "South Atlantic", 37: "South Atlantic",
    45: "South Atlantic", 51: "South Atlantic", 54: "South Atlantic",

    # East South Central
    1: "East South Central", 21: "East South Central",
    28: "East South Central", 47: "East South Central",

    # West South Central
    5: "West South Central", 22: "West South Central",
    40: "West South Central", 48: "West South Central",

    # Mountain
    4: "Mountain", 8: "Mountain", 16: "Mountain", 30: "Mountain",
    32: "Mountain", 35: "Mountain", 49: "Mountain", 56: "Mountain",

    # Pacific
    2: "Pacific", 6: "Pacific", 15: "Pacific",
    41: "Pacific", 53: "Pacific",
}


def setup_directories() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)


def _require_columns(
    df: pd.DataFrame,
    columns: list[str],
    source_name: str,
) -> None:
    missing = sorted(set(columns) - set(df.columns))
    if missing:
        raise ValueError(
            f"{source_name} is missing required columns: {missing}"
        )


def load_qcew_state_private(
    path_pattern: str,
    start_year: int,
    end_year: int,
) -> pd.DataFrame:
    """
    Load statewide, private, NAICS 722 QCEW observations.

    QCEW aggregation level 55:
        Statewide, NAICS 3-digit, by ownership sector.

    QCEW ownership code 5:
        Private.
    """
    files = sorted(
        f for f in glob.glob(path_pattern)
        if "merged" not in f.lower() and "panel" not in f.lower()
    )

    if not files:
        raise FileNotFoundError(
            f"No QCEW files found matching: {path_pattern}"
        )

    frames: list[pd.DataFrame] = []
    file_audit: list[dict] = []

    required = [
        "area_fips",
        "own_code",
        "industry_code",
        "agglvl_code",
        "year",
        "qtr",
        "disclosure_code",
        "month1_emplvl",
        "month2_emplvl",
        "month3_emplvl",
        "avg_wkly_wage",
    ]

    for file_name in files:
        path = Path(file_name)

        try:
            raw = pd.read_csv(
                path,
                dtype=str,
                keep_default_na=False,
                low_memory=False,
            )
        except Exception as exc:
            raise RuntimeError(f"Could not read {path}: {exc}") from exc

        raw.columns = raw.columns.str.strip()
        _require_columns(raw, required, str(path))

        selected = raw.loc[
            (raw["agglvl_code"].str.strip() == "55")
            & (raw["own_code"].str.strip() == "5")
            & (raw["industry_code"].str.strip() == "722")
        ].copy()

        selected["year"] = pd.to_numeric(
            selected["year"], errors="coerce"
        )
        selected["qtr"] = pd.to_numeric(
            selected["qtr"], errors="coerce"
        )

        selected = selected.loc[
            selected["year"].between(start_year, end_year)
            & selected["qtr"].between(1, 4)
        ].copy()

        before_disclosure = len(selected)

        # "N" means employment/wages are not disclosed.
        selected = selected.loc[
            selected["disclosure_code"].str.strip().ne("N")
        ].copy()

        file_audit.append(
            {
                "file": path.name,
                "selected_before_disclosure_filter": before_disclosure,
                "selected_after_disclosure_filter": len(selected),
            }
        )

        frames.append(selected)

    pd.DataFrame(file_audit).to_csv(
        TABLE_DIR / "qcew_file_audit.csv",
        index=False,
    )

    qcew = pd.concat(frames, ignore_index=True)

    if qcew.empty:
        raise ValueError(
            "No statewide private NAICS 722 observations were extracted."
        )

    numeric_columns = [
        "year",
        "qtr",
        "month1_emplvl",
        "month2_emplvl",
        "month3_emplvl",
        "avg_wkly_wage",
    ]
    for column in numeric_columns:
        qcew[column] = pd.to_numeric(qcew[column], errors="coerce")

    # Require all three monthly employment values. Do not treat suppressed
    # zero-filled observations as real zero employment.
    qcew = qcew.dropna(
        subset=[
            "year",
            "qtr",
            "month1_emplvl",
            "month2_emplvl",
            "month3_emplvl",
        ]
    ).copy()

    qcew = qcew.loc[
        (qcew["month1_emplvl"] > 0)
        & (qcew["month2_emplvl"] > 0)
        & (qcew["month3_emplvl"] > 0)
    ].copy()

    qcew["area_fips"] = (
        qcew["area_fips"].astype(str).str.strip().str.zfill(5)
    )
    qcew["state_fips"] = pd.to_numeric(
        qcew["area_fips"].str[:2],
        errors="coerce",
    )

    qcew = qcew.loc[
        qcew["state_fips"].isin(VALID_STATE_FIPS)
    ].copy()

    qcew["state_fips"] = qcew["state_fips"].astype(int)
    qcew["year"] = qcew["year"].astype(int)
    qcew["qtr"] = qcew["qtr"].astype(int)

    duplicate_mask = qcew.duplicated(
        subset=["state_fips", "year", "qtr"],
        keep=False,
    )
    if duplicate_mask.any():
        duplicates = qcew.loc[
            duplicate_mask,
            ["state_fips", "year", "qtr", "area_title"],
        ].sort_values(["state_fips", "year", "qtr"])

        duplicates.to_csv(
            TABLE_DIR / "qcew_duplicate_state_quarters.csv",
            index=False,
        )
        raise ValueError(
            "Duplicate QCEW state-quarter observations found. "
            "See output/tables/qcew_duplicate_state_quarters.csv."
        )

    years_present = sorted(qcew["year"].unique().tolist())
    expected_years = list(range(start_year, end_year + 1))
    missing_years = sorted(set(expected_years) - set(years_present))

    if missing_years:
        raise ValueError(
            "The QCEW input is incomplete. Missing years: "
            f"{missing_years}. The analysis requires {start_year}-{end_year}."
        )

    qcew["quarterly_date"] = (
        qcew["year"].astype(str)
        + "q"
        + qcew["qtr"].astype(str)
    )

    qcew["avg_emplvl"] = qcew[
        ["month1_emplvl", "month2_emplvl", "month3_emplvl"]
    ].mean(axis=1)

    qcew["log_emp"] = np.log(qcew["avg_emplvl"])

    keep_columns = [
        "state_fips",
        "year",
        "qtr",
        "quarterly_date",
        "area_title",
        "avg_emplvl",
        "log_emp",
        "avg_wkly_wage",
    ]
    return qcew[keep_columns].copy()


def load_minimum_wage_data(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Minimum-wage file not found: {path}")

    mw = pd.read_excel(path)
    mw.columns = mw.columns.str.strip()

    rename_map = {
        "State FIPS Code": "state_fips",
        "Name": "state_name",
        "State Abbreviation": "state_abbr",
        "Quarterly Date": "quarterly_date",
        "Quarterly Federal Minimum": "federal_mw_start",
        "Quarterly State Minimum": "state_mw_start",
        "Quarterly Federal Average": "federal_mw_avg",
        "Quarterly State Average": "state_mw_avg",
        "Quarterly Federal Maximum": "federal_mw_end",
        "Quarterly State Maximum": "state_mw_end",
    }
    mw = mw.rename(columns=rename_map)

    required = [
        "state_fips",
        "state_name",
        "state_abbr",
        "quarterly_date",
        "federal_mw_avg",
        "state_mw_avg",
    ]
    _require_columns(mw, required, str(path))

    mw["state_fips"] = pd.to_numeric(
        mw["state_fips"], errors="coerce"
    )
    mw = mw.loc[mw["state_fips"].isin(VALID_STATE_FIPS)].copy()
    mw["state_fips"] = mw["state_fips"].astype(int)

    mw["quarterly_date"] = (
        mw["quarterly_date"].astype(str).str.strip().str.lower()
    )

    wage_columns = [
        "federal_mw_start",
        "state_mw_start",
        "federal_mw_avg",
        "state_mw_avg",
        "federal_mw_end",
        "state_mw_end",
    ]
    for column in wage_columns:
        if column in mw.columns:
            mw[column] = pd.to_numeric(mw[column], errors="coerce")

    duplicate_mask = mw.duplicated(
        subset=["state_fips", "quarterly_date"],
        keep=False,
    )
    if duplicate_mask.any():
        mw.loc[duplicate_mask].to_csv(
            TABLE_DIR / "minimum_wage_duplicates.csv",
            index=False,
        )
        raise ValueError(
            "Duplicate state-quarter minimum-wage observations found."
        )

    return mw


def retain_balanced_panel(panel: pd.DataFrame) -> pd.DataFrame:
    expected_periods = (
        (END_YEAR - START_YEAR + 1) * 4
    )

    coverage = (
        panel.groupby("state_fips")
        .agg(
            observed_quarters=("quarterly_date", "nunique"),
            first_year=("year", "min"),
            last_year=("year", "max"),
            state_name=("state_name", "first"),
            state_abbr=("state_abbr", "first"),
        )
        .reset_index()
    )

    coverage["balanced"] = (
        (coverage["observed_quarters"] == expected_periods)
        & (coverage["first_year"] == START_YEAR)
        & (coverage["last_year"] == END_YEAR)
    )

    coverage.to_csv(
        TABLE_DIR / "state_panel_coverage.csv",
        index=False,
    )

    valid_states = coverage.loc[
        coverage["balanced"], "state_fips"
    ]

    dropped_states = coverage.loc[
        ~coverage["balanced"],
        ["state_fips", "state_name", "observed_quarters"],
    ]

    if not dropped_states.empty:
        warnings.warn(
            "Dropping incomplete states from the estimation panel. "
            "See output/tables/state_panel_coverage.csv."
        )

    balanced = panel.loc[
        panel["state_fips"].isin(valid_states)
    ].copy()

    if balanced["state_fips"].nunique() < 45:
        raise ValueError(
            "Fewer than 45 complete state panels remain. "
            "Check the QCEW files and merge audit."
        )

    return balanced


def add_treatment_variables(panel: pd.DataFrame) -> pd.DataFrame:
    """
    Define an absorbing treatment based on the first substantial,
    state-initiated, binding minimum-wage increase.

    This treatment supports the estimand:
        Effect of entering the post-treatment period after the state's
        first qualifying increase during the analysis window.

    It does not estimate the marginal effect of every subsequent increase.
    """
    panel = panel.sort_values(
        ["state_fips", "year", "qtr"]
    ).copy()

    grouped = panel.groupby("state_fips", sort=False)

    panel["state_mw_change"] = grouped["state_mw_avg"].diff()
    panel["effective_mw_change"] = grouped["effective_mw"].diff()
    panel["log_effective_mw_change"] = grouped[
        "log_effective_mw"
    ].diff()

    panel["state_rate_binding"] = (
        panel["state_mw_avg"]
        > panel["federal_mw_avg"] + WAGE_TOLERANCE
    )

    panel["state_initiated_increase"] = (
        (panel["year"] >= TREATMENT_START_YEAR)
        & (panel["state_mw_change"] > WAGE_TOLERANCE)
        & (panel["effective_mw_change"] > WAGE_TOLERANCE)
        & panel["state_rate_binding"]
    )

    # Any economically meaningful increase: approximately >= 1%.
    panel["qualifying_hike_any"] = (
        panel["state_initiated_increase"]
        & (panel["log_effective_mw_change"] >= np.log(1.01))
    )

    # Main treatment: >= 5% increase in the effective quarterly average.
    panel["qualifying_hike_05"] = (
        panel["state_initiated_increase"]
        & (
            panel["log_effective_mw_change"]
            >= np.log(1.0 + MAJOR_HIKE_THRESHOLD)
        )
    )

    for suffix, indicator in [
        ("any", "qualifying_hike_any"),
        ("05", "qualifying_hike_05"),
    ]:
        first_treatment = (
            panel.loc[panel[indicator]]
            .groupby("state_fips")["t"]
            .min()
        )

        panel[f"G_hike_{suffix}"] = panel[
            "state_fips"
        ].map(first_treatment)

        panel[f"cohort_hike_{suffix}"] = (
            panel[f"G_hike_{suffix}"].fillna(0).astype(int)
        )

        panel[f"D_hike_{suffix}"] = (
            panel[f"G_hike_{suffix}"].notna()
            & (panel["t"] >= panel[f"G_hike_{suffix}"])
        ).astype(int)

        panel[f"rel_time_hike_{suffix}"] = np.where(
            panel[f"G_hike_{suffix}"].notna(),
            panel["t"] - panel[f"G_hike_{suffix}"],
            np.nan,
        )

    return panel


def build_panel() -> pd.DataFrame:
    setup_directories()

    qcew = load_qcew_state_private(
        QCEW_GLOB,
        START_YEAR,
        END_YEAR,
    )
    mw = load_minimum_wage_data(MW_PATH)

    merged = qcew.merge(
        mw,
        on=["state_fips", "quarterly_date"],
        how="left",
        validate="one_to_one",
        indicator=True,
    )

    merge_audit = (
        merged["_merge"]
        .value_counts(dropna=False)
        .rename_axis("merge_status")
        .reset_index(name="rows")
    )
    merge_audit.to_csv(
        TABLE_DIR / "merge_audit.csv",
        index=False,
    )

    unmatched = merged.loc[
        merged["_merge"] == "left_only",
        ["state_fips", "year", "qtr", "quarterly_date"],
    ]

    if not unmatched.empty:
        unmatched.to_csv(
            TABLE_DIR / "unmatched_qcew_state_quarters.csv",
            index=False,
        )
        raise ValueError(
            f"{len(unmatched)} QCEW state-quarters have no wage match. "
            "See output/tables/unmatched_qcew_state_quarters.csv."
        )

    panel = merged.drop(columns="_merge").copy()
    panel = panel.sort_values(
        ["state_fips", "year", "qtr"]
    ).reset_index(drop=True)

    # One-based period index. Zero remains reserved for never-treated cohorts.
    panel["t"] = (
        4 * (panel["year"] - START_YEAR)
        + panel["qtr"]
    ).astype(int)

    panel["effective_mw"] = panel[
        ["federal_mw_avg", "state_mw_avg"]
    ].max(axis=1)

    if panel["effective_mw"].isna().any():
        raise ValueError("Missing effective minimum-wage values found.")

    if (panel["effective_mw"] <= 0).any():
        raise ValueError("Non-positive effective minimum-wage values found.")

    panel["log_effective_mw"] = np.log(panel["effective_mw"])

    # Convert to 2012 dollars.
    panel["cpi_u"] = panel["year"].map(CPI_U)
    if panel["cpi_u"].isna().any():
        raise ValueError("Missing CPI values for one or more years.")

    panel["real_mw_2012"] = (
        panel["effective_mw"] * CPI_2012 / panel["cpi_u"]
    )
    panel["log_real_mw"] = np.log(panel["real_mw_2012"])

    panel["census_division"] = panel["state_fips"].map(
        CENSUS_DIVISION
    )
    if panel["census_division"].isna().any():
        raise ValueError("Census division mapping is incomplete.")

    panel = retain_balanced_panel(panel)
    panel = add_treatment_variables(panel)

    if panel.duplicated(["state_fips", "t"]).any():
        raise ValueError("Final panel has duplicate state-period rows.")

    expected_quarters = (END_YEAR - START_YEAR + 1) * 4
    state_counts = panel.groupby("state_fips")["t"].nunique()

    if not state_counts.eq(expected_quarters).all():
        raise ValueError("The final estimation panel is not balanced.")

    panel.to_csv(OUT_PATH, index=False)

    treatment_summary = (
        panel.groupby("state_fips", as_index=False)
        .agg(
            state_name=("state_name", "first"),
            state_abbr=("state_abbr", "first"),
            census_division=("census_division", "first"),
            G_hike_any=("G_hike_any", "first"),
            G_hike_05=("G_hike_05", "first"),
        )
    )
    treatment_summary["treated_any"] = (
        treatment_summary["G_hike_any"].notna()
    )
    treatment_summary["treated_05"] = (
        treatment_summary["G_hike_05"].notna()
    )

    treatment_summary.to_csv(
        TABLE_DIR / "state_treatment_summary.csv",
        index=False,
    )

    metadata = {
        "start_year": START_YEAR,
        "end_year": END_YEAR,
        "quarters_per_state": expected_quarters,
        "states_in_balanced_panel": int(
            panel["state_fips"].nunique()
        ),
        "observations": int(len(panel)),
        "treatment_start_year": TREATMENT_START_YEAR,
        "major_hike_threshold": MAJOR_HIKE_THRESHOLD,
        "states_treated_any": int(
            treatment_summary["treated_any"].sum()
        ),
        "states_treated_05": int(
            treatment_summary["treated_05"].sum()
        ),
        "states_never_treated_05": int(
            (~treatment_summary["treated_05"]).sum()
        ),
        "main_estimand": (
            "Effect of entering the post-treatment period after the first "
            "state-initiated, binding effective minimum-wage increase of "
            "at least 5 percent during 2011-2019."
        ),
    }

    with open(
        OUTPUT_DIR / "panel_metadata.json",
        "w",
        encoding="utf-8",
    ) as file:
        json.dump(metadata, file, indent=2)

    print("\nPanel successfully created")
    print(f"Path: {OUT_PATH}")
    print(f"Observations: {len(panel):,}")
    print(f"States/DC: {panel['state_fips'].nunique()}")
    print(
        "States with >=5% qualifying hike: "
        f"{treatment_summary['treated_05'].sum()}"
    )
    print(
        "Never-treated states under main definition: "
        f"{(~treatment_summary['treated_05']).sum()}"
    )

    return panel


if __name__ == "__main__":
    build_panel()
