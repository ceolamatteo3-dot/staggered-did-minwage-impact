from __future__ import annotations

import json
import shutil
import subprocess
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from linearmodels.panel import PanelOLS


PROJECT_DIR =  Path('C:/Users/Federico/Desktop/causal_inference_wage_employment')

PANEL_PATH = PROJECT_DIR / "output" / "merged_qcew_mw_panel.csv"
OUTPUT_DIR = PROJECT_DIR / "output"
TABLE_DIR = OUTPUT_DIR / "tables"
FIGURE_DIR = OUTPUT_DIR / "figures"

OUTCOME = "log_emp"
CLUSTER = "state_fips"

REQUIRED_COLUMNS = [
    "state_fips",
    "state_name",
    "year",
    "qtr",
    "t",
    "avg_emplvl",
    "log_emp",
    "effective_mw",
    "real_mw_2012",
    "log_real_mw",
    "cohort_hike_05",
    "D_hike_05",
    "rel_time_hike_05",
    "census_division",
]


def setup_directories() -> None:
    for path in [OUTPUT_DIR, TABLE_DIR, FIGURE_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def load_panel() -> pd.DataFrame:
    if not PANEL_PATH.exists():
        raise FileNotFoundError(
            f"Panel not found: {PANEL_PATH}\n"
            "Run build_panel.py first."
        )

    df = pd.read_csv(PANEL_PATH)

    missing = sorted(set(REQUIRED_COLUMNS) - set(df.columns))
    if missing:
        raise ValueError(f"Panel is missing columns: {missing}")

    numeric_columns = [
        "state_fips",
        "year",
        "qtr",
        "t",
        "avg_emplvl",
        "log_emp",
        "effective_mw",
        "real_mw_2012",
        "log_real_mw",
        "cohort_hike_05",
        "D_hike_05",
        "rel_time_hike_05",
    ]

    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.dropna(
        subset=[
            "state_fips",
            "t",
            "log_emp",
            "log_real_mw",
            "cohort_hike_05",
        ]
    ).copy()

    df["state_fips"] = df["state_fips"].astype(int)
    df["t"] = df["t"].astype(int)
    df["cohort_hike_05"] = df["cohort_hike_05"].astype(int)
    df["D_hike_05"] = df["D_hike_05"].astype(int)

    if df.duplicated(["state_fips", "t"]).any():
        raise ValueError("Duplicate state-period observations found.")

    periods_per_state = df.groupby("state_fips")["t"].nunique()
    if periods_per_state.nunique() != 1:
        raise ValueError(
            "The estimation panel is unbalanced. Re-run build_panel.py."
        )

    expected_n = (
        df["state_fips"].nunique()
        * df["t"].nunique()
    )
    if len(df) != expected_n:
        raise ValueError(
            "The panel does not contain every state-period combination."
        )

    return df.sort_values(["state_fips", "t"]).reset_index(drop=True)


def run_continuous_twfe(df: pd.DataFrame) -> None:
    """
    Secondary continuous-treatment benchmark.

    The coefficient is an employment-minimum-wage elasticity under the
    strong assumption that, conditional on state and time fixed effects,
    changes in log minimum wage are exogenous to state-specific employment
    shocks.

    This is not the same estimand as the binary first-hike DiD.
    """
    model_data = df[
        ["state_fips", "t", OUTCOME, "log_real_mw"]
    ].dropna().copy()

    panel = model_data.set_index(
        ["state_fips", "t"]
    ).sort_index()

    model = PanelOLS(
        dependent=panel[OUTCOME],
        exog=panel[["log_real_mw"]],
        entity_effects=True,
        time_effects=True,
        drop_absorbed=True,
    ).fit(
        cov_type="clustered",
        cluster_entity=True,
        debiased=True,
    )

    estimate = float(model.params["log_real_mw"])
    std_error = float(model.std_errors["log_real_mw"])

    result = pd.DataFrame(
        [
            {
                "model": "Continuous TWFE benchmark",
                "estimand": (
                    "Elasticity of food-services employment with respect "
                    "to the effective minimum wage"
                ),
                "estimate": estimate,
                "std_error": std_error,
                "ci_low": estimate - 1.96 * std_error,
                "ci_high": estimate + 1.96 * std_error,
                "nobs": int(model.nobs),
                "clusters": int(df["state_fips"].nunique()),
                "note": (
                    "Secondary associational/causal benchmark under "
                    "strong continuous-treatment exogeneity assumptions."
                ),
            }
        ]
    )

    result.to_csv(
        TABLE_DIR / "continuous_twfe_python.csv",
        index=False,
    )

    with open(
        TABLE_DIR / "continuous_twfe_python_summary.txt",
        "w",
        encoding="utf-8",
    ) as file:
        file.write(str(model.summary))

    print("\nContinuous TWFE benchmark")
    print(result.to_string(index=False))


def export_design_diagnostics(df: pd.DataFrame) -> None:
    state_level = (
        df.groupby("state_fips", as_index=False)
        .agg(
            state_name=("state_name", "first"),
            state_abbr=("state_abbr", "first"),
            census_division=("census_division", "first"),
            cohort_hike_05=("cohort_hike_05", "first"),
            mean_employment=("avg_emplvl", "mean"),
            first_effective_mw=("effective_mw", "first"),
            last_effective_mw=("effective_mw", "last"),
        )
    )

    state_level["treated"] = (
        state_level["cohort_hike_05"] > 0
    )

    state_level.to_csv(
        TABLE_DIR / "analysis_state_diagnostics.csv",
        index=False,
    )

    cohort_sizes = (
        state_level.loc[state_level["treated"]]
        .groupby("cohort_hike_05")
        .size()
        .rename("states")
        .reset_index()
    )

    if not cohort_sizes.empty:
        cohort_sizes["small_cohort_warning"] = (
            cohort_sizes["states"] < 3
        )

    cohort_sizes.to_csv(
        TABLE_DIR / "cohort_sizes_main_treatment.csv",
        index=False,
    )

    diagnostics = {
        "observations": int(len(df)),
        "states": int(df["state_fips"].nunique()),
        "periods": int(df["t"].nunique()),
        "treated_states": int(state_level["treated"].sum()),
        "never_treated_states": int(
            (~state_level["treated"]).sum()
        ),
        "treatment_cohorts": int(
            state_level.loc[
                state_level["treated"], "cohort_hike_05"
            ].nunique()
        ),
        "smallest_treatment_cohort": (
            int(cohort_sizes["states"].min())
            if not cohort_sizes.empty
            else None
        ),
    }

    with open(
        OUTPUT_DIR / "analysis_diagnostics.json",
        "w",
        encoding="utf-8",
    ) as file:
        json.dump(diagnostics, file, indent=2)

    if diagnostics["never_treated_states"] == 0:
        warnings.warn(
            "No never-treated states exist under the main definition. "
            "Callaway-Sant'Anna must rely entirely on not-yet-treated units."
        )

    if (
        diagnostics["smallest_treatment_cohort"] is not None
        and diagnostics["smallest_treatment_cohort"] < 3
    ):
        warnings.warn(
            "At least one treatment cohort contains fewer than three states. "
            "Interpret cohort-specific estimates cautiously."
        )


def find_rscript() -> str | None:
    rscript = shutil.which("Rscript")
    if rscript:
        return rscript

    windows_r_directory = Path(r"C:\Program Files\R")
    if windows_r_directory.exists():
        candidates = sorted(
            windows_r_directory.glob("R-*/bin/Rscript.exe"),
            reverse=True,
        )
        if candidates:
            return str(candidates[0])

    return None


def run_r_modern_did() -> None:
    script = PROJECT_DIR / "modern_did.R"

    if not script.exists():
        warnings.warn(f"R script not found: {script}")
        return

    rscript = find_rscript()
    if not rscript:
        warnings.warn(
            "Rscript was not found. Install R and required R packages, "
            "or run modern_did.R manually."
        )
        return

    command = [
        rscript,
        str(script),
        str(PANEL_PATH),
    ]

    print(f"\nRunning modern DiD analysis with: {rscript}")

    try:
        subprocess.run(
            command,
            cwd=PROJECT_DIR,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"modern_did.R failed with exit code {exc.returncode}."
        ) from exc

    print("R modern-DiD analysis completed successfully.")


def main() -> None:
    setup_directories()
    df = load_panel()

    export_design_diagnostics(df)
    run_continuous_twfe(df)
    run_r_modern_did()


if __name__ == "__main__":
    main()
