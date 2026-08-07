# Causal Effect of State Minimum Wage Increases on Food-Services Employment (2009–2019)

## Executive Summary

This repository contains the replication package evaluating the causal impact of state-initiated minimum wage increases **≥ 5%** on employment in the **Food Services and Drinking Places industry (NAICS 722)**[cite: 1, 2]. 

Using a balanced quarterly panel of 51 U.S. state entities over 44 quarters (2009 Q1 – 2019 Q4, **N = 2,244**), we employ the heterogeneity-robust **Callaway & Sant'Anna (2021)** Difference-in-Differences (DiD) estimator[cite: 1, 3], supplemented by **Goodman-Bacon (2021)** decompositions[cite: 1, 3] and **Rambachan & Roth (2023)** `HonestDiD` sensitivity analysis[cite: 1, 3].

### Key Causal Findings
* **First-Year Balanced ATT ($e \in [0,3]$):** -0.08 log points (SE = 0.34, 95% CI: [-0.0076, +0.0059]). Statistically indistinguishable from zero (p = 0.81).
* **Overall Group-Weighted ATT:** -0.42 log points (SE = 0.67, 95% CI: [-0.0174, +0.0090]). Statistically indistinguishable from zero (p = 0.54).
* **On-Impact Effect ($e = 0$):** A small, transitory drop of -1.21% (SE = 0.42%, p < 0.01), which completely dissipates by quarter $e = 1$ (+0.21%).
* **Conclusion:** Moderate, binding minimum wage increases (≥ 5%) do not lead to structural job losses in the food-services sector[cite: 1, 2].

---

## Headline Results Table

| Estimand / Model | Point Estimate | Std. Error | 95% Confidence Interval | Stat. Sig. ($p < 0.05$) |
| :--- | :---: | :---: | :---: | :---: |
| **C&S First-Year Balanced ATT ($e \in [0,3]$)** | **-0.00083** | **0.00344** | **[-0.0076, +0.0059]** | **No** |
| **C&S Group-Weighted Overall ATT** | **-0.00417** | **0.00674** | **[-0.0174, +0.0090]** | **No** |
| **C&S On-Impact Quarter Effect ($e = 0$)** | **-0.01215** | **0.00415** | **[-0.0216, -0.0027]** | **Yes** |
| Continuous TWFE (State + Quarter FE) | +0.04307 | 0.06178 | [-0.0780, +0.1642] | No |
| Continuous TWFE (State Seasonality FE) | +0.04416 | 0.06222 | [-0.0777, +0.1660] | No |
| Continuous TWFE (Census Division-Time FE) | +0.07147 | 0.05415 | [-0.0347, +0.1776] | No |
| Continuous TWFE (+ State Trends) | -0.06706 | 0.03831 | [-0.1421, +0.0080] | No |

---

## Methodological Justifications

### 1. Relegating Two-Way Fixed Effects (TWFE) to Benchmark Status
Under staggered adoption and dynamic treatment effects, linear TWFE regressions suffer from the **Goodman-Bacon (2021)** negative weighting problem[cite: 1, 3]. TWFE compares newly treated units against previously treated units, introducing severe bias and potential sign-reversals (Sun & Abraham 2021; de Chaisemartin & D'Haultfœuille 2020)[cite: 1].

### 2. Choice of Control Group (`control_group = "notyettreated"`)
Over the 2009–2019 panel, pure "never-treated" states under $D_{\text{hike}, 05}$ are limited and non-randomly geographically clustered[cite: 1, 3]. Incorporating not-yet-treated states expands the control pool at each time $t$, significantly improving statistical power while preserving non-parametric identification[cite: 1, 3].

### 3. Selection of Base Period (`base_period = "universal"`)
Setting `base_period = "universal"` in `att_gt()` anchors pre-treatment comparisons to $t = g - 1$[cite: 1, 3]. This is an essential mathematical prerequisite for calculating variance-covariance matrices compatible with **Rambachan & Roth (2023)** sensitivity bounds[cite: 1, 3].

### 4. Cohort Composition Balancing (`balance_e = 3`)
In dynamic event studies, younger cohorts exit the sample at higher relative event times ($e$)[cite: 1, 3]. Unbalanced event-study summaries suffer from **composition bias**[cite: 1]. Restricting aggregation to $e \in \{0,1,2,3\}$ across a balanced set of cohorts isolates a clean 1-year causal horizon[cite: 1, 3].

### 5. Rejection of State Linear Trends as Primary Specification
As proved by **Meer & West (2016)**, if minimum wage increases affect the *growth rate* rather than just the *level* of employment, state linear trends absorb the post-treatment outcome trajectory itself (trend-soaking bias), driving true causal effects toward zero or artificially reversing their sign[cite: 1].

### 6. Rejection of Division-Time & State-Quarter Fixed Effects
Saturating linear TWFE models with division-quarter FEs (396 parameters) or state-quarter FEs (204 parameters) consumes excessive degrees of freedom and forces estimation off high-frequency real-wage inflation noise without resolving underlying staggered TWFE negative weighting[cite: 1].

### 7. Implementation of HonestDiD Sensitivity Analysis
Pre-trend Wald tests have low power against non-linear trend deviations and introduce pre-test selection bias (Roth 2022)[cite: 1]. **Rambachan & Roth (2023)** replaces binary parallel trend tests with robust confidence bounds under bounded violations[cite: 1, 3]:
* **Relative Magnitude Breakdown Value ($\bar{M}^* \approx 0.45$):** The on-impact effect loses statistical significance if post-treatment trend violations exceed 45% of the maximum pre-treatment fluctuation[cite: 1, 3].
* **Smoothness Breakdown Value ($M^* \approx 0.005$):** The effect loses significance if trend slope changes exceed 0.5% per quarter[cite: 1, 3].

---

## Repository Structure & Execution

```text
causal_inference_wage_employment/
├── data/
│   ├── qcew/                      # BLS QCEW industry CSV files (NAICS 722)
│   └── mw/
│       └── mw_state_quarterly.xlsx # UKCPR quarterly state minimum wage database
├── output/
│   ├── tables/                    # CSV & TXT summaries (C&S ATTs, TWFE, diagnostics)
│   └── figures/                   # Bacon decomposition & HonestDiD sensitivity plots
├── build_panel.py                 # Panel assembly, CPI deflation, treatment cohort assignment
├── causal_analysis.py             # Pipeline orchestrator & continuous TWFE benchmark runner
└── modern_did.R                   # C&S DiD, Goodman-Bacon, & HonestDiD sensitivity analysis
