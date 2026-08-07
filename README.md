# Causal Effect of State Minimum Wage Increases on Food-Services Employment (2009–2019)

## Executive Summary

This repository contains the replication package evaluating the causal impact of state-initiated minimum wage increases ($\ge 5 % \%$) on employment in the **Food Services and Drinking Places industry (NAICS 722)**. 

Using a balanced quarterly panel of 51 U.S. state entities over 44 quarters (2009 Q1 – 2019 Q4, $N = 2,244$), we employ the heterogeneity-robust **Callaway & Sant'Anna (2021)** Difference-in-Differences (DiD) estimator, supplemented by **Goodman-Bacon (2021)** decompositions and **Rambachan & Roth (2023)** `HonestDiD` sensitivity analysis.

### Key Causal Findings
* **First-Year Balanced ATT ($e \in [0,3]$):** $-0.08\%$ log points ($SE = 0.34\%$, $95\%\text{ CI}: [-0.75\%, +0.59\%]$). Statistically indistinguishable from zero ($p = 0.81$).
* **Overall Group-Weighted ATT:** $-0.42\%$ log points ($SE = 0.67\%$, $95\%\text{ CI}: [-1.72\%, +0.91\%]$). Statistically indistinguishable from zero ($p = 0.54$).
* **On-Impact Effect ($e = 0$):** A small, transitory drop of $-1.21\%$ ($SE = 0.42\%$, $p < 0.01$), which completely dissipates by quarter $e = 1$ ($+0.21\%$).
* **Conclusion:** Moderate, binding minimum wage increases ($\ge 5\%$) do not lead to structural job losses in the food-services sector.

---

## Headline Results Table

| Estimand / Model | Point Estimate | Std. Error | 95% Confidence Interval | Stat. Sig. ($p < 0.05$) |
| :--- | :---: | :---: | :---: | :---: |
| **C&S First-Year Balanced ATT ($e \in [0,3]$)** | **$-0.00083$** | **$0.00344$** | **$[-0.0076, +0.0059]$** | **No** |
| **C&S Group-Weighted Overall ATT** | **$-0.00417$** | **$0.00674$** | **$[-0.0174, +0.0090]$** | **No** |
| **C&S On-Impact Quarter Effect ($e = 0$)** | **$-0.01215$** | **$0.00415$** | **$[-0.0216, -0.0027]$** | **Yes** |
| Continuous TWFE (State + Quarter FE) | $+0.04307$ | $0.06178$ | $[-0.0780, +0.1642]$ | No |
| Continuous TWFE (State Seasonality FE) | $+0.04416$ | $0.06222$ | $[-0.0777, +0.1660]$ | No |
| Continuous TWFE (Census Division-Time FE) | $+0.07147$ | $0.05415$ | $[-0.0347, +0.1776]$ | No |
| Continuous TWFE (+ State Trends) | $-0.06706$ | $0.03831$ | $[-0.1421, +0.0080]$ | No |

---

## Methodological Justifications

### 1. Relegating Two-Way Fixed Effects (TWFE) to Benchmark Status
Under staggered adoption and dynamic treatment effects, linear TWFE regressions suffer from the **Goodman-Bacon (2021)** negative weighting problem. TWFE compares newly treated units against previously treated units, introducing severe bias and potential sign-reversals (Sun & Abraham 2021; de Chaisemartin & D'Haultfœuille 2020).

### 2. Choice of Control Group (`control_group = "notyettreated"`)
Over the 2009–2019 panel, pure "never-treated" states under $D_{\text{hike}, 05}$ are limited and non-randomly geographically clustered. Incorporating not-yet-treated states expands the control pool at each time $t$, significantly improving statistical power while preserving non-parametric identification.

### 3. Selection of Base Period (`base_period = "universal"`)
Setting `base_period = "universal"` in `att_gt()` anchors pre-treatment comparisons to $t = g - 1$. This is an essential mathematical prerequisite for calculating variance-covariance matrices compatible with **Rambachan & Roth (2023)** sensitivity bounds.

### 4. Cohort Composition Balancing (`balance_e = 3`)
In dynamic event studies, younger cohorts exit the sample at higher relative event times ($e$). Unbalanced event-study summaries suffer from **composition bias**. Restricting aggregation to $e \in \{0,1,2,3\}$ across a balanced set of cohorts isolates a clean 1-year causal horizon.

### 5. Rejection of State Linear Trends as Primary Specification
As proved by **Meer & West (2016)**, if minimum wage increases affect the *growth rate* rather than just the *level* of employment, state linear trends absorb the post-treatment outcome trajectory itself (trend-soaking bias), driving true causal effects toward zero or artificially reversing their sign.

### 6. Rejection of Division-Time & State-Quarter Fixed Effects
Saturating linear TWFE models with division-quarter FEs ($396$ parameters) or state-quarter FEs ($204$ parameters) consumes excessive degrees of freedom and forces estimation off high-frequency real-wage inflation noise without resolving underlying staggered TWFE negative weighting.

### 7. Implementation of HonestDiD Sensitivity Analysis
Pre-trend Wald tests have low power against non-linear trend deviations and introduce pre-test selection bias (Roth 2022). **Rambachan & Roth (2023)** replaces binary parallel trend tests with robust confidence bounds under bounded violations:
* **Relative Magnitude Breakdown Value ($\bar{M}^* \approx 0.45$):** The on-impact effect loses statistical significance if post-treatment trend violations exceed $45\%$ of the maximum pre-treatment fluctuation.
* **Smoothness Breakdown Value ($M^* \approx 0.005$):** The effect loses significance if trend slope changes exceed $0.5\%$ per quarter.

---

## Repository Structure & Execution
