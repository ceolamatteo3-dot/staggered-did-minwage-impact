args <- commandArgs(trailingOnly = TRUE)

panel_path <- if (length(args) >= 1) {
  args[1]
} else {
  "output/merged_qcew_mw_panel.csv"
}

required_packages <- c(
  "did",
  "fixest",
  "bacondecomp",
  "ggplot2",
  "dplyr",
  "readr",
  "tibble",
  "HonestDiD"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing R packages: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them with:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))\nOr from GitHub:\nremotes::install_github('asheshrambachan/HonestDiD')"
    )
  )
}

suppressPackageStartupMessages({
  library(did)
  library(fixest)
  library(bacondecomp)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(HonestDiD)
})

dir.create("output", showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

set.seed(20260807)

bootstrap_iterations <- as.integer(
  Sys.getenv("DID_BITERS", unset = "1999")
)

if (is.na(bootstrap_iterations) || bootstrap_iterations < 999) {
  bootstrap_iterations <- 1999
}

message("Reading panel: ", panel_path)

panel <- read_csv(panel_path, show_col_types = FALSE) |>
  mutate(
    state_fips = as.integer(state_fips),
    year = as.integer(year),
    qtr = as.integer(qtr),
    t = as.integer(t),
    log_emp = as.numeric(log_emp),
    avg_emplvl = as.numeric(avg_emplvl),
    effective_mw = as.numeric(effective_mw),
    log_real_mw = as.numeric(log_real_mw),
    cohort_hike_05 = as.integer(cohort_hike_05),
    D_hike_05 = as.integer(D_hike_05),
    rel_time_hike_05 = as.numeric(rel_time_hike_05),
    census_division = as.factor(census_division)
  ) |>
  arrange(state_fips, t)

required_columns <- c(
  "state_fips",
  "year",
  "qtr",
  "t",
  "log_emp",
  "avg_emplvl",
  "log_real_mw",
  "cohort_hike_05",
  "D_hike_05",
  "census_division"
)

missing_columns <- setdiff(required_columns, names(panel))

if (length(missing_columns) > 0) {
  stop(
    "Panel is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(panel[c("state_fips", "t")]) > 0) {
  stop("Duplicate state-period observations found.")
}

panel_counts <- panel |>
  count(state_fips, name = "periods")

if (n_distinct(panel_counts$periods) != 1) {
  stop("The panel is unbalanced across states.")
}

expected_rows <- n_distinct(panel$state_fips) * n_distinct(panel$t)

if (nrow(panel) != expected_rows) {
  stop("The panel does not contain every state-period combination.")
}

if (any(is.na(panel$log_emp))) {
  stop("Missing log employment values found.")
}

if (any(is.na(panel$cohort_hike_05))) {
  stop("Missing treatment-cohort values found.")
}

state_treatment <- panel |>
  group_by(state_fips) |>
  summarise(
    state_name = first(state_name),
    state_abbr = first(state_abbr),
    cohort_hike_05 = first(cohort_hike_05),
    treated = first(cohort_hike_05) > 0,
    .groups = "drop"
  )

cohort_sizes <- state_treatment |>
  filter(treated) |>
  count(cohort_hike_05, name = "states")

write_csv(
  cohort_sizes,
  "output/tables/cohort_sizes_r.csv"
)

if (nrow(cohort_sizes) > 0 && min(cohort_sizes$states) < 3) {
  warning(
    "At least one treatment cohort has fewer than three states. ",
    "Cohort-specific inference may be unstable."
  )
}

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

fixest_scc <- ssc(
  K.adj = TRUE,
  K.fixef = "nonnested",
  G.adj = TRUE,
  t.df = "min"
)

extract_fixest_coefficient <- function(fit, model_name) {
  ct <- coeftable(
    fit,
    vcov = ~state_fips,
    ssc = fixest_scc
  )
  
  if (!"log_real_mw" %in% rownames(ct)) {
    stop("log_real_mw was not estimated in model: ", model_name)
  }
  
  estimate <- unname(ct["log_real_mw", 1])
  std_error <- unname(ct["log_real_mw", 2])
  p_value <- unname(ct["log_real_mw", 4])
  n_obs <- nobs(fit)
  
  tibble(
    model = model_name,
    estimand = "Employment-minimum-wage elasticity",
    estimate = estimate,
    std_error = std_error,
    ci_low = estimate - 1.96 * std_error,
    ci_high = estimate + 1.96 * std_error,
    p_value = p_value,
    nobs = n_obs,
    clusters = n_distinct(panel$state_fips),
    inference = paste(
      "State-clustered SE with fixest finite-sample corrections;",
      "not CR2/CRV2"
    )
  )
}

extract_aggte_overall <- function(object, model_name) {
  estimate <- object$overall.att
  std_error <- object$overall.se
  
  tibble(
    model = model_name,
    estimate_log_points = estimate,
    std_error = std_error,
    ci_low = estimate - 1.96 * std_error,
    ci_high = estimate + 1.96 * std_error,
    approximate_percent_effect = 100 * (exp(estimate) - 1),
    approximate_percent_ci_low = 100 * (
      exp(estimate - 1.96 * std_error) - 1
    ),
    approximate_percent_ci_high = 100 * (
      exp(estimate + 1.96 * std_error) - 1
    )
  )
}

extract_dynamic_table <- function(object) {
  critical_value <- object$crit.val.egt
  
  if (is.null(critical_value) || length(critical_value) == 0) {
    critical_value <- 1.96
  }
  
  tibble(
    event_time = object$egt,
    estimate = object$att.egt,
    std_error = object$se.egt,
    simultaneous_ci_low = (
      object$att.egt - critical_value * object$se.egt
    ),
    simultaneous_ci_high = (
      object$att.egt + critical_value * object$se.egt
    ),
    critical_value = critical_value
  )
}

# ------------------------------------------------------------------
# HonestDiD wrapper for Callaway & Sant'Anna (did::aggte)
# ------------------------------------------------------------------

honest_did <- function(...) UseMethod("honest_did")

honest_did.AGGTEobj <- function(
    es,
    e = 0,
    type = c("smoothness", "relative_magnitude"),
    gridPoints = 100,
    ...
) {
  type <- match.arg(type)
  
  if (es$type != "dynamic") {
    stop("honest_did requires a dynamic event study object from aggte()")
  }
  
  if (es$DIDparams$base_period != "universal") {
    stop("Use a universal base period for honest_did (base_period = 'universal')")
  }
  
  es_inf_func <- es$inf.function$dynamic.inf.func.e
  n <- nrow(es_inf_func)
  V <- t(es_inf_func) %*% es_inf_func / n / n
  
  referencePeriod <- -1
  hasReference <- any(es$egt == referencePeriod)
  
  if (hasReference) {
    referencePeriodIndex <- which(es$egt == referencePeriod)
    V <- V[-referencePeriodIndex, -referencePeriodIndex]
    beta <- es$att.egt[-referencePeriodIndex]
  } else {
    beta <- es$att.egt
  }
  
  nperiods <- nrow(V)
  npre <- sum(1 * (es$egt < referencePeriod))
  npost <- nperiods - npre
  
  if (e >= 0) {
    e_index <- e + 1
  } else {
    stop("e must be non-negative (post-treatment event period)")
  }
  
  baseVec1 <- HonestDiD::basisVector(index = e_index, size = npost)
  
  orig_ci <- HonestDiD::constructOriginalCS(
    betahat = beta,
    sigma = V,
    numPrePeriods = npre,
    numPostPeriods = npost,
    l_vec = baseVec1
  )
  
  if (type == "relative_magnitude") {
    robust_ci <- HonestDiD::createSensitivityResults_relativeMagnitudes(
      betahat = beta,
      sigma = V,
      numPrePeriods = npre,
      numPostPeriods = npost,
      l_vec = baseVec1,
      gridPoints = gridPoints,
      ...
    )
  } else if (type == "smoothness") {
    robust_ci <- HonestDiD::createSensitivityResults(
      betahat = beta,
      sigma = V,
      numPrePeriods = npre,
      numPostPeriods = npost,
      l_vec = baseVec1,
      ...
    )
  }
  
  return(list(robust_ci = robust_ci, orig_ci = orig_ci, type = type))
}

# ------------------------------------------------------------------
# 1. Continuous-treatment TWFE benchmarks
# ------------------------------------------------------------------

message("=== 1. CONTINUOUS-TREATMENT BENCHMARKS ===")

continuous_twfe <- feols(
  log_emp ~ log_real_mw | state_fips + t,
  data = panel
)

continuous_state_seasonality <- feols(
  log_emp ~ log_real_mw | state_fips^qtr + t,
  data = panel
)

continuous_division_time <- feols(
  log_emp ~ log_real_mw | state_fips + census_division^t,
  data = panel
)

continuous_state_trends <- feols(
  log_emp ~ log_real_mw | state_fips[t] + t,
  data = panel
)

continuous_results <- bind_rows(
  extract_fixest_coefficient(
    continuous_twfe,
    "State FE + national quarter FE"
  ),
  extract_fixest_coefficient(
    continuous_state_seasonality,
    "State-by-calendar-quarter seasonality + national quarter FE"
  ),
  extract_fixest_coefficient(
    continuous_division_time,
    "State FE + Census-division-by-quarter FE"
  ),
  extract_fixest_coefficient(
    continuous_state_trends,
    "State FE + national quarter FE + state-specific linear trends"
  )
)

write_csv(
  continuous_results,
  "output/tables/continuous_twfe_fixest.csv"
)

capture.output(
  etable(
    continuous_twfe,
    continuous_state_seasonality,
    continuous_division_time,
    continuous_state_trends,
    vcov = ~state_fips,
    ssc = fixest_scc,
    headers = c(
      "TWFE",
      "State seasonality",
      "Division-time FE",
      "State trends"
    )
  ),
  file = "output/tables/continuous_twfe_fixest.txt"
)

print(continuous_results)

# ------------------------------------------------------------------
# 2. Goodman-Bacon decomposition
# ------------------------------------------------------------------

message("=== 2. GOODMAN-BACON DECOMPOSITION ===")

if (n_distinct(panel$D_hike_05) > 1) {
  bacon_result <- bacon(
    log_emp ~ D_hike_05,
    data = panel,
    id_var = "state_fips",
    time_var = "t"
  )
  
  bacon_table <- as_tibble(bacon_result)
  
  write_csv(
    bacon_table,
    "output/tables/goodman_bacon_main_treatment.csv"
  )
  
  bacon_plot <- ggplot(
    bacon_table,
    aes(
      x = weight,
      y = estimate,
      color = type
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_point(
      size = 2.6,
      alpha = 0.8
    ) +
    labs(
      title = "Goodman-Bacon decomposition",
      subtitle = paste(
        "Legacy TWFE for first state-initiated binding hike",
        "of at least 5%"
      ),
      x = "TWFE weight",
      y = "2x2 DiD estimate",
      color = "Comparison type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom"
    )
  
  ggsave(
    "output/figures/goodman_bacon_main_treatment.png",
    bacon_plot,
    width = 8.5,
    height = 5.5,
    dpi = 300
  )
} else {
  warning(
    "D_hike_05 has no treatment variation. ",
    "Goodman-Bacon decomposition was skipped."
  )
}

# ------------------------------------------------------------------
# 3. Callaway-Sant'Anna group-time ATT
# ------------------------------------------------------------------

message("=== 3. CALLAWAY-SANT'ANNA GROUP-TIME ATT ===")

treated_states <- sum(state_treatment$treated)
never_treated_states <- sum(!state_treatment$treated)

if (treated_states == 0) {
  stop("No treated states exist under the main treatment definition.")
}

if (never_treated_states == 0) {
  warning(
    "No never-treated states exist. Identification will use only ",
    "not-yet-treated states."
  )
}

cs <- att_gt(
  yname = "log_emp",
  tname = "t",
  idname = "state_fips",
  gname = "cohort_hike_05",
  xformla = ~1,
  data = panel,
  panel = TRUE,
  allow_unbalanced_panel = FALSE,
  control_group = "notyettreated",
  anticipation = 0,
  base_period = "universal",
  est_method = "reg",
  bstrap = TRUE,
  cband = TRUE,
  biters = bootstrap_iterations,
  clustervars = "state_fips",
  print_details = TRUE
)

capture.output(
  summary(cs),
  file = "output/tables/cs_group_time_summary.txt"
)

group_time_table <- tibble(
  cohort = cs$group,
  time = cs$t,
  event_time = cs$t - cs$group,
  estimate = cs$att,
  std_error = cs$se
)

write_csv(
  group_time_table,
  "output/tables/cs_group_time_att.csv"
)

# ------------------------------------------------------------------
# 4. Native Callaway-Sant'Anna pretrend diagnostic
# ------------------------------------------------------------------

message("=== 4. PRE-TREATMENT DIAGNOSTIC ===")

pretrend_p_value <- if (!is.null(cs$Wpval)) {
  cs$Wpval
} else {
  NA_real_
}

pretrend_table <- tibble(
  test = "Native did-package Wald pretest",
  p_value = pretrend_p_value,
  interpretation = paste(
    "Failure to reject is supporting evidence only;",
    "it does not prove post-treatment parallel trends."
  )
)

write_csv(
  pretrend_table,
  "output/tables/cs_pretrend_test.csv"
)

print(pretrend_table)

# ------------------------------------------------------------------
# 5. Modern dynamic event study
# ------------------------------------------------------------------

message("=== 5. MODERN DYNAMIC EVENT STUDY ===")

cs_dynamic <- aggte(
  cs,
  type = "dynamic",
  min_e = -8,
  max_e = 8,
  na.rm = TRUE
)

capture.output(
  summary(cs_dynamic),
  file = "output/tables/cs_dynamic_summary.txt"
)

dynamic_table <- extract_dynamic_table(cs_dynamic)

write_csv(
  dynamic_table,
  "output/tables/cs_dynamic_event_study.csv"
)

dynamic_plot <- ggdid(cs_dynamic) +
  labs(
    title = "Callaway-Sant'Anna dynamic effects",
    subtitle = paste(
      "First state-initiated binding effective MW increase >= 5%;",
      "simultaneous 95% confidence bands"
    ),
    x = "Quarters relative to first qualifying increase",
    y = "ATT on log food-services employment"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  "output/figures/cs_dynamic_event_study.png",
  dynamic_plot,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

# ------------------------------------------------------------------
# 6. First-year ATT with balanced event-time composition
# ------------------------------------------------------------------

message("=== 6. COHORT-BALANCED FIRST-YEAR ATT ===")

cs_first_year <- aggte(
  cs,
  type = "dynamic",
  min_e = 0,
  max_e = 3,
  balance_e = 3,
  na.rm = TRUE
)

capture.output(
  summary(cs_first_year),
  file = "output/tables/cs_first_year_summary.txt"
)

first_year_result <- extract_aggte_overall(
  cs_first_year,
  paste(
    "C&S cohort-balanced first-year ATT:",
    "event times 0-3"
  )
)

write_csv(
  first_year_result,
  "output/tables/cs_first_year_att.csv"
)

first_year_dynamic <- extract_dynamic_table(cs_first_year)

write_csv(
  first_year_dynamic,
  "output/tables/cs_first_year_by_event_time.csv"
)

print(first_year_result)

# ------------------------------------------------------------------
# 7. Overall group-weighted ATT
# ------------------------------------------------------------------

message("=== 7. OVERALL GROUP-WEIGHTED ATT ===")

cs_group <- aggte(
  cs,
  type = "group",
  na.rm = TRUE
)

capture.output(
  summary(cs_group),
  file = "output/tables/cs_group_aggregation_summary.txt"
)

group_result <- extract_aggte_overall(
  cs_group,
  "C&S group-weighted overall ATT"
)

write_csv(
  group_result,
  "output/tables/cs_group_weighted_overall_att.csv"
)

group_effects <- tibble(
  cohort = cs_group$egt,
  estimate = cs_group$att.egt,
  std_error = cs_group$se.egt
)

write_csv(
  group_effects,
  "output/tables/cs_effects_by_cohort.csv"
)

print(group_result)

# ------------------------------------------------------------------
# 8. Compact headline results
# ------------------------------------------------------------------

headline_results <- bind_rows(
  first_year_result |>
    transmute(
      model,
      estimate_log_points,
      std_error,
      ci_low,
      ci_high,
      approximate_percent_effect,
      approximate_percent_ci_low,
      approximate_percent_ci_high
    ),
  group_result |>
    transmute(
      model,
      estimate_log_points,
      std_error,
      ci_low,
      ci_high,
      approximate_percent_effect,
      approximate_percent_ci_low,
      approximate_percent_ci_high
    )
)

write_csv(
  headline_results,
  "output/tables/headline_causal_results.csv"
)

# ------------------------------------------------------------------
# 9. Rambachan & Roth (2023) HonestDiD sensitivity analysis
# ------------------------------------------------------------------

message("=== 9. RAMBACHAN & ROTH (2023) HONESTDID SENSITIVITY ANALYSIS ===")

# 9.1 Relative Magnitude Restrictions (Mbar)
message("Evaluating Relative Magnitude Restrictions (Mbar)...")

hd_rm <- honest_did(
  es = cs_dynamic,
  e = 0,
  type = "relative_magnitude",
  Mbarvec = seq(from = 0, to = 2, by = 0.5)
)

rm_results_table <- as_tibble(hd_rm$robust_ci)

write_csv(
  rm_results_table,
  "output/tables/honestdid_relative_magnitudes.csv"
)

rm_plot <- HonestDiD::createSensitivityPlot_relativeMagnitudes(
  hd_rm$robust_ci,
  hd_rm$orig_ci
) +
  labs(
    title = "HonestDiD: Relative Magnitude Bounds (Rambachan & Roth, 2023)",
    subtitle = expression("Sensitivity of ATT at e=0 to relative violations of parallel trends (" * bar(M) * ")"),
    x = expression("Relative magnitude bound (" * bar(M) * ")"),
    y = "95% Robust CI for ATT (e=0)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  "output/figures/honestdid_relative_magnitudes.png",
  rm_plot,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

# 9.2 Smoothness / Slope Restriction (M)
message("Evaluating Smoothness / Slope Restrictions (M)...")

hd_sd <- honest_did(
  es = cs_dynamic,
  e = 0,
  type = "smoothness",
  Mvec = seq(from = 0, to = 0.05, by = 0.01)
)

sd_results_table <- as_tibble(hd_sd$robust_ci)

write_csv(
  sd_results_table,
  "output/tables/honestdid_smoothness.csv"
)

sd_plot <- HonestDiD::createSensitivityPlot(
  hd_sd$robust_ci,
  hd_sd$orig_ci
) +
  labs(
    title = "HonestDiD: Smoothness Restrictions (Rambachan & Roth, 2023)",
    subtitle = "Sensitivity of ATT at e=0 to maximum period-to-period trend slope changes",
    x = "Maximum slope change parameter (M)",
    y = "95% Robust CI for ATT (e=0)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  "output/figures/honestdid_smoothness.png",
  sd_plot,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

print(rm_results_table)

message("Analysis complete.")
message(
  "Main estimand: post-treatment effect of the first state-initiated, ",
  "binding effective minimum-wage increase of at least 5%."
)
message(
  "Rambachan & Roth (2023) HonestDiD sensitivity analysis successfully completed."
)