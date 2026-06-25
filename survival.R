# ============================================================================
# Survival Analysis of Event-Level Expression
# ============================================================================
# Purpose: Identify survival-associated events using cutpoint-based
#          dichotomization, Cox regression, and log-rank testing.
# ============================================================================

rm(list = ls())

library(tidyverse)
library(survival)
library(survminer)

# ---------------------------------------------------------------------------
# Input files and output settings
# ---------------------------------------------------------------------------
expression_file <- 'subs_normexpr_filt_blacklist.csv'
survival_file <- 'survival.csv'

# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------
expression_data <- read.csv(expression_file, row.names = NULL, sep = ',', header = TRUE)
survival_data <- read.csv(survival_file, row.names = NULL, sep = ',', header = TRUE)

prepare_survival_matrix <- function(expression_df, survival_df) {
  expression_df %>%
    select(error, contains('Tumor')) %>%
    mutate(across(where(is.numeric), ~ na_if(.x, 0))) %>%
    pivot_longer(cols = contains('Tumor'), names_to = 'Sample.ID', values_to = 'Value') %>%
    pivot_wider(names_from = error, values_from = Value) %>%
    mutate(
      Sample.ID = str_remove(Sample.ID, '\\.Primary\\.Tumor'),
      Sample.ID = str_replace_all(Sample.ID, '\\.', '-'),
      Sample.ID = str_replace_all(Sample.ID, ' ', '')
    ) %>%
    rename_with(~ str_replace_all(., '-', '.')) %>%
    select(where(~ mean(is.na(.)) <= 0.8)) %>%
    left_join(survival_df, by = 'Sample.ID') %>%
    column_to_rownames(var = 'Sample.ID')
}

survival_matrix <- prepare_survival_matrix(expression_data, survival_data)
error_variables <- names(survival_matrix)[seq_len(ncol(survival_matrix) - 2)]

# ---------------------------------------------------------------------------
# Cutpoint calculation
# ---------------------------------------------------------------------------
cutpoint_result <- surv_cutpoint(
  survival_matrix,
  time = 'OS.day',
  event = 'OS',
  variables = error_variables
)
summary(cutpoint_result)
categorized_data <- surv_categorize(cutpoint_result)

# ---------------------------------------------------------------------------
# Helper functions for survival tests
# ---------------------------------------------------------------------------
run_cox_regression <- function(data_frame, variables) {
  result <- map_dfr(variables, function(variable_name) {
    fit <- coxph(Surv(OS.day, OS) ~ get(variable_name), data = data_frame)
    fit_summary <- summary(fit)

    tibble(
      Error = variable_name,
      HR = fit_summary$coefficients[1, 2],
      P = fit_summary$coefficients[1, 5]
    )
  })

  arrange(result, P)
}

run_log_rank_test <- function(data_frame, variables) {
  p_values <- map_dbl(variables, function(variable_name) {
    fit <- survfit(Surv(OS.day, OS) ~ get(variable_name), data = data_frame)
    surv_pvalue(fit)[1, 2]
  })

  tibble(Error = variables, KM.cutpoint.P = p_values) %>%
    arrange(KM.cutpoint.P)
}

# ---------------------------------------------------------------------------
# Survival analyses
# ---------------------------------------------------------------------------
cox_results <- run_cox_regression(categorized_data, error_variables)
write.table(
  cox_results,
  'survival_cox.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

log_rank_results <- run_log_rank_test(categorized_data, error_variables)
write.table(
  log_rank_results,
  'survival.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

cat('Survival analysis completed.\n')
