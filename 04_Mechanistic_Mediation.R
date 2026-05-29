# ==============================================================================
# Script Name: 04_Mechanistic_Mediation.R
# Description: Causal Mediation Analysis to evaluate the mediating effect of 
#              post-treatment lipid rebound on psychiatric outcomes.
# ==============================================================================

# Note: The raw dataset 'final_cohort_df' utilized in this script is restricted 
# due to patient privacy regulations governed by the Shanghai Hospital Link 
# Database (SHLD). This script serves as a methodological record and assumes 
# the de-identified dataset has been securely pre-loaded in the local environment.

# 1. Define Causal Mediation Function
library(dplyr)
library(survival)
library(mediation)

analyze_marker_mediation <- function(data, marker_prefix, sims = 1000, quantile_prob = 0.80) {
  
  col_P2 <- paste0(marker_prefix, "_P2")
  col_P3 <- paste0(marker_prefix, "_P3")
  
  df_clean <- data %>%
    filter(!is.na(!!sym(col_P2)) & !is.na(!!sym(col_P3))) %>%
    mutate(
      delta_lipid = !!sym(col_P3) - !!sym(col_P2),
      treat_num = ifelse(treatment == "GLP1RA", 1, 0),
      Gender = as.factor(Gender),
      Age = as.numeric(Age)
    ) %>%
    filter(!is.na(treat_num) & !is.na(Age) & !is.na(Gender) &
             !is.na(event2) & !is.na(S2_followup_years) & S2_followup_years > 0)
  
  if(nrow(df_clean) < 100) return(NULL)
  
  threshold_val <- quantile(df_clean$delta_lipid[df_clean$treat_num == 0], probs = quantile_prob, na.rm = TRUE)
  
  df_clean <- df_clean %>%
    mutate(is_high_rebound = ifelse(delta_lipid > threshold_val, 1, 0))
  
  formula_m <- as.formula(paste("is_high_rebound ~ treat_num +", col_P2, "+ Age + Gender"))
  fitM <- glm(formula_m, data = df_clean, family = binomial(link = "logit"))
  
  formula_y <- as.formula(paste("Surv(S2_followup_years, event2) ~ treat_num + is_high_rebound +", col_P2, "+ Age + Gender"))
  fitY <- survreg(formula_y, data = df_clean, dist = "weibull")
  
  set.seed(2026)
  med_out <- mediate(
    model.m = fitM,
    model.y = fitY,
    treat = "treat_num",
    mediator = "is_high_rebound",
    control.value = 0,
    treat.value = 1,
    sims = sims
  )
  
  res <- data.frame(
    Biomarker = marker_prefix,
    Total_Effect = med_out$tau.coef,
    Total_P = med_out$tau.p,
    Direct_Effect_ADE = med_out$z.avg,
    Direct_P = med_out$z.avg.p,
    Indirect_Effect_ACME = med_out$d.avg,
    Indirect_P = med_out$d.avg.p,
    Proportion_Mediated = med_out$n.avg
  )
  
  return(res)
}

# 2. Execute Analysis on Re-anchored Cohort
res_TC  <- analyze_marker_mediation(final_cohort_df, "TC", sims = 1000)
res_TG  <- analyze_marker_mediation(final_cohort_df, "TG", sims = 1000)
res_LDL <- analyze_marker_mediation(final_cohort_df, "LDL", sims = 1000)

mediation_summary <- bind_rows(res_TC, res_TG, res_LDL)
print(mediation_summary)