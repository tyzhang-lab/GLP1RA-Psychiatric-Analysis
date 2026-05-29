# ==============================================================================
# Script Name: 01_PSM_and_Primary_Analysis.R
# Description: Propensity Score Matching (PSM) and Primary Cox Proportional 
#              Hazards Models for On-Treatment and Post-Treatment Periods.
# ==============================================================================

# Note: The raw dataset 'combined_data' utilized in this study is restricted due 
# to patient privacy regulations governed by the Shanghai Hospital Link Database (SHLD). 
# This script serves as a methodological record and assumes the de-identified 
# dataset has been securely pre-loaded in the local environment.

# 1. Propensity Score Matching (1:1)
library(dplyr)
library(MatchIt)
library(survival)

set.seed(12345)

psm_model <- matchit(
  treatment ~ Age + Gender + EntryYear + A1c_new + co_renal + co_neuro + 
    co_ophthalmic + co_peripheral + co_ketoacidosis + co_myocardial_infarction + 
    co_ischemic_stroke + co_epilepsy + co_pain + co_hypertension + co_cancer + 
    co_fat + co_cog_impair + co_eat_disorder + co_sleep_disorder + drug_metformin + 
    drug_insulin + drug_sulfonylureas + drug_glinide + drug_Glyco_i + drug_TZD + 
    drug_GKA + drug_PPAR + drug_SGLT2i + drug_ACEI + drug_ARB + drug_betaR + 
    drug_CCB + drug_diuretic + drug_antiplatelet + drug_PPI + drug_statins + 
    drug_fibrates + drug_NSAIDs + drug_immuno + drug_opiates + drug_antiepileptic + 
    drug_orlistat + hospital_all + outpatient,
  data = combined_data,
  method = "nearest",
  caliper = 0.1,
  ratio = 1
)

matched_df <- match.data(psm_model)

# 2. Primary Analysis: On-Treatment Period
matched_df_S1 <- matched_df %>%
  mutate(
    first_date = as.Date(first_date),
    Endpoint_phychia = as.Date(Endpoint_phychia),
    S1_followup_end = as.Date(S1_followup_end),
    event = ifelse(!is.na(Endpoint_phychia) & Endpoint_phychia <= S1_followup_end, 1, 0),
    S1_followup_years = round(S1_followup_dur / 365.25, 3),
    treatment = factor(treatment, levels = c("DPP4i", "GLP1RA"))
  )

cox_on_treatment <- coxph(Surv(S1_followup_years, event) ~ treatment, data = matched_df_S1, robust = TRUE)
summary(cox_on_treatment)

# 3. Primary Analysis: Post-Treatment Period
matched_df_S2 <- matched_df %>%
  mutate(
    swrq = as.Date(swrq), 
    Endpoint_phychia = as.Date(Endpoint_phychia),
    drug_shift_time = as.Date(drug_shift_time),
    S1_drugduration = as.Date(S1_drugduration)
  ) %>%
  filter(S1_drugduration < as.Date("2024-12-31")) %>%
  filter(is.na(drug_shift_time) | S1_drugduration < drug_shift_time) %>%
  filter(is.na(swrq) | S1_drugduration < swrq) %>% 
  filter(is.na(Endpoint_phychia) | S1_drugduration < Endpoint_phychia) %>%
  mutate(
    S2_start = S1_drugduration,
    S2_followup_end = pmin(as.Date("2024-12-31"), swrq, Endpoint_phychia, drug_shift_time, na.rm = TRUE),
    S2_followup_days = as.numeric(difftime(S2_followup_end, S2_start, units = "days"))
  ) %>%
  filter(S2_followup_days > 0) %>%
  mutate(
    event2 = ifelse(!is.na(Endpoint_phychia) & Endpoint_phychia <= S2_followup_end, 1, 0),
    S2_followup_years = round(S2_followup_days / 365.25, 3),
    treatment = factor(treatment, levels = c("DPP4i", "GLP1RA"))
  )

cox_post_treatment <- coxph(Surv(S2_followup_years, event2) ~ treatment, data = matched_df_S2, robust = TRUE)
summary(cox_post_treatment)