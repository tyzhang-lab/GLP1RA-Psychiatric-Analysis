# ==============================================================================
# Script Name: 02_Sensitivity_Analysis_Reanchored.R
# Description: Re-anchored Propensity Score Matching at the time of treatment 
#              discontinuation for Post-Treatment Analysis.
# ==============================================================================

# Note: The raw dataset 'combined_data' utilized in this study is restricted due 
# to patient privacy regulations governed by the Shanghai Hospital Link Database (SHLD). 
# This script serves as a methodological record and assumes the de-identified 
# dataset has been securely pre-loaded in the local environment.

# 1. Prepare Re-anchored Cohort at Discontinuation
library(dplyr)
library(MatchIt)
library(survival)

post_tx_cohort <- combined_data %>%
  mutate(
    swrq = as.Date(swrq), 
    Endpoint_phychia = as.Date(Endpoint_phychia),
    drug_shift_time = as.Date(drug_shift_time),
    S1_drugduration = as.Date(S1_drugduration)
  ) %>%
  filter(S1_drugduration < as.Date("2024-12-31")) %>%
  filter(is.na(drug_shift_time) | S1_drugduration < drug_shift_time) %>%
  filter(is.na(swrq) | S1_drugduration < swrq) %>% 
  filter(is.na(Endpoint_phychia) | S1_drugduration < Endpoint_phychia)

# 2. Re-anchored Propensity Score Matching
# Note: All covariates were re-assessed/updated using S1_drugduration as the new baseline (Time Zero).
set.seed(12345)

reanchored_psm <- matchit(
  treatment ~ Age + Gender + EntryYear + A1c_new + co_renal + co_neuro + 
    co_ophthalmic + co_peripheral + co_ketoacidosis + co_myocardial_infarction + 
    co_ischemic_stroke + co_epilepsy + co_pain + co_hypertension + co_cancer + 
    co_fat + co_cog_impair + co_eat_disorder + co_sleep_disorder + drug_metformin + 
    drug_insulin + drug_sulfonylureas + drug_glinide + drug_Glyco_i + drug_TZD + 
    drug_GKA + drug_PPAR + drug_SGLT2i + drug_ACEI + drug_ARB + drug_betaR + 
    drug_CCB + drug_diuretic + drug_antiplatelet + drug_PPI + drug_statins + 
    drug_fibrates + drug_NSAIDs + drug_immuno + drug_opiates + drug_antiepileptic + 
    drug_orlistat + hospital_all + outpatient + S1_followup_dur,
  data = post_tx_cohort,
  method = "nearest",
  caliper = 0.1,
  ratio = 1
)

matched_reanchored_df <- match.data(reanchored_psm)

# 3. Post-Treatment Survival Analysis
matched_reanchored_df <- matched_reanchored_df %>%
  mutate(
    S2_followup_end = pmin(as.Date("2024-12-31"), swrq, Endpoint_phychia, drug_shift_time, na.rm = TRUE),
    S2_followup_days = as.numeric(difftime(S2_followup_end, S1_drugduration, units = "days"))
  ) %>%
  filter(S2_followup_days > 0) %>%
  mutate(
    event2 = ifelse(!is.na(Endpoint_phychia) & Endpoint_phychia <= S2_followup_end, 1, 0),
    S2_followup_years = round(S2_followup_days / 365.25, 3),
    treatment = factor(treatment, levels = c("DPP4i", "GLP1RA"))
  )

cox_reanchored <- coxph(Surv(S2_followup_years, event2) ~ treatment, data = matched_reanchored_df, robust = TRUE)
summary(cox_reanchored)