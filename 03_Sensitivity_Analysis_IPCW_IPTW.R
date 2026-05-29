# ==============================================================================
# Script Name: 03_Sensitivity_Analysis_IPCW_IPTW.R
# Description: Independent Sensitivity Analyses for Retrospective Cohort Study 
#              (IPTW, IPCW, and Competing Risk Models).
# ==============================================================================

# Note: The raw dataset 'combined_data' and the baseline matched cohort 'matched_df' 
# utilized in this script are restricted due to patient privacy regulations governed 
# by the Shanghai Hospital Link Database (SHLD). This script serves as a methodological 
# record and assumes these datasets have been securely pre-loaded in the local 
# environment (e.g., 'matched_df' generated from Script 01).

library(dplyr)
library(survival)
library(WeightIt)
library(cmprsk)

# 1. Independent Stabilized IPTW Model
set.seed(12345)

iptw_model <- weightit(
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
  method = "ps",
  estimand = "ATE",
  stabilize = TRUE
)

iptw_df <- combined_data %>%
  mutate(
    weights_iptw = iptw_model$weights,
    first_date = as.Date(first_date),
    Endpoint_phychia = as.Date(Endpoint_phychia),
    S1_followup_end = as.Date(S1_followup_end),
    event = ifelse(!is.na(Endpoint_phychia) & Endpoint_phychia <= S1_followup_end, 1, 0),
    S1_followup_years = round(S1_followup_dur / 365.25, 3),
    treatment = factor(treatment, levels = c("DPP4i", "GLP1RA"))
  )

cox_iptw <- coxph(Surv(S1_followup_years, event) ~ treatment, data = iptw_df, weights = weights_iptw, robust = TRUE)
summary(cox_iptw)

# 2. Independent IPCW Model (Adjusting for drug discontinuation or switching)
# Note: Applied to the baseline matched cohort (matched_df) from Script 01
ipcw_df <- matched_df %>%
  mutate(
    first_date = as.Date(first_date),
    S1_drugduration = as.Date(S1_drugduration),
    drug_shift_time = as.Date(drug_shift_time),
    Endpoint_phychia = as.Date(Endpoint_phychia),
    S1_followup_end = as.Date(S1_followup_end),
    event = ifelse(!is.na(Endpoint_phychia) & Endpoint_phychia <= S1_followup_end, 1, 0),
    S1_followup_years = round(S1_followup_dur / 365.25, 3),
    treatment = factor(treatment, levels = c("DPP4i", "GLP1RA")),
    censored_informative = ifelse(event == 0 & 
                                    (S1_followup_end == S1_drugduration | S1_followup_end == drug_shift_time) &
                                    S1_followup_end < as.Date("2024-12-31"), 1, 0)
  )

set.seed(12345)

ipcw_model <- weightit(
  censored_informative ~ treatment + Age + Gender + EntryYear + A1c_new + co_renal + 
    co_neuro + co_hypertension + co_fat + hospital_all, 
  data = ipcw_df,
  method = "ps",
  estimand = "ATE",
  stabilize = TRUE,
  is.censoring = TRUE
)

ipcw_df <- ipcw_df %>%
  mutate(weights_ipcw = ipcw_model$weights)

cox_ipcw <- coxph(Surv(S1_followup_years, event) ~ treatment, data = ipcw_df, weights = weights_ipcw, robust = TRUE)
summary(cox_ipcw)

# 3. Competing Risk Model (Fine-Gray using cmprsk)
cmprsk_df <- matched_df %>%
  mutate(
    Endpoint_phychia = as.Date(Endpoint_phychia),
    swrq = as.Date(swrq),
    S1_followup_end = as.Date(S1_followup_end),
    S1_followup_years = round(S1_followup_dur / 365.25, 3),
    status_cr = case_when(
      !is.na(Endpoint_phychia) & Endpoint_phychia <= S1_followup_end ~ 1, 
      !is.na(swrq) & swrq <= S1_followup_end ~ 2,                         
      TRUE ~ 0                                                            
    ),
    treatment = factor(treatment, levels = c("DPP4i", "GLP1RA"))
  )

cov_matrix <- model.matrix(~ treatment, data = cmprsk_df)[, -1, drop = FALSE]

fg_model <- crr(
  ftime = cmprsk_df$S1_followup_years,
  fstatus = cmprsk_df$status_cr,
  cov1 = cov_matrix,
  failcode = 1,
  cencode = 0
)
summary(fg_model)