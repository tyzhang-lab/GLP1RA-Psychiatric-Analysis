# Risk of Depressive and Anxiety Disorders Following Discontinuation of GLP-1 Receptor Agonists in Type 2 Diabetes: A Retrospective Cohort Study

This repository contains the standardized R analytical code and analytical framework for the replication of the retrospective active-controlled cohort study described in our paper.

---

## Project Overview
Our study established a retrospective active-controlled cohort design to investigate the association between glucagon-like peptide-1 receptor agonist (GLP-1RA) therapy and the risk of depressive and anxiety disorders using a large-scale electronic health record database.

Within this cohort study framework:
- Primary Analysis (As-Treated): Compares psychiatric event risks among GLP-1RA initiators and active comparators (DPP-4is or SGLT-2is) according to actual drug exposure during the active on-treatment period.
- Post-Treatment Analysis: Evaluates outcomes during the post-discontinuation window to detect potential psychiatric event variations following treatment cessation.
- Sensitivity analyses: Includes an Intention-to-Treat (ITT) analysis, Inverse Probability of Treatment Weighting (IPTW), Inverse Probability of Censoring Weighting (IPCW), Fine-Gray competing risk models (accounting for mortality), and a Re-anchored Propensity Score Matching at the time of medication discontinuation to isolate legacy effects.
- Mechanistic Evaluation: Applies causal mediation analysis to evaluate the potential mediating effects of post-treatment metabolic rebound.

---

## Repository Structure & Scripts

The analytical pipeline is modularized into four core executable R scripts. For reproducibility, it is highly recommended to review or run them sequentially:

- 01_PSM_and_Primary_Analysis.R
  Performs the baseline 1:1 propensity score matching (PSM), executes primary Cox proportional hazards regression models with robust standard errors, and generates data frames for Kaplan-Meier estimates for both the on-treatment and post-treatment phases.

- 02_Sensitivity_Analysis_Reanchored.R
  Constructs the re-anchored post-treatment cohort, moving Time Zero to the exact point of treatment discontinuation, re-assesses time-varying covariates, and fits re-anchored Cox survival models.

- 03_Sensitivity_Analysis_IPCW_IPTW.R
  Implements alternative robust sensitivity measures including stabilized IPTW, informative censoring adjustments via IPCW, and Fine-Gray subdistribution hazard models to address competing risks of mortality.

- 04_Mechanistic_Mediation.R
  Constructs the casual mediation modeling function using parametric Weibull survival links to calculate Average Causal Mediation Effects (ACME) for lipid trajectory fluctuations.

---

## Data Notice
Raw electronic health records from the SHLD are restricted due to patient privacy regulations and are not hosted in this repository. Please refer to the "Data Availability" section in the manuscript for data access protocols and the provided Source Data files.

---

## Environment & Dependencies
All mathematical tracking and modeling were executed using the open-source platform R software (version 4.5.1).

To execute the scripts in this repository, ensure that the following core libraries and dependencies are pre-installed:
- tidyverse (v2.0.0): For data wrangling, preprocessing, and tabular harmonization.
- MatchIt (v4.7.2): For propensity score estimation, caliper matching, and balance diagnostics.
- survival (v3.8-3): For fitting Cox proportional hazards, estimating hazard ratios, and handling right-censored survival outcomes.
- survminer (v0.5.0): For extracting survival plot data tables.
- cmprsk (v3.1): For calculating cumulative incidence under a multi-state competing risk framework (Fine-Gray model).
- WeightIt (v1.2.0): For deriving stabilized weights for IPTW and IPCW modeling.
- mediation (v4.5.0): For causal mediation computations and bootstrapping of indirect mediation proportions.
