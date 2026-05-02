Project title: Risk of Depressive and Anxiety Disorders Following Discontinuation of GLP-1 Receptor Agonists in Type 2 Diabetes: A Target Trial Emulation

Our study applied a target trial emulation (TTE) (1, 2) framework to investigate the association between glucagon-like peptide-1 receptor agonist (GLP-1RA) therapy and the risk of depressive and anxiety disorders. Within this TTE framework, we conducted a primary as-treated analysis to compare psychiatric event risks among GLP-1RA initiators and patients starting active comparators (DPP-4is or SGLT-2is) according to actual drug exposure. We also performed an intention-to-treat (ITT) analysis as a sensitivity measure. Risks were estimated across distinct clinical phases by separating the active on-treatment period from the post-treatment period. As a further sensitivity analysis, a re-anchored propensity score matching was performed for the post-treatment period to isolate legacy effects following medication discontinuation.

All analyses were performed in R software (version 4.5.1).

The following R packages were used for the analysis:
MatchIt: Propensity score matching (PSM) for baseline covariate balance.
survival: Cox proportional hazards regression and Kaplan-Meier estimates.
cmprsk: Fine-Gray subdistribution hazard models for competing risk of mortality.
WeightIt: Inverse probability of treatment weighting (IPTW) and inverse probability of censoring weighting (IPCW).
mediation: Causal mediation analysis for evaluating post-treatment lipid rebound.

Reference:
1. Hubbard, R. A. et al. ‘Target Trial Emulation’ for Observational Studies - Potential and Pitfalls. N Engl J Med 391, 1975–1977 (2024).
2. Cashin, A. G. et al. Transparent Reporting of Observational Studies Emulating a Target Trial-The TARGET Statement. JAMA 334, 1084–1093 (2025).
