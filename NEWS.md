
\<!-NEWS.md is generated from NEWS.Rmd. Please edit that file –\>

# PSweight 2.1.1 (Release Date: 2025-03-01)

Version 2.1.0 of PSweight extends its application to survey observational data, enabling its use in real-world studies that incorporate complex survey designs. This version introduces several enhancements to facilitate population-level inference under survey settings.

1.The PSmethod and PStrim functions have been expanded to incorporate survey weights for population-level propensity score estimation. This update allows for the estimation of propensity scores using a survey-weighted regression model, ensuring that population-level effects are appropriately accounted for in the analysis.

2.Additionally, the new release introduces survey-based propensity score weighting estimators and extends the augmented estimators to survey settings. Three augmented estimators are now implemented: Moment Estimator (MOM), Clever Covariate Estimator (CVR), and Weighted Regression Estimator (WET, default). These estimators are applicable under three common survey sampling settings, including Retrospective, where sampling depends on both treatment assignment and covariates; Prospective, where sampling occurs before treatment assignment; and Independent, where the sampling process is independent of treatment assignment.

3.Enhancements have also been made to balance assessment and propensity score visualization functions. The SumStat function and related diagnostic tools have been extended to incorporate survey-weighted analysis, and a new function, SumStat_SW, has been introduced to support survey-specific balance diagnostics and propensity score visualization. These updates ensure that users can assess covariate balance and evaluate the effectiveness of propensity score weighting within survey-based observational studies.

4.Certain limitations remain in this release. Currently, multigroup and cluster-based propensity score weighting functionalities have not been extended to survey data settings and remain unchanged. Furthermore, external propensity score estimates (ps.estimate) are not supported under survey settings to ensure the integrity of population-level estimation. Bootstrap variance estimation is also not supported for survey observational data, and users are encouraged to use sandwich variance estimators instead for variance estimation.

These updates significantly enhance the applicability of PSweight for survey-based causal inference, enabling researchers to account for complex sampling mechanisms while maintaining robust estimation procedures. The theoretical framework behind these enhancements is detailed in Zeng, Li, and Tong (2025), which provides best practices for incorporating survey weights into propensity score weighting methods to enhance the applicability and accuracy of causal inference in survey observational studies.

Version 2.1.1 adds additional examples with more detailed output for each core function to better demonstrate usage under survey settings.