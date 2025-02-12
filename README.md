
<!-- README.md is generated from README.Rmd. Please edit that file -->


# PSweight

This `PSweight` Package is to perform propensity score weighting
analysis for causal inference. Two main modules are included to assist
the design and analysis of observational studies. In the design module,
the `SumStat` function is used to generate distributional plots of the
estimated propensity scores and balance diagnostics after propensity
score weighting. The `summary` and `plot` functions are available to
tabulate and plot weighted balance statistics for visual comparisons. In
the analysis module, the `PSweight` function, the average potential
outcomes for each treatment group is estimated using weighting, and the
`summary` function generates point estimates, standard errors and
confidence intervals for the desired causal contrasts of interest. The
current version of `PSweight` package includes the following types of
weights: the overlap weights (ATO), the inverse probability of treatment
weights (ATE), the average treatment effect among the treated weights
(ATT), the matching weights (ATM) and the entropy weights (ATEN), and
allows for binary and multiple (categorical) treatments. In addition to
the simple weighting estimator, the package also implements the
augmented weighting estimator that combines weighting and outcome
regression. For binary outcomes, both the additive and ratio estimands
(causal relative risk and odds ratio) are considered, and variance is
estimated by either the sandwich method or nonparametric bootstrap. To
allow for additional flexibility in specifying the propensity score and
outcome models, the package can also work with user-supplied propensity
score estimates and outcome predictions through `ps.estimate` and
`out.estimate`, and provide a sandwich standard error that ignores the
variability in estimating these nuisances. Additionally, `PSweight` package has been extended its application to survey observational data, providing comprehensive implementations for propensity score weighting and augmented estimators under common survey settings, making it a powerful tool for causal inference in complex survey datasets.


## Installation

You can install the released version of PSweight from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("PSweight")
```

## Maintainer

yukang.zeng@yale.edu

## Update

Version 2.1.0 of PSweight extends its application to survey observational data, enabling its use in real-world studies that incorporate complex survey designs. This version introduces several enhancements to facilitate population-level inference under survey settings.

1.The PSmethod and PStrim functions have been expanded to incorporate survey weights for population-level propensity score estimation. This update allows for the estimation of propensity scores using a survey-weighted regression model, ensuring that population-level effects are appropriately accounted for in the analysis.

2.Additionally, the new release introduces survey-based propensity score weighting estimators and extends the augmented estimators to survey settings. Three augmented estimators are now implemented: Moment Estimator (MOM), Clever Covariate Estimator (CVR), and Weighted Regression Estimator (WET, default). These estimators are applicable under three common survey sampling settings, including Retrospective, where sampling depends on both treatment assignment and covariates; Prospective, where sampling occurs before treatment assignment; and Independent, where the sampling process is independent of treatment assignment.

3.Enhancements have also been made to balance assessment and propensity score visualization functions. The SumStat function and related diagnostic tools have been extended to incorporate survey-weighted analysis, and a new function, SumStat_SW, has been introduced to support survey-specific balance diagnostics and propensity score visualization. These updates ensure that users can assess covariate balance and evaluate the effectiveness of propensity score weighting within survey-based observational studies.

4.Certain limitations remain in this release. Currently, multigroup and cluster-based propensity score weighting functionalities have not been extended to survey data settings and remain unchanged. Furthermore, external propensity score estimates (ps.estimate) are not supported under survey settings to ensure the integrity of population-level estimation. Bootstrap variance estimation is also not supported for survey observational data, and users are encouraged to use sandwich variance estimators instead for variance estimation.

These updates significantly enhance the applicability of PSweight for survey-based causal inference, enabling researchers to account for complex sampling mechanisms while maintaining robust estimation procedures. The theoretical framework behind these enhancements is detailed in Zeng, Li, and Tong (2025), which provides best practices for incorporating survey weights into propensity score weighting methods to enhance the applicability and accuracy of causal inference in survey observational studies.

## Downloads

![Downloads
Status](https://cranlogs.r-pkg.org/badges/grand-total/PSweight)

## Example

This is a basic example on design:

``` r
library(PSweight)
#> Warning: replacing previous import 'lifecycle::last_warnings' by
#> 'rlang::last_warnings' when loading 'tibble'
#> Warning: replacing previous import 'lifecycle::last_warnings' by
#> 'rlang::last_warnings' when loading 'pillar'
example("SumStat")
#> 
#> SumStt> data("psdata")
#> 
#> SumStt> # the propensity model
#> SumStt> ps.formula<-trt~cov1+cov2+cov3+cov4+cov5+cov6
#> 
#> SumStt> # using SumStat to estimate propensity scores
#> SumStt> msstat <- SumStat(ps.formula, trtgrp="2", data=psdata,
#> SumStt+    weight=c("IPW","overlap","treated","entropy","matching"))
#> 
#> SumStt> #summary(msstat)
#> SumStt> 
#> SumStt> # importing user-supplied propensity scores "e.h"
#> SumStt> # fit <- nnet::multinom(formula=ps.formula, data=psdata, maxit=500, trace=FALSE)
#> SumStt> # e.h <- fit$fitted.values
#> SumStt> # varname <- c("cov1","cov2","cov3","cov4","cov5","cov6")
#> SumStt> # msstat0 <- SumStat(zname="trt", xname=varname, data=psdata, ps.estimate=e.h,
#> SumStt> #  trtgrp="2",  weight=c("IPW","overlap","treated","entropy","matching"))
#> SumStt> # summary(msstat0)
#> SumStt> 
#> SumStt> 
#> SumStt> 
#> SumStt>
```

This is a basic example on analysis:

``` r
example("PSweight")
#> 
#> PSwght> data("psdata")
#> 
#> PSwght> # the propensity and outcome models
#> PSwght> ps.formula<-trt~cov1+cov2+cov3+cov4+cov5+cov6
#> 
#> PSwght> out.formula<-Y~cov1+cov2+cov3+cov4+cov5+cov6
#> 
#> PSwght> # without augmentation
#> PSwght> ato1<-PSweight(ps.formula = ps.formula,yname = 'Y',data = psdata,weight = 'overlap')
#> 
#> PSwght> summary(ato1)
#> 
#> Closed-form inference: 
#> 
#> Original group value:  1, 2, 3 
#> 
#> Contrast: 
#>             1  2 3
#> Contrast 1 -1  1 0
#> Contrast 2 -1  0 1
#> Contrast 3  0 -1 1
#> 
#>            Estimate Std.Error      lwr      upr  Pr(>|z|)    
#> Contrast 1 -1.24161   0.16734 -1.56960 -0.91362 1.177e-13 ***
#> Contrast 2  1.12482   0.17099  0.78968  1.45996 4.764e-11 ***
#> Contrast 3  2.36643   0.25854  1.85970  2.87315 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> PSwght> # augmented weighting estimator, takes longer time to calculate sandwich variance
#> PSwght> # ato2<-PSweight(ps.formula = ps.formula,yname = 'Y',data = psdata,
#> PSwght> #              augmentation = TRUE,out.formula = out.formula,family = 'gaussian',weight = 'overlap')
#> PSwght> # summary(ato2)
#> PSwght> 
#> PSwght> 
#> PSwght> 
#> PSwght>
```
