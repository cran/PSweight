% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SumStat_SW.R
\name{SumStat_SW}
\alias{SumStat_SW}
\title{Calculate summary statistics for propensity score weighting}
\usage{
SumStat_SW(
  ps.formula = NULL,
  ps.estimate = NULL,
  trtgrp = NULL,
  Z = NULL,
  covM = NULL,
  zname = NULL,
  xname = NULL,
  data = NULL,
  weight = "overlap",
  survey.indicator = FALSE,
  survey.design = "Independent",
  svywtname = NULL,
  delta = 0,
  method = "glm",
  ps.control = list()
)
}
\arguments{
\item{ps.formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the propensity score model to be fitted. Additional details of model specification are given under "Details". This argument is optional if \code{ps.estimate} is not \code{NULL}.}

\item{ps.estimate}{an optional matrix or data frame containing estimated (generalized) propensity scores for each observation. Typically, this is an N by J matrix, where N is the number of observations and J is the total number of treatment levels. Preferably, the column names of this matrix should match the names of treatment level, if column names are missing or there is a mismatch, the column names would be assigned according to the alphabatic order of treatment levels. A vector of propensity score estimates is also allowed in \code{ps.estimate}, in which case a binary treatment is implied and the input is regarded as the propensity to receive the last category of treatment by alphabatic order, unless otherwise stated by \code{trtgrp}.}

\item{trtgrp}{an optional character defining the "treated" population for estimating the average treatment effect among the treated (ATT). Only necessary if \code{weight = "treated"}. This option can also be used to specify the treatment (in a two-treatment setting) when a vector argument is supplied for \code{ps.estimate}. Default value is the last group in the alphebatic order.}

\item{Z}{an optional vector specifying the values of treatment, only necessary when the covariate matrix \code{covM} is provided instead of \code{data}.}

\item{covM}{an optional covariate matrix or data frame including covariates, their interactions and higher-order terms. When the covariate matrix \code{covM} is provided, the balance statistics are generated according to each column of this matrix.}

\item{zname}{an optional character specifying the name of the treatment variable in \code{data}.}

\item{xname}{an optional vector of characters including the names of covariates in \code{data}.}

\item{data}{an optional data frame containing the variables in the propensity score model. If not found in data, the variables are taken from \code{environment(formula)}.}

\item{weight}{a character or vector of characters including the types of weights to be used. \code{"IPW"} specifies the inverse probability weights for estimating the average treatment effect among the combined population (ATE). \code{"treated"} specifies the weights for estimating the average treatment effect among the treated (ATT). \code{"overlap"} specifies the (generalized) overlap weights for estimating the average treatment effect among the overlap population (ATO), or population at clinical equipoise. \code{"matching"} specifies the matching weights for estimating the average treatment effect among the matched population (ATM). \code{"entropy"} specifies the entropy weights for the average treatment effect of entropy weighted population (ATEN). Default is \code{"overlap"}.}

\item{survey.indicator}{logical. Indicates whether survey weights are used in the estimation. 
Default is \code{FALSE}. If \code{TRUE}, survey weights specified in \code{svywtname} will be used.}

\item{survey.design}{character. Specifies the survey design scenario for estimation. 
Acceptable values are "Retrospective", "Independent", and "Prospective". 
"Retrospective" indicates that the sampling process depends on both treatment assignment and covariates, 
"Independent" (the default) means that the sampling process is independent of treatment assignment, 
and "Prospective" signifies that sampling is conducted prior to treatment assignment, 
although treatment may later be influenced by the sampling process.}

\item{svywtname}{an optional character specifying the name of the survey weight variable in \code{data}.
Default is \code{NULL}. Only required if \code{survey.indicator} is \code{TRUE}. If \code{survey.indicator} is \code{TRUE}
and \code{svywtname} is not provided, a default survey weight of 1 will be applied to all samples.}

\item{delta}{trimming threshold for estimated (generalized) propensity scores. Should be no larger than 1 / number of treatment groups. Default is 0, corresponding to no trimming.}

\item{method}{a character to specify the method for estimating propensity scores. \code{"glm"} is default, and \code{"gbm"} and \code{"SuperLearner"} are also allowed.}

\item{ps.control}{a list to specify additional options when \code{method} is set to \code{"gbm"} or \code{"SuperLearner"}.}
}
\value{
SumStat returns a \code{SumStat} object including a list of the following value:
treatment group, propensity scores, fitted propensity model, propensity score weights, effective sample sizes,
and balance statistics. A summary of \code{SumStat} can be obtained with \code{\link{summary.SumStat}}.

\describe{
\item{\code{ trtgrp}}{a character indicating the treatment group.}

\item{\code{ propensity}}{a data frame of estimated population level propensity scores. When \code{survey.indicator = TRUE}, it is population level propensity score estimated by survey-weighted regression.}

\item{\code{ ps.population.fitObjects}}{the fitted population level propensity model details}

\item{\code{ propensity.sample}}{a data frame of estimated sample level propensity scores.}

\item{\code{ ps.sample.fitObjects}}{the fitted sample level propensity model details}

\item{\code{ ps.weights}}{a data frame of propensity score weights.}

\item{\code{ ess}}{a table of effective sample sizes. This serves as a conservative measure to
characterize the variance inflation or precision loss due to weighting, see Li and Li (2019).}

\item{\code{ unweighted.sumstat}}{A list of tables including covariate means and variances
by treatment group and standardized mean differences.}

\item{\code{ ATE.sumstat}}{If \code{"IPW"} is included in \code{weight}, this is a list of summary statistics using inverse probability weighting.}

\item{\code{ ATT.sumstat}}{If \code{"treated"} is included in \code{weight}, this is a list of summary statistics using the ATT weights.}

\item{\code{ ATO.sumstat}}{If \code{"overlap"} is included in \code{weight}, this is a list of summary statistics using the overlap weights.}

\item{\code{ ATM.sumstat}}{If \code{"matching"} is included in \code{weight}, this is a list of summary statistics using the matching weights.}

\item{\code{ ATEN.sumstat}}{If \code{"entropy"} is included in \code{weight}, this is a list of summary statistics using the entropy weights.}

\item{\code{ trim}}{If \code{delta > 0}, this is a table summarizing the number of observations before and after trimming.}

}
}
\description{
\code{SumStat} is used to generate distributional plots of the estimated propensity scores and balance
diagnostics after propensity score weighting.
}
\details{
A typical form for \code{ps.formula} is \code{treatment ~ terms} where \code{treatment} is the treatment
variable (identical to the variable name used to specify \code{zname}) and \code{terms} is a series of terms
which specifies a linear predictor for \code{treatment}. \code{ps.formula} specifies logistic or multinomial logistic
models for estimating the propensity scores, when \code{ps.estimate} is \code{NULL}.
The function now includes support for the survey setting, mainly focusing on binary treatments by incorporating survey weights into propensity score estimation and both point and augmented estimators for the outcome stage. In survey settings,
external propensity score estimates are not supported; both population-level and sample-level propensity scores are estimated using the internal routines.

When comparing two treatments, \code{ps.estimate} can either be a vector or a two-column matrix of estimated
propensity scores. If a vector is supplied, it is assumed to be the propensity scores to receive the treatment, and
the treatment group corresponds to the last group in the alphebatic order, unless otherwise specified by \code{trtgrp}.
When comparing multiple (J>=3) treatments, \code{ps.estimate} needs to be specified as an N by J matrix,
where N indicates the number of observations, and J indicates the total number of treatments.
This matrix specifies the estimated generalized propensity scores to receive each of the J treatments.
In general, \code{ps.estimate} should have column names that indicate the level of the treatment variable,
which should match the levels given in \code{Z}.
If column names are empty or there is a mismatch, the column names will be created following
the alphebatic order of treatmentlevels. The rightmost coulmn of \code{ps.estimate} is then assumed
to be the treatment group when estimating ATT (\code{"treated"}). \code{trtgrp} can also be used to specify the treatment
group for estimating ATT.

To generate balance statistics, one can directly specify \code{Z} and \code{covM} to indicate the treatment levels and
covariate matrix. Alternatively, one can supply \code{data}, \code{zname}, and \code{xname} to indicate the
same information. When both are specified, the function will prioritize inputs from \code{Z} and \code{covM}.
When \code{ps.estimate} is not \code{NULL}, argument \code{zname}.

In survey settings, when \code{survey.indicator} is \code{TRUE}, the argument \code{svywtname} (which specifies the survey weight variable in \code{data}) is required; 
if \code{svywtname} is not provided, a default survey weight of 1 is applied to all observations. The argument 
\code{survey.design} must be specified to reflect the sampling mechanism: for example, \code{"Retrospective"} indicates 
that the sampling process depends on both treatment assignment and covariates, \code{"Independent"} assumes that sampling 
is independent of treatment assignment, and \code{"Prospective"} signifies that sampling is conducted prior to treatment assignment, 
although treatment may later be influenced by the sampling results.

Current version of \code{PSweight} allows for five types of propensity score weights used to estimate ATE (\code{"IPW"}), ATT (\code{"treated"}), and
ATO(\code{"overlap"}), ATM (\code{"matching"}) and ATEN (\code{"entropy"}). These weights are members of a larger class of balancing weights defined in Li, Morgan, and Zaslavsky (2018).
When there is a practical violation of the positivity assumption, \code{delta} defines the symmetric
propensity score trimming rule following Crump et al. (2009). With multiple treatments, \code{delta} defines the
multinomial trimming rule introduced in Yoshida et al. (2019). The overlap weights can also be considered as
a data-driven continuous trimming strategy without specifying trimming rules, see Li, Thomas and Li (2019).
Additional details on balancing weights and generalized overlap weights for multiple treatment groups are provided in
Li and Li (2019). For details about matching weights and entropy weights, please refer to Li and Greene (2013) and Zhou, Matsouaka and Thomas (2020).
Zeng, Li, and Tong (2025) further specify how the survey weights can be incorporated into propensity score weighting under both retrospective and prospective scenarios.
Their approach supports both a weighting-only estimator and all three augmented estimators(MOM, CVR and WET), with corresponding sandwich variance estimators developed.
These enhancements are implemented in the current version of \code{PSweight}.


\code{"glm"} is the default method for propensity score estimation. Logistic regression will be used for binary outcomes,
and multinomial logistic regression will be used for outcomes with more than two categories. The alternative method option of \code{"gbm"} serves as an API to call the \code{gbm()} function from the
\code{gbm} package. Additional argument in the \code{gbm()} function can be supplied through the \code{ps.control=list()} argument in \code{SumStat()}. Please refer to the user manual of the gbm package for all the
allowed arguments. Currently, models for binary or multinomial treatment will be automatically chosen based on the number of treatment categories.
\code{"SuperLearner"} is also allowed in the \code{method} argument to pass the propensity score estimation to the \code{SuperLearner()} function in SuperLearner package.
Currently, the SuperLearner method only supports binary treatment with the default method set to \code{"SL.glm"}. The estimation approach is default to \code{"method.NNLS"} in the \code{SumStat()} function.
Prediction algorithm and other tuning parameters can also be passed through \code{ps.control=list()} to \code{SumStat()}. Please refer to the user manual of the \code{SuperLearner} package for all the allowed specifications.
}
\examples{

data("psdata")
data("psdata_bin_prospective_fp")
data("psdata_bin_retrospective_fp")

# Define the common propensity score formula
ps.formula <- trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6

# The weight choices we want to test:
myWeights <- c("IPW", "overlap", "treated", "entropy", "matching", "svywt")

## Example 1: Prospective design, survey indicator TRUE, using all weight types
sumstat_pros <- SumStat_SW(ps.formula = ps.formula, trtgrp = "2", 
                           data = psdata_bin_prospective_fp,
                           weight = myWeights, delta = 0.1, method = "glm", ps.control = list(),
                           survey.indicator = TRUE, svywtname = "survey_weight", 
                           survey.design = "Prospective")
summary(sumstat_pros)
print(sumstat_pros$ps.weights)
print(sumstat_pros$ess)
print(sumstat_pros$unweighted.sumstat)

## Example 2: Retrospective design, survey indicator TRUE, using all weight types
sumstat_retro <- SumStat_SW(ps.formula = ps.formula, trtgrp = "2", 
                            data = psdata_bin_retrospective_fp,
                            weight = myWeights, delta = 0.1, method = "glm", ps.control = list(),
                            survey.indicator = TRUE, svywtname = "survey_weight", 
                            survey.design = "Retrospective")
summary(sumstat_retro)
print(sumstat_retro$ps.weights)
print(sumstat_retro$ess)
print(sumstat_retro$unweighted.sumstat)

}
\references{
Crump, R. K., Hotz, V. J., Imbens, G. W., Mitnik, O. A. (2009).
Dealing with limited overlap in estimation of average treatment effects. Biometrika, 96(1), 187-199.

Greenwell B., Boehmke B.,Cunningham J, GBM Developers (2020) gbm: Generalized Boosted Regression Models. Cran: https://cran.r-project.org/web/packages/gbm/index.html

Li, L., Greene, T. (2013).
A weighting analogue to pair matching in propensity score analysis. The International Journal of Biostatistics, 9(2), 215-234.

Li, F., Morgan, K. L., Zaslavsky, A. M. (2018).
Balancing covariates via propensity score weighting.
Journal of the American Statistical Association, 113(521), 390-400.

Li, F., Thomas, L. E., Li, F. (2019).
Addressing extreme propensity scores via the overlap weights. American Journal of Epidemiology, 188(1), 250-257.

Polley E., LeDell E., Kennedy C., Lendle S., van der Laan M. (2019) SuperLearner: Super Learner Prediction. Cran: https://cran.r-project.org/web/packages/SuperLearner/index.html

Yoshida, K., Solomon, D.H., Haneuse, S., Kim, S.C., Patorno, E., Tedeschi, S.K., Lyu, H.,
Franklin, J.M., Stürmer, T., Hernández-Díaz, S. and Glynn, R.J. (2019).
Multinomial extension of propensity score trimming methods: A simulation study.
American Journal of Epidemiology, 188(3), 609-616.

Li, F., Li, F. (2019). Propensity score weighting for causal inference with multiple treatments.
The Annals of Applied Statistics, 13(4), 2389-2415.

Zhou, Y., Matsouaka, R. A., Thomas, L. (2020).
Propensity score weighting under limited overlap and model misspecification. Statistical Methods in Medical Research (Online)

Zeng, Y., Li, F., & Tong, G. (2025). 
Moving toward best practice when using propensity score weighting in survey observational studies.
arXiv preprint arXiv:2501.16156.
}
