#' Retrospective Superpopulation Dataset for PSweight Analysis
#'
#' A simulated dataset representing the retrospective superpopulation (i.e. the full dataset)
#' for propensity score weighting analysis under a retrospective design. In this design, the
#' sampling probability is influenced by both covariates and treatment assignment.
#'
#' @docType data
#'
#' @usage data(psdata_bin_retrospective_sp)
#'
#' @format A data frame with 1500 rows and 10 variables:
#' \describe{
#'   \item{Y}{Outcome variable.}
#'   \item{trt}{Binary treatment indicator (coded as 1 and 2).}
#'   \item{cov1, cov2, cov3, cov4, cov5, cov6}{Pre-treatment covariates.}
#'   \item{sampling_prob}{Sampling probability computed as a function of covariates and treatment (e.g., with an additional term depending on \code{trt}).}
#'   \item{survey_weight}{Survey weight computed as the inverse of \code{sampling_prob}.}
#'   \item{sample_indicator}{Indicator variable (0/1) denoting whether an observation was selected into the sample.}
#' }
#'
#' @details In the retrospective design, the sampling probability is determined by both
#' covariates and treatment assignment. This dataset represents the full superpopulation
#' in which the sampling mechanism is retrospective.
#'
#' @keywords datasets
#'
#' @examples
#' data(psdata_bin_retrospective_sp)
#' head(psdata_bin_retrospective_sp)
"psdata_bin_retrospective_sp"