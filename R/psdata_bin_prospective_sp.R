#' Prospective Superpopulation Dataset for PSweight Analysis
#'
#' A simulated dataset representing the prospective superpopulation (i.e. the full dataset)
#' for propensity score weighting analysis under a prospective design. In this design, the
#' sampling probability is determined solely by covariates.
#'
#' @docType data
#'
#' @usage data(psdata_bin_prospective_sp)
#'
#' @format A data frame with 1500 rows and 10 variables:
#' \describe{
#'   \item{Y}{Outcome variable.}
#'   \item{trt}{Binary treatment indicator (coded as 1 and 2).}
#'   \item{cov1, cov2, cov3, cov4, cov5, cov6}{Pre-treatment covariates.}
#'   \item{sampling_prob}{Sampling probability computed solely from covariates (e.g., using a logistic function).}
#'   \item{survey_weight}{Survey weight computed as the inverse of \code{sampling_prob}.}
#'   \item{sample_indicator}{Indicator variable (0/1) denoting whether an observation was selected into the sample.}
#' }
#'
#' @details In the prospective design, sampling is conducted before treatment assignment,
#' and the probability of an observation being selected depends only on its covariate values.
#'
#' @keywords datasets
#'
#' @examples
#' data(psdata_bin_prospective_sp)
#' head(psdata_bin_prospective_sp)
"psdata_bin_prospective_sp"