#' Retrospective Factual Sample Dataset for PSweight Analysis
#'
#' A simulated dataset representing the observed (factual) sample drawn from the retrospective
#' superpopulation. This dataset is obtained by selecting only those observations with
#' \code{sample_indicator == 1}, where the sampling probability is a function of both covariates and treatment.
#'
#' @docType data
#'
#' @usage data(psdata_bin_retrospective_fp)
#'
#' @format A data frame with the subset of rows from \code{psdata_bin_retrospective_sp} (e.g. 600 rows)
#' and the same 10 variables as \code{psdata_bin_retrospective_sp}.
#'
#' @details This dataset is derived from \code{psdata_bin_retrospective_sp} by retaining only
#' observations that were selected (i.e. those with \code{sample_indicator == 1}). It represents
#' the observed sample used for analysis under a retrospective design.
#'
#' @keywords datasets
#'
#' @examples
#' data(psdata_bin_retrospective_fp)
#' head(psdata_bin_retrospective_fp)
"psdata_bin_retrospective_fp"