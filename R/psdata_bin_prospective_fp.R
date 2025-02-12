#' Prospective Factual Sample Dataset for PSweight Analysis
#'
#' A simulated dataset representing the observed (factual) sample drawn from the prospective
#' superpopulation. This dataset is obtained by selecting only those observations with
#' \code{sample_indicator == 1}, where sampling depends solely on covariates.
#'
#' @docType data
#'
#' @usage data(psdata_bin_prospective_fp)
#'
#' @format A data frame with the subset of rows from \code{psdata_bin_prospective_sp} (e.g. 600 rows)
#' and the same 10 variables as \code{psdata_bin_prospective_sp}.
#'
#' @details This dataset is derived from \code{psdata_bin_prospective_sp} by retaining only
#' observations that were selected (i.e. those with \code{sample_indicator == 1}). It represents
#' the observed sample used for analysis under a prospective design.
#'
#' @keywords datasets
#'
#' @examples
#' data(psdata_bin_prospective_fp)
#' head(psdata_bin_prospective_fp)
"psdata_bin_prospective_fp"