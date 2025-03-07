% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_SumStat_SW.R
\name{plot.SumStat_SW}
\alias{plot.SumStat_SW}
\title{Plot the distribution of propensity scores and balance statistics}
\usage{
\method{plot}{SumStat_SW}(
  x,
  type = "balance",
  weighted.var = TRUE,
  threshold = 0.1,
  metric = "ASD",
  breaks = 50,
  ...
)
}
\arguments{
\item{x}{a \code{SumStat} object obtained with \code{\link{SumStat}} function.}

\item{type}{a character indicating the type of plot to produce, including histogram of estimated propensity scores (\code{"hist"}), density of estimated propensity scores (\code{"density"}), and plot of balance statistics (\code{"balance"}).}

\item{weighted.var}{logical. Indicating whether weighted variance should be used in calculating the balance statistics. Default is \code{TRUE}.}

\item{threshold}{an optional numeric value indicating the balance threshold for the balance plot. Default is 0.1. Only valid when \code{type = "balance"}.}

\item{metric}{a character indicating the type of metric used in balance plot. Only \code{"ASD"} or \code{"PSD"} is allowed. If not specified, the default is \code{"ASD"}. See  \code{\link{summary.SumStat}} for additional details on balance metrics.}

\item{breaks}{a single number giving the number of cells for the histogram. Default is 50.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
Plot of the indicated type.
}
\description{
Summarize the SumStat object, generate histogram or density of estimated propensity scores and
plot the balance statistics under weighting versus no weighting.
}
\details{
For the balance plot, a vertical line at \code{threshold} is used to define balance on covariates.
The default value is \code{threshold = 0.1} following Austin and Stuart (2015). If more than 2 treatments
are considered, only density of the estimated generalized propensity scores will be produced, regardless of
whether \code{type = "density"} or \code{type = "hist"}.
}
\examples{
data("psdata")
data("psdata_bin_prospective_fp")
data("psdata_bin_retrospective_fp")

ps.formula<-trt~cov1+cov2+cov3+cov4+cov5+cov6
msstat <- SumStat_SW(ps.formula = ps.formula, trtgrp = "2", 
                     data = psdata_bin_prospective_fp,
                     weight=c("IPW","overlap","treated","entropy","matching"), 
                     survey.indicator = TRUE, svywtname = "survey_weight", 
                     survey.design = "Prospective")

plot(msstat, type="hist")
plot(msstat, type="balance", weighted.var=TRUE, threshold=0.1, metric="ASD")

}
\references{
Austin, P.C. and Stuart, E.A. (2015). Moving towards best practice when using inverse probability of treatment weighting (IPTW) using the propensity score to estimate causal treatment effects in observational studies.
Statistics in Medicine, 34(28), 3661-3679.
}
