#' Fitting potential outcome regression with different methods in survey observational data
#'
#' The function \code{OUTmethod_SW} is an internal function to estimate the potential outcomes given a specified model through formula.
#' It is built into function \code{PSweight}, and is used for constructing the augmented estimators.
#'
#' @param out.formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the outcome model to be fitted.
#' @param out.weights A numeric vector of weights to be used in the weighted regression estimator. 
#' For the moment estimator(MOM) and clever covariates estimator(CVR), this parameter defaults to NULL; 
#' however, in the weighted regression estimator(WET) it should be set to the balancing weights used for the weighted outcome regression model.
#' @param y a vector of the observed outcome in the training data (\code{datain}).
#' @param out.method a character to specify the method for estimating the outcome regression model. \code{"glm"} is default, and \code{"gbm"} and \code{"SuperLearner"} are also allowed.
#' @param family a description of the error distribution and link function to be used in the outcome model. Supported distributional families include
#' \code{"gaussian" (link = identity)}, \code{"binomial" (link = logit)} and \code{"poisson" (link = log)}. Default is \code{"gaussian"}.
#' @param datain The training data for the outcome model. In the context of \code{PSweight}, it refers to the data observed for each treatment group.
#' @param dataout The prediction data for the outcome model. In the context of \code{PSweight}, it refers to the full data.
#' @param out.control a list to specify additional options when \code{out.method} is set to \code{"gbm"} or \code{"SuperLearner"}.
#'
#' @details  A typical form for \code{out.formula} is \code{y ~ terms} where \code{y} is the outcome
#' variable and \code{terms} is a series of terms which specifies linear predictors (on the link function scale). \code{out.formula} by default specifies generalized
#' linear model given the gaussian error through the default arguments \code{method = "glm"} and \code{family = "gaussian"}.  It fits the logistic regression when \code{family = "binomal"},and poisson
#' regression when \code{family = "poisson"}. The argument \code{out.method} allows user to choose
#' model other than glm to fit the outcome regression models for constructing the augmented estimator. We have included \code{gbm} and \code{SuperLearner} as alternative machine learning estimators.
#' Additional argument in them can be supplied through the \code{...} argument. Please refer to the user manual of the \code{gbm} and \code{SuperLearner} packages for all the
#' allowed arguments.
#' 
#' When \code{out.weights} is provided, weighted regression is applied in all supported methods. In \code{"glm"}, 
#' weighted least squares is used for Gaussian outcomes, and weighted likelihood for binomial or Poisson outcomes. 
#' If \code{family = "poisson"} and non-integer outcomes are detected (due to continuous-valued responses or weighting effects), 
#' the function automatically switches to \code{"gaussian"} to ensure numerical stability. In \code{"gbm"}, outcome 
#' weights are passed via the \code{weights} argument, and if non-integer outcomes exist under \code{family = "poisson"}, 
#' the model switches to \code{"gaussian"}. For \code{"SuperLearner"}, weighted estimation is supported through the 
#' \code{obsWeights} argument. Since \code{"SuperLearner"} does not support Poisson regression, a warning is issued, 
#' and the model defaults to \code{"gaussian"}. These adjustments ensure compatibility and stability across different 
#' estimation approaches.
#'
#' @return
#'
#' \describe{
#'
#' \item{\code{ m.est}}{a vector of predicted outcome on the \code{dataout}.}
#'
#' \item{\code{ gamma.h}}{estimated coefficient of the outcome model when \code{method = "glm"}.}
#'
#' }
#'
#' @export
#'
#' @examples
#' # Load example datasets
#' data("psdata")
#' data("psdata_bin_prospective_fp")
#' data("psdata_bin_retrospective_fp")
#' 
#' # Define the outcome model formula.
#' out.formula <- Y ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6
#' 
#' # Extract the outcome vector from the retrospective data.
#' y <- psdata_bin_retrospective_fp$Y
#' 
#' # Use only the observations in treatment group 1 as the training data.
#' datain <- psdata_bin_retrospective_fp[psdata_bin_retrospective_fp$trt == 1, ]
#' 
#' # Fit the outcome regression model via OUTmethod_SW.
#' # By default, out.method = "glm" and family = "gaussian" are used.
#' outfit <- OUTmethod_SW(out.formula = out.formula, y = y, datain = datain,
#'                        dataout = psdata_bin_retrospective_fp)
#' 
#' # Print the predicted outcome vector on dataout.
#' cat("Predicted outcomes (first 10 values):\n")
#' print(head(outfit$m.est, 10))
#' 
#' # Print the estimated coefficient vector from the GLM.
#' cat("\nEstimated coefficients (gamma.h):\n")
#' print(outfit$gamma.h)
#'
#' @importFrom stats quasibinomial update
#'
OUTmethod_SW <- function(out.formula = out.formula, out.weights = NULL, y = y,
                         out.method = "glm", family = "gaussian",
                         datain = datain, dataout = dataout, out.control = list()) {
  
  gamma.h <- NULL  # Only used in "glm" method
  datain$out.weights <- out.weights
  
  if (out.method == "glm" && family == "poisson") {
    if (!is.null(out.weights) || any(y != round(y))) {
      warning("GLM with Poisson family: non-integer outcomes detected (due to weights and/or non-integer y values); defaulting to gaussian")
      family <- "gaussian"
    }
  }
  if (out.method == "gbm" && family == "poisson") {
    if (!is.null(out.weights) || any(y != round(y))) {
      warning("GBM with Poisson distribution: non-integer outcomes detected (due to weights and/or non-integer y values); defaulting to gaussian")
      out.control$distribution <- "gaussian"
      family <- "gaussian"
    }
  }
  
  if (out.method == "glm") {
    # Fit weighted regression if out.weights is provided; otherwise, fit unweighted
    if (!is.null(out.weights)) {
      if (family %in% c("gaussian")) {
        fitglm <- lm(out.formula, data = datain, weights = out.weights)
      } else if (family == "binomial") {
        fitglm <- glm(out.formula, family = binomial(link = "logit"), data = datain, weights = out.weights)
      }
    } else {
      if (family == "gaussian") {
        fitglm <- lm(out.formula, data = datain)
      } else if (family == "binomial") {
        fitglm <- glm(out.formula, family = binomial(link = "logit"), data = datain)
      } else if (family == "poisson") {
        fitglm <- glm(out.formula, family = poisson(), data = datain)
      }
    }
    m.est <- predict(fitglm, newdata = dataout, type = "response")
    gamma.h <- as.numeric(coef(fitglm))
    
  } else if (out.method == "gbm") {
    # Set default control parameters and incorporate out.weights if provided
    if (!("weights" %in% names(out.control))) out.control$weights <- out.weights
    if (!("var.monotone" %in% names(out.control))) out.control$var.monotone <- NULL
    if (!("n.trees" %in% names(out.control))) out.control$n.trees <- 100
    if (!("interaction.depth" %in% names(out.control))) out.control$interaction.depth <- 1
    if (!("n.minobsinnode" %in% names(out.control))) out.control$n.minobsinnode <- 10
    if (!("shrinkage" %in% names(out.control))) out.control$shrinkage <- 0.1
    if (!("bag.fraction" %in% names(out.control))) out.control$bag.fraction <- 0.5
    if (!("train.fraction" %in% names(out.control))) out.control$train.fraction <- 1
    if (!("cv.folds" %in% names(out.control))) out.control$cv.folds <- 0
    if (!("class.stratify.cv" %in% names(out.control))) out.control$class.stratify.cv <- NULL
    if (!("n.cores" %in% names(out.control))) out.control$n.cores <- NULL
    out.control$verbose <- FALSE
    
    if (family == "gaussian") {
      out.control$distribution <- "gaussian"
    } else if (family == "binomial") {
      out.control$distribution <- "bernoulli"
    } else if (family == "poisson") {
      out.control$distribution <- "poisson"
    }
    fitgbm <- gbm::gbm(formula = out.formula, data = datain,
                       distribution = out.control$distribution,
                       var.monotone = out.control$var.monotone,
                       weights = out.control$weights,  # Incorporate outcome weights if provided
                       n.trees = out.control$n.trees,
                       interaction.depth = out.control$interaction.depth,
                       n.minobsinnode = out.control$n.minobsinnode,
                       shrinkage = out.control$shrinkage,
                       bag.fraction = out.control$bag.fraction,
                       train.fraction = out.control$train.fraction,
                       cv.folds = out.control$cv.folds,
                       keep.data = TRUE,
                       verbose = FALSE,
                       class.stratify.cv = out.control$class.stratify.cv,
                       n.cores = out.control$n.cores)
    m.est <- predict(fitgbm, newdata = dataout, type = "response", n.trees = out.control$n.trees)
    
  } else if (out.method == "SuperLearner") {
    # If out.weights is provided, set obsWeights in out.control accordingly.
    if (!is.null(out.weights)) {
      out.control$obsWeights <- out.weights
    }
    out.control$newX <- NULL
    out.control$id <- NULL
    out.control$verbose <- FALSE
    if (!("control" %in% names(out.control))) out.control$control <- list()
    if (!("cvControl" %in% names(out.control))) out.control$cvControl <- list()
    if (!("env" %in% names(out.control))) out.control$env <- parent.frame()
    if (family == "poisson") {
      warning("poisson regression not supported in SuperLearner; defaulting to gaussian")
      family <- "gaussian"
    }
    SL.all <- c("SL.bartMachine", "SL.bayesglm", "SL.biglasso", "SL.caret", "SL.caret.rpart",
                "SL.cforest", "SL.earth", "SL.extraTrees", "SL.gam", "SL.gbm", "SL.glm",
                "SL.glm.interaction", "SL.glmnet", "SL.ipredbagg", "SL.kernelKnn", "SL.knn",
                "SL.ksvm", "SL.lda", "SL.leekasso", "SL.lm", "SL.loess", "SL.logreg", "SL.mean",
                "SL.nnet", "SL.nnls", "SL.polymars", "SL.qda", "SL.randomForest", "SL.ranger",
                "SL.ridge", "SL.rpart", "SL.rpartPrune", "SL.speedglm", "SL.speedlm", "SL.step",
                "SL.step.forward", "SL.step.interaction", "SL.stepAIC", "SL.svm", "SL.template",
                "SL.xgboost")
    if ("SL.library" %in% names(out.control)) {
      if (length(unlist(out.control$SL.library)) > 1) {
        out.control$SL.library <- unlist(out.control$SL.library)[1]
        warning("only one method allowed in SL.library; the first element is used")
      }
      if (!out.control$SL.library %in% SL.all) {
        stop("SL.library unrecognized; please use listWrappers() from SuperLearner to see supported methods")
      }
    } else {
      if (family == "binomial") {
        out.control$SL.library <- "SL.glm"
      } else {
        out.control$SL.library <- "SL.lm"
      }
    }
    covM <- model.matrix(out.formula, datain)
    fitSL <- SuperLearner::SuperLearner(Y = y, X = data.frame(covM), newX = NULL,
                                        family = family, SL.library = out.control$SL.library,
                                        method = "method.NNLS", id = NULL, verbose = FALSE,
                                        control = out.control$control, cvControl = out.control$cvControl,
                                        obsWeights = out.control$obsWeights, env = out.control$env)
    covout <- model.matrix(out.formula, dataout)
    m.est <- predict(fitSL, newdata = data.frame(covout))$pred
  }
  
  m.est <- as.numeric(m.est)
  return(list(m.est = m.est, gamma.h = gamma.h))
}
