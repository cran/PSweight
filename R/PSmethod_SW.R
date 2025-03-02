#' Fitting propensity score models with different methods in survey observational data
#'
#' The function \code{PSmethod_SW} is an internal function to estimate the propensity scores given a specified model through formula.
#' It is built into function \code{Sumstat_SW}, \code{PStrim_SW} and \code{PSweight_SW}.
#'
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the propensity score model to be fitted. Additional details of model specification
#' are given under "Details". This argument is optional if \code{ps.estimate} is not \code{NULL}.
#' @param method a character to specify the method for estimating propensity scores. \code{"glm"} is default, and \code{"gbm"} and \code{"SuperLearner"} are also allowed.
#' @param data an optional data frame containing the variables in the propensity score model.
#' @param survey.weight an numeric vector specifying survey weights for each observation.
#' @param ncate a numeric to specify the number of treatment groups present in the given data.
#' @param ps.control a list to specify additional options when \code{method} is set to \code{"gbm"} or \code{"SuperLearner"}.
#'
#' @details  A typical form for \code{ps.formula} is \code{treatment ~ terms} where \code{treatment} is the treatment
#' variable and \code{terms} is a series of terms which specifies a linear predictor. \code{ps.formula} by default specifies generalized
#' linear models given the default argument \code{method = "glm"}.  It fits the logistic regression when \code{ncate = 2},and multinomial
#' logistic regression when \code{ncate > 2}. The argument \code{method} allows user to choose
#' model other than glm to fit the propensity score models. We have included \code{gbm} and \code{SuperLearner} as two alternative machine learning methods.
#' Additional arguments of the machine learning estimators can be supplied through the \code{...} argument. Note that SuperLearner does not handle multiple groups and the current version of multinomial
#' logistic regression is not supported by gbm. We suggest user to use them with extra caution. Please refer to the user manual of the \code{gbm} and \code{SuperLearner} packages for all the
#' allowed arguments.
#'
#'
#' @return
#'
#' \describe{
#'
#' \item{\code{ e.h}}{a data frame of estimated propensity scores.}
#'
#' \item{\code{ ps.fitObjects}}{the fitted propensity model details}
#' 
#' \item{\code{ beta.h}}{estimated coefficient of the propensity model when \code{method = "glm"}.}
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
#' # Define the propensity score model.
#' ps.formula <- trt ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6
#' 
#' # Extract the survey weights from the retrospective data.
#' survey.weight <- psdata_bin_retrospective_fp$survey_weight
#' 
#' # Specify the number of treatment groups (for binary treatment, ncate = 2)
#' ncate <- 2
#' 
#' # Fit the propensity score model using PSmethod_SW.
#' psfit <- PSmethod_SW(ps.formula = ps.formula,
#'                      data = psdata_bin_retrospective_fp,
#'                      survey.weight = survey.weight,
#'                      ncate = ncate)
#' 
#' # Print the first 10 rows of the estimated propensity scores.
#' cat("Estimated propensity scores (first 10 observations):\n")
#' print(head(psfit$e.h, 10))
#' 
#' # For the 'glm' method, print the estimated coefficient vector.
#' cat("\nEstimated coefficients (beta.h):\n")
#' print(psfit$beta.h)
#' 
#' # Users can also inspect the full fitted model object.
#' cat("\nFitted propensity model object summary:\n")
#' print(summary(psfit$ps.fitObjects))
#' 
#' 
PSmethod_SW<-function(ps.formula=ps.formula, method="glm", data=data, survey.weight = survey.weight, ncate=ncate, ps.control=list()){

  zname<-all.vars(ps.formula)[1]

  facz<-as.factor(data[,zname])
  #creat a dictionary for the original and recoded values in Z
  dic<-levels(facz)

  ps.fitObjects <- NULL  # Initialize an empty list to store fit objects

  ############## logistic #############################################################
  beta.h<-NULL #only return coefficient when gbm
  if (method=="glm"){
    if (exists("distribution")){
      warning("distribution argument only necessary for gbm method; set to bernoulli in glm method with logit link function by default")
    }

    if(ncate==2){
      #change z to 0/1 if
      dataps<-data

      dataps[,zname]<- as.numeric(facz)-1
      dataps$survey.weight <- survey.weight
      
      fitglm <- glm(formula = ps.formula, data=dataps,weights = survey.weight,family = quasibinomial(link = "logit"))
      
      e.h <- fitglm$fitted.values
      e.h <- cbind(1-e.h,e.h)
    }else{
      fitglm <- nnet::multinom(formula = ps.formula, data=data, weights = survey.weight,maxit = 500, Hess = TRUE, trace = FALSE)
      e.h <- fitglm$fitted.values
    }
    beta.h<-as.numeric(t(coef(fitglm)))
    ps.fitObjects <- fitglm  # Store the glm fit object
  }

  ############## gbm ###################################################################
  if (method=="gbm"){
    if ("distribution" %in% names(ps.control)){
      if (!ps.control$distribution %in% c("bernoulli","adaboost","multinomial")){
        stop("only bernoulli, adaboost, or multinomial distributions in 'gbm' are supported in propensity score models of PSweight")
      }
    }
    if ("weights" %in% names(ps.control)) warning("weights argument set to the survey weight input when estimating the population level propensity score.")
    ps.control$weights <- if (!is.null(survey.weight)) survey.weight else NULL


    if (!("var.monotone" %in% names(ps.control))) ps.control$var.monotone=NULL
    if (!("n.trees" %in% names(ps.control))) ps.control$n.trees=100
    if (!("interaction.depth" %in% names(ps.control))) ps.control$interaction.depth=1
    if (!("n.minobsinnode" %in% names(ps.control))) ps.control$n.minobsinnode=10
    if (!("shrinkage" %in% names(ps.control))) ps.control$shrinkage=0.1
    if (!("bag.fraction" %in% names(ps.control))) ps.control$bag.fraction=0.5
    if (!("train.fraction" %in% names(ps.control))) ps.control$train.fraction=1
    if (!("cv.folds" %in% names(ps.control))) ps.control$cv.folds=0
    if (!("class.stratify.cv" %in% names(ps.control))) ps.control$class.stratify.cv=NULL
    if (!("n.cores" %in% names(ps.control))) ps.control$n.cores=NULL
    if ("verbose" %in% names(ps.control)) warning("verbose argument set to F for gbm in PSweight")
    ps.control$verbose <- FALSE


    if(ncate==2){
      #change z to 0/1
      dataps<-data
      dataps[,zname]<- as.numeric(facz)-1

      if ("distribution" %in% names(ps.control)) {
        if (!ps.control$distribution %in% c("adaboost","bernoulli")) {
          ps.control$distribution<-"bernoulli"
          warning("supplied unsupported distribution for binary outcome in gbm; reset to bernoulli")
        }
      }else{
        ps.control$distribution<-"bernoulli"
      }

      fitgbm <- gbm::gbm(formula = ps.formula, data=dataps, distribution=ps.control$distribution,weights = ps.control$weights, var.monotone=ps.control$var.monotone, n.trees=ps.control$n.trees,
                       interaction.depth = ps.control$interaction.depth, n.minobsinnode = ps.control$n.minobsinnode, shrinkage = ps.control$shrinkage, bag.fraction = ps.control$bag.fraction,
                       train.fraction = ps.control$train.fraction, cv.folds = ps.control$cv.folds, keep.data = T, verbose = F,
                       class.stratify.cv = ps.control$class.stratify.cv, n.cores = ps.control$n.cores)



      e.h<-exp(fitgbm$fit)/(1+exp(fitgbm$fit))
      e.h<-cbind(1-e.h,e.h)
    }else if (ncate>2){
      if ("distribution" %in% names(ps.control)) {
        if (!ps.control$distribution=="multinomial"){
          warning("distribution for multi-category outcome reset to multinomial")
        }
      }

      ps.control$distribution<-"multinomial"
      warning("current multinomial distribution is broken in gbm; fitted results are rescaled to have rowsums of 1")

      fitgbm <- gbm::gbm(formula = ps.formula, data=data,distribution=ps.control$distribution, weights = ps.control$weights, var.monotone=ps.control$var.monotone, n.trees=ps.control$n.trees,
                       interaction.depth = ps.control$interaction.depth, n.minobsinnode = ps.control$n.minobsinnode, shrinkage = ps.control$shrinkage, bag.fraction = ps.control$bag.fraction,
                       train.fraction = ps.control$train.fraction, cv.folds = ps.control$cv.folds, keep.data = T, verbose = F,
                       class.stratify.cv = ps.control$class.stratify.cv, n.cores = ps.control$n.cores)


      #standardize
      e.h<-predict(fitgbm, newdata = data, type = "response", n.trees=ps.control$n.trees)[,,1] # stop warning with 'n.trees'

    }
    ps.fitObjects <- fitgbm  # Store the gbm fit object
  }

############## super learner #############################################################
  if (method=="SuperLearner"){
    if ("distribution" %in% names(ps.control)){
      warning("distribution argument not supported by SuperLearner; only family argument is supported")
    }
    if ("obsWeights" %in% names(ps.control))  warning("obsWeights argument set to the survey weight input when estimating the population level propensity score.")
    ps.control$obsWeights <- if (!is.null(survey.weight)) survey.weight else NULL
  

    if ("newX" %in% names(ps.control)) warning("newX argument set to NULL for SuperLearner in PSweight; please use method argument")
    ps.control$newX <- NULL
    if ("id" %in% names(ps.control)) warning("id argument set to NULL for SuperLearner in PSweight")
    ps.control$id <- NULL
    if ("verbose" %in% names(ps.control)) warning("verbose argument set to F for SuperLearner in PSweight")
    ps.control$verbose <- FALSE
    if (!("control" %in% names(ps.control))) ps.control$control<-list()
    if (!("cvControl" %in% names(ps.control))) ps.control$cvControl<-list()
    if (!("env" %in% names(ps.control))) ps.control$env<-parent.frame()


    family="quasibinomial"
    zvalue<-as.numeric(facz)-1

    if(ncate>2){
      stop("only binary outcomes are supported in SuperLearner for propensity score model in PSweight")
    }else{
      SL.all<-c("SL.bartMachine","SL.bayesglm","SL.biglasso","SL.caret","SL.caret.rpart","SL.cforest",
                "SL.earth","SL.extraTrees", "SL.gam","SL.gbm","SL.glm","SL.glm.interaction","SL.glmnet","SL.ipredbagg",
                "SL.kernelKnn","SL.knn","SL.ksvm","SL.lda","SL.leekasso","SL.lm", "SL.loess","SL.logreg","SL.mean","SL.nnet",
                "SL.nnls","SL.polymars","SL.qda","SL.randomForest","SL.ranger","SL.ridge","SL.rpart","SL.rpartPrune","SL.speedglm",
                "SL.speedlm","SL.step","SL.step.forward","SL.step.interaction","SL.stepAIC","SL.svm","SL.template","SL.xgboost")

        if ("SL.library" %in% names(ps.control)){
          if (length(unlist(ps.control$SL.library))>1) {
            ps.control$SL.library<-unlist(ps.control$SL.library)[1]
            warning("only one method allowed in SL.library argument of SuperLearner in PSweight; the first element in SL.library is taken")
          }
          if (!ps.control$SL.library %in% SL.all) {
            stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values")
          }
          }else{ #no SL.library specified
        ps.control$SL.library="SL.glm"
      }

      covM<-model.matrix(ps.formula, data)
      
      #method is fixed to "method.NNLS"
      fitsl <-SuperLearner::SuperLearner(Y=zvalue, X=data.frame(covM) , newX = NULL, family = family, SL.library= ps.control$SL.library,
                                         method = "method.NNLS", id = NULL, verbose = FALSE, control =  ps.control$control, cvControl =  ps.control$cvControl,
                                         obsWeights =  ps.control$obsWeights, env = ps.control$env)



      e.h<-cbind((1-fitsl$SL.predict),fitsl$SL.predict)
    }
    ps.fitObjects <- fitsl  # Store the SuperLearner fit object

  }

  #relabel the propensity name
  colnames(e.h)<-dic

  return(list(e.h=e.h,beta.h=beta.h, ps.fitObjects = ps.fitObjects))
}

