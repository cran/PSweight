#fit propensity with different models

PSmethod<-function(ps.formula=ps.formula, method="glm", data=data, ncate=ncate,...){

  zname<-all.vars(ps.formula)[1]

  facz<-as.factor(data[,zname])
  #creat a dictionary for the original and recoded values in Z
  dic<-levels(facz)

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

      fitglm <- glm(formula = ps.formula, data=dataps,family = binomial(link = "logit"))
      e.h <- fitglm$fitted.values
      e.h <- cbind(1-e.h,e.h)
    }else{
      fitglm <- nnet::multinom(formula = ps.formula, data=data, maxit = 500, Hess = TRUE, trace = FALSE)
      e.h <- fitglm$fitted.values
    }
    beta.h<-as.numeric(t(coef(fitglm)))
  }

  ############## gbm ###################################################################
  if (method=="gbm"){
    if (exists("distribution")){
      if (!distribution %in% c("bernoulli","adaboost","multinomial")){
        stop("only bernoulli, adaboost, or multinomial distributions in 'gbm' are supported in propensity score models of PSweight")
      }
    }

    if (exists("var.monotone")) var.monotone=var.monotone else var.monotone=NULL
    if (exists("weights")) weights=weights else weights=NULL
    if (exists("n.trees")) n.trees=n.trees else n.trees=100
    if (exists("interaction.depth")) interaction.depth=interaction.depth else interaction.depth=1
    if (exists("n.minobsinnode")) n.minobsinnode=n.minobsinnode else n.minobsinnode=10
    if (exists("shrinkage")) shrinkage=shrinkage else shrinkage=0.1
    if (exists("bag.fraction")) bag.fraction=bag.fraction else bag.fraction=0.5
    if (exists("train.fraction")) train.fraction=train.fraction else train.fraction=1
    if (exists("cv.folds")) cv.folds=cv.folds else cv.folds=0
    if (exists("class.stratify.cv ")) class.stratify.cv =class.stratify.cv  else class.stratify.cv=NULL
    if (exists("n.cores ")) n.cores =n.cores  else n.cores=NULL
    if (exists("verbose")) warning("verbose argument set to F for SuperLearner in PSweight")

    if(ncate==2){
      #change z to 0/1
      dataps<-data
      dataps[,zname]<- as.numeric(facz)-1

      if (exists("distribution")) {
        if (!distribution %in% c("adaboost","bernoulli")) {
          distribution<-"bernoulli"
          warning("supplied unsupported distribution for binary outcome in gbm; reset to bernoulli")
        }
      }else{
        distribution<-"bernoulli"
      }

      fitgbm <- gbm::gbm(formula = ps.formula, data=dataps,distribution=distribution, var.monotone=var.monotone, n.trees=n.trees,
                       interaction.depth = interaction.depth, n.minobsinnode = n.minobsinnode, shrinkage = shrinkage, bag.fraction = bag.fraction,
                       train.fraction = train.fraction, cv.folds = cv.folds, keep.data = T, verbose = F,
                       class.stratify.cv = class.stratify.cv, n.cores = n.cores)

      e.h<-exp(fitgbm$fit)/(1+exp(fitgbm$fit))
      e.h<-cbind(1-e.h,e.h)
    }else if (ncate>2){
      if (exists("distribution")) {
        if (!distribution=="multinomial"){
          warning("distribution for multi-category outcome reset to multinomial")
        }
      }

      distribution<-"multinomial"

      warning("current multinomial distribution is broken in gbm; fitted results are rescaled to have rowsums of 1")

      fitgbm <- gbm::gbm(formula = ps.formula, data=data,distribution=distribution, var.monotone=var.monotone, n.trees=n.trees,
                       interaction.depth = interaction.depth, n.minobsinnode = n.minobsinnode, shrinkage = shrinkage, bag.fraction = bag.fraction,
                       train.fraction = train.fraction, cv.folds = cv.folds, keep.data = T, verbose = F,
                       class.stratify.cv = class.stratify.cv, n.cores = n.cores)
      #standardize
      e.hpre<-apply(fitgbm$fit,2,function(x){exp(x)/(1+exp(x))})
      e.h<-apply(e.hpre, 2, function(x) x/rowSums(e.hpre))
    }
  }

############## super learner #############################################################
  if (method=="SuperLearner"){
    if (exists("distribution")){
      warning("distribution argument not supported by SuperLearner; only family argument is supported")
    }

    if (exists("newX")) warning("newX argument set to NULL for SuperLearner in PSweight; please use method argument")
    if (exists("id")) warning("id argument set to NULL for SuperLearner in PSweight")
    if (exists("verbose")) warning("verbose argument set to F for SuperLearner in PSweight")
    if (exists("obsWeights")) obsWeights<-obsWeights else obsWeights=NULL
    if (exists("control")) control<-control else control<-list()
    if (exists("cvControl")) cvControl<-cvControl else cvControl<-list()
    if (exists("env")) env<-env else env<-parent.frame()


    family="binomial"
    zvalue<-as.numeric(facz)-1

    if(ncate>2){
      stop("only binary outcomes are supported in SuperLearner for propensity score model in PSweight")
    }else{
      SL.all<-c("SL.bartMachine","SL.bayesglm","SL.biglasso","SL.caret","SL.caret.rpart","SL.cforest",
                "SL.earth","SL.extraTrees", "SL.gam","SL.gbm","SL.glm","SL.glm.interaction","SL.glmnet","SL.ipredbagg",
                "SL.kernelKnn","SL.knn","SL.ksvm","SL.lda","SL.leekasso","SL.lm", "SL.loess","SL.logreg","SL.mean","SL.nnet",
                "SL.nnls","SL.polymars","SL.qda","SL.randomForest","SL.ranger","SL.ridge","SL.rpart","SL.rpartPrune","SL.speedglm",
                "SL.speedlm","SL.step","SL.step.forward","SL.step.interaction","SL.stepAIC","SL.svm","SL.template","SL.xgboost")

      if (exists("SL.library")){
        if (length(unlist(SL.library))>1) {
          SL.library<-unlist(SL.library)[1]
          warning("only one method allowed in SL.library argument of SuperLearner in PSweight; the first element in SL.library is taken")
        }
        if (!SL.library %in% SL.all) stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values")
      }else{ #no SL.library specified
        SL.library="SL.glm"
      }

      covM<-model.matrix(ps.formula, data)
      #method is fixed to "method.NNLS"
      fitsl <-SuperLearner::SuperLearner(Y=zvalue, X=data.frame(covM), newX = NULL, family = family, SL.library=SL.library,
                                         method = "method.NNLS", id = NULL, verbose = FALSE, control = control, cvControl = cvControl,
                                         obsWeights = obsWeights, env = env)

      e.h<-cbind((1-fitsl$SL.predict),fitsl$SL.predict)
    }

  }

  #relabel the propensity name
  colnames(e.h)<-dic

  return(list(e.h=e.h,beta.h=beta.h))
}





















