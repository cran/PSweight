# Outcome regression method

outmethod<-function(out.formula=out.formula, y=y, methodout="glm", family='gaussian', datain=datain, dataout=dataout,...){


################################################################################################
  gamma.h<-NULL #only useful in glm
  if(methodout=='glm'){

    if(family=='gaussian'){
      fitglm<-lm(out.formula,data=datain)
    }else if(family=='binomial'){
      fitglm<-glm(out.formula,family = binomial(link = "logit"),data=datain)
    }else if(family=='poisson'){
      fitglm<-glm(out.formula,family = poisson(),data=datain)
    }
    m.est<-predict(fitglm,type = "response",dataout)
    gamma.h<-as.numeric(coef(fitglm))

  }else if(methodout=='gbm'){
    if (exists("distribution")){
      if (!distribution %in% c("bernoulli","adaboost","gaussian","poisson")){
        stop("only bernoulli, adaboost, gaussian, poisson, supported for augmentation")
      }
    }

    if (exists("var.monotone")) var.monotone=var.monotone else var.monotone=NULL
    if (exists("weights")) weights=weights else weights=NULL
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

    if(family=='gaussian'){
      distribution<-"gaussian"
    }else if(family=='binomial'){
      if (!exists("distribution")){
        distribution<-"bernoulli"
      }
    }else if(family=='poisson'){
      distribution<-"poisson"
    }

    fitgbm <- gbm::gbm(formula = out.formula, data=datain,distribution=distribution, var.monotone=var.monotone, n.trees=n.trees,
                     interaction.depth = interaction.depth, n.minobsinnode = n.minobsinnode, shrinkage = shrinkage, bag.fraction = bag.fraction,
                     train.fraction = train.fraction, cv.folds = cv.folds, keep.data = T, verbose = F,
                     class.stratify.cv = class.stratify.cv, n.cores = n.cores)


    m.est<-predict(fitgbm,type = "response",newdata=dataout)

  }else if(methodout=='SuperLearner'){


    if (exists("newX")) warning("newX argument set to NULL for SuperLearner in PSweight; please use methodout argument")
    if (exists("id")) warning("id argument set to NULL for SuperLearner in PSweight")
    if (exists("verbose")) warning("verbose argument set to F for SuperLearner in PSweight")
    if (exists("obsWeights")) obsWeights=obsWeights else obsWeights=NULL
    if (exists("control")) control=control else control=list()
    if (exists("cvControl")) cvControl=cvControl else cvControl=list()
    if (exists("env")) env=env else env=parent.frame()

    if(family=='poisson'){
      warning("poisson regression not supported in SuperLearner")
      family<-'gaussian'
    }

    SL.all<-c("SL.bartMachine","SL.bayesglm","SL.biglasso","SL.caret","SL.caret.rpart","SL.cforest",
              "SL.earth","SL.extraTrees", "SL.gam","SL.gbm","SL.glm","SL.glm.interaction","SL.glmnet","SL.ipredbagg",
              "SL.kernelKnn","SL.knn","SL.ksvm","SL.lda","SL.leekasso","SL.lm", "SL.loess","SL.logreg","SL.mean","SL.nnet",
              "SL.nnls","SL.polymars","SL.qda","SL.randomForest","SL.ranger","SL.ridge","SL.rpart","SL.rpartPrune","SL.speedglm",
              "SL.speedlm","SL.step","SL.step.forward","SL.step.interaction","SL.stepAIC","SL.svm","SL.template","SL.xgboost")


      if (exists("SL.library")){
        if (length(unlist(SL.library))>1) {
          SL.library=unlist(SL.library)[1]
          warning("only one methodout allowed in SL.library argument of SuperLearner in PSweight; the first element in SL.library is taken")
        }else {
          if (!SL.library %in% SL.all){
            stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values")
          }
        }
      }else{ #no SL.library specified
        if(family=='binomial'){
          SL.library="SL.glm"
        }else{
          SL.library="SL.lm"
        }
      }

    covM<-model.matrix(out.formula, datain)
    #methodout is fixed to "method.NNLS"
    fitSL <-SuperLearner::SuperLearner(Y=y, X=data.frame(covM), newX = NULL, family = family, SL.library=SL.library,
                                     method = "method.NNLS", id = NULL, verbose = FALSE, control = control, cvControl = cvControl,
                                     obsWeights = obsWeights, env = env)

    covout<-model.matrix(out.formula,dataout)

    m.est <- predict(fitSL, newdata=data.frame(covout))$pred
  }



  return(list(m.est=m.est,gamma.h=gamma.h))
}





















