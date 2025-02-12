#### Point and variance for binary treatment in survey observational data##############################################
###########################################################################################

binest_SW<-function(ps.formula=NULL,ps.estimate=NULL,zname=NULL,yname,data,trtgrp=NULL, 
                    survey.design = 'Retrospective', svywtname = NULL, 
                    augmentation=FALSE, augmentation.type = 'WET',
                    bootstrap=FALSE,R=50,out.formula=NULL,out.estimate=NULL,family=NULL,
                    weight="overlap",ps.method='glm',ps.control=list(),out.method='glm',out.control=list()){

  #preprocess formula and extract y
  out.formula<-as.formula(out.formula)
  y<-unlist(data[yname])
  
  #extract survey weight
  survey.weight <- NULL
  if(is.null(svywtname)){
    data$survey_weight <- 1
    survey.weight <- data$survey_weight
  }else{
    survey.weight <- as.numeric(unlist(data[svywtname]))
  }
  
  #set ordered group
  facz<-as.factor(unlist(data[zname]))

  #set the treatment label
  oldlevel<-levels(facz)
  newlevel<-oldlevel

  z<-as.numeric(facz)-1
  n<-length(z) #total obs

  #set group for ATT
  if(weight=="treated"){
    if(is.null(trtgrp)){
      trtgrp<-oldlevel[2]
    }else{
      newlevel<-rev(unique(c(trtgrp,oldlevel)))
      facz<-factor(unlist(data[zname]),levels = newlevel)
      z<-as.numeric(facz)-1
    }
  }

  matchlevel<-match(oldlevel,newlevel) #match level for ATT

  #reassign value as numeric
  data[zname]<-z
  data_p<-data[,colnames(data)!=yname] #data without outcome

  # obtain ps estimation
  # estimate with formula
  if(is.null(ps.estimate)){
    #estimate population level propensity score
    fit<-do.call(PSmethod_SW,list(ps.formula = ps.formula, method=ps.method, data=data_p,survey.weight = survey.weight,ncate=2,ps.control=ps.control))
    W<- model.matrix(ps.formula,data_p)                      # design matrix

    e.h <- as.numeric(fit$e.h[,2])
    beta.h<-as.numeric(fit$beta.h)
    
    #estimate sample level propensity score
    fit.sample<-do.call(PSmethod,list(ps.formula = ps.formula, method=ps.method, data=data_p,ncate=2,ps.control=ps.control))
    e.h.sample <- as.numeric(fit.sample$e.h[,2])
    beta.h.sample<-as.numeric(fit.sample$beta.h)
    
    #construct ratio for transformation between different types of survey weights
    r.z <- ifelse(data_p$trt==1, e.h/e.h.sample, (1-e.h)/(1-e.h.sample))
    
    if(ps.method!='glm'){
      ps.estimate<-fit$e.h
      ps.estimate.sample<-fit.sample$e.h
    }
    
  }else{
    stop("External propensity score estimates not supported for population level estimation.")
  }

  if(weight=="entropy"){
    stop("Entropy weight not supported in survey setting.")
  }
 
  #tilting function
  ftilt<-tiltbin(weight = weight)

  if (survey.design == "Retrospective") { 
    data$balancing.wt <- ifelse(z==1, ftilt(e.h)/e.h, ftilt(e.h)/(1-e.h)) * survey.weight
    data$h.p <- ftilt(e.h) * survey.weight * (1/r.z)
  } else if (survey.design == "Independent") { 
    data$balancing.wt <- ifelse(z==1, ftilt(e.h.sample)/e.h.sample, ftilt(e.h.sample)/(1-e.h.sample)) * survey.weight
    data$h.p <- ftilt(e.h.sample) * survey.weight
  } else if (survey.design == "Prospective") { 
    data$balancing.wt <- ifelse(z==1, ftilt(e.h)/e.h.sample, ftilt(e.h)/(1-e.h.sample)) * survey.weight
    data$h.p <- ftilt(e.h) * survey.weight
  } else {
    stop("survey.design must be 'Retrospective', 'Independent' or 'Prospective'.")
  }

  #compute outcome regression for augmentation
  if(augmentation){
    if(is.null(out.estimate) && length(out.formula) == 0){
      stop("When augmentation = TRUE and out.estimate is not provided, a valid out.formula must be supplied.")
    }
    #no outcome estimation provided
    if(is.null(out.estimate)){
      #fit two outcome regression model for different treatment groups
      offset.e<-rep(1,n) #for poisson regression
      dataaug<-data[,colnames(data)!=zname]
      dataaug0<-dataaug[z==0,]
      dataaug1<-dataaug[z==1,]
      out.weights0 <- as.numeric(dataaug0$balancing.wt)
      out.weights1 <- as.numeric(dataaug1$balancing.wt)
      

      if(augmentation.type == "MOM"){
        
        XY<-model.matrix(formula(out.formula),data=dataaug)
        
        #predict outcome
        fitout0<-do.call(OUTmethod_SW,list(out.formula=out.formula,y=y[z==0], out.method=out.method, family=family, datain=dataaug0, dataout=dataaug,out.control=out.control))
        m0.h<-fitout0$m.est
        gamma0.h<-fitout0$gamma.h
        fitout1<-do.call(OUTmethod_SW,list(out.formula=out.formula,y=y[z==1], out.method=out.method, family=family, datain=dataaug1, dataout=dataaug,out.control=out.control))
        m1.h<-fitout1$m.est
        gamma1.h<-fitout1$gamma.h

      }else if (augmentation.type == "CVR"){
        
        out.formula.cvr <- update(out.formula, . ~ . + balancing.wt)
        
        XY<-model.matrix(formula(out.formula.cvr),data=dataaug)
        
        #predict outcome
        fitout0<-do.call(OUTmethod_SW,list(out.formula=out.formula.cvr,y=y[z==0], out.method=out.method, family=family, datain=dataaug0, dataout=dataaug,out.control=out.control))
        m0.h<-fitout0$m.est
        gamma0.h<-fitout0$gamma.h
        fitout1<-do.call(OUTmethod_SW,list(out.formula=out.formula.cvr,y=y[z==1], out.method=out.method, family=family, datain=dataaug1, dataout=dataaug,out.control=out.control))
        m1.h<-fitout1$m.est
        gamma1.h<-fitout1$gamma.h
        

      }else if (augmentation.type == "WET"){
        
        XY<-model.matrix(formula(out.formula),data=dataaug)

        #predict outcome
        fitout0<-do.call(OUTmethod_SW,list(out.formula=out.formula, out.weights = out.weights0, y=y[z==0], out.method=out.method, family=family, datain=dataaug0, dataout=dataaug, out.control=out.control))
        m0.h<-fitout0$m.est
        gamma0.h<-fitout0$gamma.h
        fitout1<-do.call(OUTmethod_SW,list(out.formula=out.formula, out.weights = out.weights1, y=y[z==1], out.method=out.method, family=family, datain=dataaug1, dataout=dataaug, out.control=out.control))
        m1.h<-fitout1$m.est
        gamma1.h<-fitout1$gamma.h
   
      }else{
        stop("augmentation.type must be 'MOM', 'CVR', or 'WET'.")
      }
      
      

      if(family=='poisson'){
        offsetlog<-model.extract(model.frame(out.formula,data = dataaug),'offset')
        if(!is.null(offsetlog)){offset.e<-exp(offsetlog)}
      }
      if(out.method!='glm'){
        out.estimate<-cbind(m0.h,m1.h)
        colnames(out.estimate)<-newlevel
      }

    }else{
      #the name for the outcome regression
      if(!setequal(colnames(out.estimate),newlevel)){
        out.estimate<-as.matrix(out.estimate)
        m0.h<-out.estimate[,1]
        m1.h<-out.estimate[,2]
        warning("wrong column name set for out.estimate, treatment set as: ",newlevel[1], " , ", newlevel[2])
      }else{
        out.estimate<-out.estimate[,match(newlevel,colnames(out.estimate))]
        m0.h<-out.estimate[,1]
        m1.h<-out.estimate[,2]
      }
    }
  }


  ##No Augmentation###############################################################
  if(augmentation==FALSE){
    muhat<-ptbin_SW(e.h,e.h.sample,z,y,ftilt,survey.design,survey.weight,r.z)
    ##No bootstrap###############################################################
    if(bootstrap==FALSE){

      conser<-1   #choose conservative or not

      #use ps formula
      if(is.null(ps.estimate)){
        tryCatch( {
          
          if (survey.design == "Retrospective") { 
            theta.h<-c(muhat,beta.h)
          } else if (survey.design == "Independent") { 
            theta.h<-c(muhat,beta.h.sample)
          } else if (survey.design == "Prospective") { 
            theta.h<-c(muhat,beta.h,beta.h.sample)
          } else {
            stop("survey.design must be 'Retrospective', 'Independent' or 'Prospective'.")
          }

          covmu<-sand_bin_SW(z,y,n,ftilt,survey.weight,r.z,theta.h,W,survey.design = survey.design,type='e')
          conser<-0 #is pd
        },error = function(w) {
          warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
        })
      }

      #use conservative
      if(conser==1){
        theta.h<-muhat
        covmu<-sand_bin_SW(z,y,n,ftilt,survey.weight,r.z,theta.h, eest=e.h, eest.sample =e.h.sample, survey.design = survey.design,type='ec')
      }

      names(muhat)<-newlevel
      colnames(covmu)<-rownames(covmu)<-newlevel
      muboot<-NULL

      
    }else{
      warning("Bootstrap not supported by PSweight in the survey setting.")
    }
  }

  ##Augmentation###############################################################
  if(augmentation==TRUE){

    #calculate point estimate
    augest <- ptbin_SW(e.h,e.h.sample,z,y,ftilt,survey.design,survey.weight,r.z,m0.h,m1.h)
    muhat <- tail(augest,2)

    ##No bootstrap###############################################################
    if(bootstrap==FALSE){
      
      if(augmentation.type %in% c('MOM','CVR')){
        out.weights <- rep(1,n)
      }else if(augmentation.type == 'WET'){
        out.weights <- data$balancing.wt
      }else{
        stop("augmentation.type must be 'MOM', 'CVR', or 'WET'.")
      }

      conser<-1       #choose conservative or not

      if(is.null(ps.estimate) & is.null(out.estimate)){

        ## both with formula
        tryCatch({

          if (survey.design == "Retrospective") { 
            theta.h<-c(augest[1:6],beta.h,beta.h.sample,gamma0.h,gamma1.h)
          } else if (survey.design == "Independent") { 
            theta.h<-c(augest[1:6],beta.h.sample,gamma0.h,gamma1.h)
          } else if (survey.design == "Prospective") { 
            theta.h<-c(augest[1:6],beta.h,beta.h.sample,gamma0.h,gamma1.h)
          } else {
            stop("survey.design must be 'Retrospective', 'Independent' or 'Prospective'.")
          }
          
          covmu<-sand_bin_SW(z,y,n,ftilt,survey.weight,r.z,theta.h,W,XY,survey.design=survey.design,out.weights=out.weights,family=family,offset.e=offset.e,type='ea')
          conser<-0 #is pd
        },error = function(w) {
          warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
        })

      }else if(!is.null(ps.estimate) & is.null(out.estimate)){
      

        ## formula on outcome regression only
        tryCatch( {
          
          if (survey.design == "Retrospective") { 
            theta.h<-c(augest[1:6],gamma0.h,gamma1.h)
          } else if (survey.design == "Independent") { 
            theta.h<-c(augest[1:6],gamma0.h,gamma1.h)
          } else if (survey.design == "Prospective") { 
            theta.h<-c(augest[1:6],gamma0.h,gamma1.h)
          } else {
            stop("survey.design must be 'Retrospective', 'Independent' or 'Prospective'.")
          }
          
          covmu<-sand_bin_SW(z,y,n,ftilt,survey.weight,r.z,theta.h,XY=XY,eest=e.h,eest.sample = e.h.sample,survey.design = survey.design, out.weights = out.weights, family=family,offset.e=offset.e,type='eca')
          conser<-0 #is pd
        },error = function(w) {
          warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
        })

      }else if(is.null(ps.estimate) & !is.null(out.estimate)){

        ## formula on ps only
        #when not pd
        tryCatch( {
          if (survey.design == "Retrospective") { 
            theta.h<-c(augest[1:6],beta.h,beta.h.sample)
          } else if (survey.design == "Independent") { 
            theta.h<-c(augest[1:6],beta.h.sample)
          } else if (survey.design == "Prospective") { 
            theta.h<-c(augest[1:6],beta.h,beta.h.sample)
          } else {
            stop("survey.design must be 'Retrospective', 'Independent' or 'Prospective'.")
          }

          covmu<-sand_bin_SW(z,y,n,ftilt,survey.weight,r.z,theta.h,W,survey.design = survey.design, m0est=m0.h,m1est=m1.h,type='eac')
          conser<-0 #is pd
        },error = function(w) {
          warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
        })
      }

      
      #if not pd then use conservative
      if(conser==1){
        theta.h<-augest[1:6]
        covmu<-sand_bin_SW(z,y,n,ftilt,survey.weight,r.z,theta.h,eest=e.h,eest.sample = e.h.sample,survey.design = survey.design,m0est=m0.h,m1est=m1.h,type='ecac')
      }

      muboot<-NULL
      names(muhat)<-newlevel
      colnames(covmu)<-rownames(covmu)<-newlevel
    }else{
      warning("Bootstrap not supported by PSweight in the survey setting.")
    }
  }
  e.h<-cbind(1-e.h,e.h)
  colnames(e.h)<-newlevel
  e.h.sample<-cbind(1-e.h.sample,e.h.sample)
  colnames(e.h.sample)<-newlevel

  #match back to original levels
  e.h<-e.h[,matchlevel]
  e.h.sample<-e.h.sample[,matchlevel]
  muhat<-muhat[matchlevel]
  covmu<-covmu[matchlevel,matchlevel]
  muboot<-muboot[,matchlevel]

  out<-list(propensity=e.h, propensity.sample=e.h.sample, muhat=muhat, covmu=covmu, muboot=muboot, group=c(oldlevel), trtgrp=trtgrp)
  class(out)<-'PSweight'
  out
}




