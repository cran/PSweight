# The point estimating function for binary group in survey observational data
ptbin_SW<-function(psest,psest.fp,z,y,ftilt,survey.design, svywt,r.z, m0est=NULL,m1est=NULL){
  
  if (survey.design == "Retrospective") { 
    balancing.wt <- ifelse(z==1, ftilt(psest)/psest, ftilt(psest)/(1-psest)) * svywt
    h.p <- ftilt(psest) * svywt * (1/r.z)
  } else if (survey.design == "Independent") { 
    balancing.wt <- ifelse(z==1, ftilt(psest.fp)/psest.fp, ftilt(psest.fp)/(1-psest.fp)) * svywt
    h.p <- ftilt(psest.fp) * svywt
  } else if (survey.design == "Prospective") { 
    balancing.wt <- ifelse(z==1, ftilt(psest)/psest.fp, ftilt(psest)/(1-psest.fp)) * svywt
    h.p <- ftilt(psest) * svywt
  } else {
    stop("survey.design must be 'Retrospective', 'Independent' or 'Prospective'.")
  }

  mu1est <- sum(z*y*balancing.wt) / sum(z*balancing.wt)
  mu0est <- sum((1-z)*y*balancing.wt) / sum((1-z)*balancing.wt)
  

  if(is.null(m0est)){
    return(c(mu0est,mu1est))
  }else{
    aug0hest<-sum(h.p*m0est)/sum(h.p)
    aug0zest<-sum((1-z)*m0est*balancing.wt)/sum((1-z)*balancing.wt)
    aug1hest<-sum(h.p*m1est)/sum(h.p)
    aug1zest<-sum(z*m1est*balancing.wt)/sum(z*balancing.wt)
    
    mu0estaug<-mu0est+aug0hest-aug0zest
    mu1estaug<-mu1est+aug1hest-aug1zest

    return(c(mu0est,aug0hest,aug0zest,mu1est,aug1hest,aug1zest,mu0estaug,mu1estaug))
  }
}

