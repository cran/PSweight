## Notation #########################
# single e stands for no augmentation
# ea stands for augmentation
# c represents conservative variance
#####################################


################## Sandwich Variance in binary group case in survey observational data####################################################
sand_bin_SW<-function(z,y,n,ftilt,svywt,r.z,thetaest,W=NULL,XY=NULL,
                      eest,eest.sample,survey.design,
                      out.weights,
                      m0est,m1est,family="gaussian",offset.e,
                      type="e"){
  

  p<-ncol(W)
  q<-ncol(XY)
  
  
  #matrix to calucalte mu0,mu1 correlation
  a<-matrix(0,2,length(thetaest))
  
  if (survey.design == "Retrospective") { 
    
    ##score function
    if(type=="e"){
      phi<-function(theta){
        mu0<-theta[1]
        mu1<-theta[2]
        beta<-theta[3:(p+2)]
        e<-plogis(c(W%*%beta))
        
        f1<-(y-mu0)*(1-z)*svywt*ftilt(e)/(1-e)
        f2<-(y-mu1)*z*svywt*ftilt(e)/e
        f3<-svywt*W*(z-e)
        
        f<-rbind(f1,f2,t(f3))
        return(f)
      }
      a[1,1]<-1
      a[2,2]<-1
    }else if(type=="ec"){
      phi<-function(theta){
        mu0<-theta[1]
        mu1<-theta[2]
        
        f1<-(y-mu0)*(1-z)*svywt*ftilt(eest)/(1-eest)
        f2<-(y-mu1)*z*svywt*ftilt(eest)/eest
        
        f<-rbind(f1,f2)
        return(f)
      }
      a[1,1]<-1
      a[2,2]<-1
    }else if(type=="ea"){
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        beta<-theta[7:(p+6)]
        beta.sample <- theta[(p+7):(2*p+6)]
        gamma0<-theta[(2*p+7):(2*p+6+q)]
        gamma1<-theta[(2*p+7+q):(2*p+2*q+6)]
        e<-plogis(c(W%*%beta))
        e.sample<-plogis(c(W%*%beta.sample))
        r.z<-ifelse(z == 1, e/e.sample, (1 - e)/(1 - e.sample))
        if (family=='gaussian'){
          m0<-c(XY%*%gamma0)
          m1<-c(XY%*%gamma1)
        }else if(family=='binomial'){
          m0<-plogis(c(XY%*%gamma0))
          m1<-plogis(c(XY%*%gamma1))
        }else if(family=='poisson'){
          m0<-exp(c(XY%*%gamma0))*offset.e
          m1<-exp(c(XY%*%gamma1))*offset.e
        }
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(e)/(1-e)
        f2<-(1/r.z)*svywt*ftilt(e)*(m0-aug0h)
        f3<-(1-z)*(m0-aug0z)*svywt*ftilt(e)/(1-e)
        
        f4<-z*(y-mu1)*svywt*ftilt(e)/e
        f5<-(1/r.z)*svywt*ftilt(e)*(m1-aug1h)
        f6<-z*(m1-aug1z)*svywt*ftilt(e)/e
        
        f7<-out.weights*XY*(y-m0)*(1-z)
        f8<-out.weights*XY*(y-m1)*z
        f9<-svywt*W*(z-e)
        f10<-W*(z-e.sample)
        
        f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8),t(f9),t(f10))
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }else if(type=='ecac'){
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(eest)/(1-eest)
        f2<-(1/r.z)*svywt*ftilt(eest)*(m0est-aug0h)
        f3<-(1-z)*(m0est-aug0z)*svywt*ftilt(eest)/(1-eest)
        
        f4<-z*(y-mu1)*svywt*ftilt(eest)/eest
        f5<-(1/r.z)*svywt*ftilt(eest)*(m1est-aug1h)
        f6<-z*(m1est-aug1z)*svywt*ftilt(eest)/eest
        
        f<-rbind(f1,f2,f3,f4,f5,f6)
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }else if(type=='eca'){
      #define the m estimator
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        gamma0<-theta[7:(6+q)]
        gamma1<-theta[(7+q):(2*q+6)]
        if (family=='gaussian'){
          m1<-c(XY%*%gamma1)
          m0<-c(XY%*%gamma0)
        }else if(family=='binomial'){
          m1<-plogis(c(XY%*%gamma1))
          m0<-plogis(c(XY%*%gamma0))
        }else if(family=='poisson'){
          m1<-exp(c(XY%*%gamma1))*offset.e
          m0<-exp(c(XY%*%gamma0))*offset.e
        }
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(eest)/(1-eest)
        f2<-(1/r.z)*svywt*ftilt(eest)*(m0-aug0h)
        f3<-(1-z)*(m0-aug0z)*svywt*ftilt(eest)/(1-eest)
        
        f4<-z*(y-mu1)*svywt*ftilt(eest)/eest
        f5<-(1/r.z)*svywt*ftilt(eest)*(m1-aug1h)
        f6<-z*(m1-aug1z)*svywt*ftilt(eest)/eest
        
        f7<-out.weights*XY*(y-m0)*(1-z)
        f8<-out.weights*XY*(y-m1)*z
        
        f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8))
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }else if(type=="eac"){
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        beta<-theta[7:(p+6)]
        beta.sample<-theta[(p+7):(2*p+6)]
        e<-plogis(c(W%*%beta))
        e.sample<-plogis(c(W%*%beta.sample))
        r.z<-ifelse(z == 1, e/e.sample, (1 - e)/(1 - e.sample))
        
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(e)/(1-e)
        f2<-(1/r.z)*svywt*ftilt(e)*(m0est-aug0h)
        f3<-(1-z)*(m0est-aug0z)*svywt*ftilt(e)/(1-e)
        
        f4<-z*(y-mu1)*svywt*ftilt(e)/e
        f5<-(1/r.z)*svywt*ftilt(e)*(m1est-aug1h)
        f6<-z*(m1est-aug1z)*svywt*ftilt(e)/e
        
        f7<-svywt*W*(z-e)
        f8<-W*(z-e.sample)
        
        f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8))
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }

    
  } else if (survey.design == "Independent") { 
    
    
    ##score function
    if(type=="e"){
      phi<-function(theta){
        mu0<-theta[1]
        mu1<-theta[2]
        beta.sample<-theta[3:(p+2)]
        e.sample<-plogis(c(W%*%beta.sample))
        
        f1<-(y-mu0)*(1-z)*svywt*ftilt(e.sample)/(1-e.sample)
        f2<-(y-mu1)*z*svywt*ftilt(e.sample)/e.sample
        f3<-W*(z-e.sample)
        
        f<-rbind(f1,f2,t(f3))
        return(f)
      }
      a[1,1]<-1
      a[2,2]<-1
    }else if(type=="ec"){
      phi<-function(theta){
        mu0<-theta[1]
        mu1<-theta[2]
        
        f1<-(y-mu0)*(1-z)*svywt*ftilt(eest.sample)/(1-eest.sample)
        f2<-(y-mu1)*z*svywt*ftilt(eest.sample)/eest.sample
        
        f<-rbind(f1,f2)
        return(f)
      }
      a[1,1]<-1
      a[2,2]<-1
    }else if(type=="ea"){
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        beta.sample<-theta[7:(p+6)]
        gamma0<-theta[(p+7):(p+6+q)]
        gamma1<-theta[(p+7+q):(p+2*q+6)]
        e.sample<-plogis(c(W%*%beta.sample))
        
        if (family=='gaussian'){
          m0<-c(XY%*%gamma0)
          m1<-c(XY%*%gamma1)
        }else if(family=='binomial'){
          m0<-plogis(c(XY%*%gamma0))
          m1<-plogis(c(XY%*%gamma1))
        }else if(family=='poisson'){
          m0<-exp(c(XY%*%gamma0))*offset.e
          m1<-exp(c(XY%*%gamma1))*offset.e
        }
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(e.sample)/(1-e.sample)
        f2<-svywt*ftilt(e.sample)*(m0-aug0h)
        f3<-(1-z)*(m0-aug0z)*svywt*ftilt(e.sample)/(1-e.sample)
        
        f4<-z*(y-mu1)*svywt*ftilt(e.sample)/e.sample
        f5<-svywt*ftilt(e.sample)*(m1-aug1h)
        f6<-z*(m1-aug1z)*svywt*ftilt(e.sample)/e.sample
        
        f7<-out.weights*XY*(y-m0)*(1-z)
        f8<-out.weights*XY*(y-m1)*z
        f9<-W*(z-e.sample)
        
        f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8),t(f9))
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }else if(type=='ecac'){
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(eest.sample)/(1-eest.sample)
        f2<-svywt*ftilt(eest.sample)*(m0est-aug0h)
        f3<-(1-z)*(m0est-aug0z)*svywt*ftilt(eest.sample)/(1-eest.sample)
        
        f4<-z*(y-mu1)*svywt*ftilt(eest.sample)/eest.sample
        f5<-svywt*ftilt(eest.sample)*(m1est-aug1h)
        f6<-z*(m1est-aug1z)*svywt*ftilt(eest.sample)/eest.sample
        
        f<-rbind(f1,f2,f3,f4,f5,f6)
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }else if(type=='eca'){
      #define the m estimator
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        gamma0<-theta[7:(6+q)]
        gamma1<-theta[(7+q):(2*q+6)]
        if (family=='gaussian'){
          m1<-c(XY%*%gamma1)
          m0<-c(XY%*%gamma0)
        }else if(family=='binomial'){
          m1<-plogis(c(XY%*%gamma1))
          m0<-plogis(c(XY%*%gamma0))
        }else if(family=='poisson'){
          m1<-exp(c(XY%*%gamma1))*offset.e
          m0<-exp(c(XY%*%gamma0))*offset.e
        }
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(eest.sample)/(1-eest.sample)
        f2<-svywt*ftilt(eest.sample)*(m0-aug0h)
        f3<-(1-z)*(m0-aug0z)*svywt*ftilt(eest.sample)/(1-eest.sample)
        
        f4<-z*(y-mu1)*svywt*ftilt(eest.sample)/eest.sample
        f5<-svywt*ftilt(eest.sample)*(m1-aug1h)
        f6<-z*(m1-aug1z)*svywt*ftilt(eest.sample)/eest.sample
        
        f7<-out.weights*XY*(y-m0)*(1-z)
        f8<-out.weights*XY*(y-m1)*z
        
        f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8))
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }else if(type=="eac"){
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        beta.sample<-theta[7:(p+6)]
        e.sample<-plogis(c(W%*%beta.sample))
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(e.sample)/(1-e.sample)
        f2<-svywt*ftilt(e.sample)*(m0est-aug0h)
        f3<-(1-z)*(m0est-aug0z)*svywt*ftilt(e.sample)/(1-e.sample)
        
        f4<-z*(y-mu1)*svywt*ftilt(e.sample)/e.sample
        f5<-svywt*ftilt(e.sample)*(m1est-aug1h)
        f6<-z*(m1est-aug1z)*svywt*ftilt(e.sample)/e.sample
        
        f7<-W*(z-e.sample)
        
        f<-rbind(f1,f2,f3,f4,f5,f6,t(f7))
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }
    
    
  } else if (survey.design == "Prospective") { 
    
    
    ##score function
    if(type=="e"){
      phi<-function(theta){
        mu0<-theta[1]
        mu1<-theta[2]
        beta<-theta[3:(p+2)]
        beta.sample<-theta[(p+3):(2*p+2)]
        e<-plogis(c(W%*%beta))
        e.sample<-plogis(c(W%*%beta.sample))
        
        f1<-(y-mu0)*(1-z)*svywt*ftilt(e)/(1-e.sample)
        f2<-(y-mu1)*z*svywt*ftilt(e)/e.sample
        f3<-svywt*W*(z-e)
        f4<-W*(z-e.sample)
        
        f<-rbind(f1,f2,t(f3),t(f4))
        return(f)
      }
      a[1,1]<-1
      a[2,2]<-1
    }else if(type=="ec"){
      phi<-function(theta){
        mu0<-theta[1]
        mu1<-theta[2]
        
        f1<-(y-mu0)*(1-z)*svywt*ftilt(eest)/(1-eest.sample)
        f2<-(y-mu1)*z*svywt*ftilt(eest)/eest.sample
        
        f<-rbind(f1,f2)
        return(f)
      }
      a[1,1]<-1
      a[2,2]<-1
    }else if(type=="ea"){
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        beta<-theta[7:(p+6)]
        beta.sample <- theta[(p+7):(2*p+6)]
        gamma0<-theta[(2*p+7):(2*p+6+q)]
        gamma1<-theta[(2*p+7+q):(2*p+2*q+6)]
        e<-plogis(c(W%*%beta))
        e.sample<-plogis(c(W%*%beta.sample))
        
        if (family=='gaussian'){
          m0<-c(XY%*%gamma0)
          m1<-c(XY%*%gamma1)
        }else if(family=='binomial'){
          m0<-plogis(c(XY%*%gamma0))
          m1<-plogis(c(XY%*%gamma1))
        }else if(family=='poisson'){
          m0<-exp(c(XY%*%gamma0))*offset.e
          m1<-exp(c(XY%*%gamma1))*offset.e
        }
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(e)/(1-e.sample)
        f2<-svywt*ftilt(e)*(m0-aug0h)
        f3<-(1-z)*(m0-aug0z)*svywt*ftilt(e)/(1-e.sample)
        
        f4<-z*(y-mu1)*svywt*ftilt(e)/e.sample
        f5<-svywt*ftilt(e)*(m1-aug1h)
        f6<-z*(m1-aug1z)*svywt*ftilt(e)/e.sample
        
        f7<-out.weights*XY*(y-m0)*(1-z)
        f8<-out.weights*XY*(y-m1)*z
        f9<-svywt*W*(z-e)
        f10<-W*(z-e.sample)
        
        f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8),t(f9),t(f10))
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }else if(type=='ecac'){
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(eest)/(1-eest.sample)
        f2<-svywt*ftilt(eest)*(m0est-aug0h)
        f3<-(1-z)*(m0est-aug0z)*svywt*ftilt(eest)/(1-eest.sample)
        
        f4<-z*(y-mu1)*svywt*ftilt(eest)/eest.sample
        f5<-svywt*ftilt(eest)*(m1est-aug1h)
        f6<-z*(m1est-aug1z)*svywt*ftilt(eest)/eest.sample
        
        f<-rbind(f1,f2,f3,f4,f5,f6)
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }else if(type=='eca'){
      #define the m estimator
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        gamma0<-theta[7:(6+q)]
        gamma1<-theta[(7+q):(2*q+6)]
        
        if (family=='gaussian'){
          m1<-c(XY%*%gamma1)
          m0<-c(XY%*%gamma0)
        }else if(family=='binomial'){
          m1<-plogis(c(XY%*%gamma1))
          m0<-plogis(c(XY%*%gamma0))
        }else if(family=='poisson'){
          m1<-exp(c(XY%*%gamma1))*offset.e
          m0<-exp(c(XY%*%gamma0))*offset.e
        }
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(eest)/(1-eest.sample)
        f2<-svywt*ftilt(eest)*(m0-aug0h)
        f3<-(1-z)*(m0-aug0z)*svywt*ftilt(eest)/(1-eest.sample)
        
        f4<-z*(y-mu1)*svywt*ftilt(eest)/eest.sample
        f5<-svywt*ftilt(eest)*(m1-aug1h)
        f6<-z*(m1-aug1z)*svywt*ftilt(eest)/eest.sample
        
        f7<-out.weights*XY*(y-m0)*(1-z)
        f8<-out.weights*XY*(y-m1)*z
        
        f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8))
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }else if(type=="eac"){
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        beta<-theta[7:(p+6)]
        beta.sample<-theta[(p+7):(2*p+6)]
        e<-plogis(c(W%*%beta))
        e.sample<-plogis(c(W%*%beta.sample))
        
        f1<-(1-z)*(y-mu0)*svywt*ftilt(e)/(1-e.sample)
        f2<-svywt*ftilt(e)*(m0est-aug0h)
        f3<-(1-z)*(m0est-aug0z)*svywt*ftilt(e)/(1-e.sample)
        
        f4<-z*(y-mu1)*svywt*ftilt(e)/e.sample
        f5<-svywt*ftilt(e)*(m1est-aug1h)
        f6<-z*(m1est-aug1z)*svywt*ftilt(e)/e.sample
        
        f7<-svywt*W*(z-e)
        f8<-W*(z-e.sample)
        
        f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8))
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
    }
    
   
  } else {
    stop("survey.design must be 'Retrospective', 'Independent' or 'Prospective'.")
  }
  


  mphi<-function(theta){
    rowMeans(phi(theta))
  }

  #define the meat B, covariance operator
  Omega<-function(theta){
    phis<-phi(theta)
    return(tcrossprod(phis)/n)
  }

  Atheta<-jacobian(mphi,thetaest)
  invAtheta <- tryCatch(
    solve(Atheta),
    error = function(w) {
      message("Atheta singular => using MASS::ginv")
      MASS::ginv(Atheta)
    }
  )

  #calculate the sandwich variance
  Vtmp<-invAtheta%*%Omega(thetaest)%*%t(invAtheta)/n

  covmu<-a%*%Vtmp%*%t(a)

  return(covmu)
}


