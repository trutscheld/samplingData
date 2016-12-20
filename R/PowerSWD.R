#' Power calculation within stepped wedge design model by Hussey et.al or Heo&Kim 
#'
#' @description Calculation of power for a lmm with cluster as random effect, fixed timepoint effects, but set to null, TP number of timepoints, I number of cluster. The design matrix has to be coded by zeros and ones.
#' @param ThetaEst expected treatment effect
#' @param alpha singificance level (by default 0.05)
#' @param Design design matrix for a given SWD model
#' @param tauq between cluster variance
#' @param sigmaq within cluster variance(between subject variance)
#' @param sigmaq.error within subject variance/error variance
#' @param noSub number of subjects within each cluster and each timepoint (for all an equal size)
#' @param type is of "cross-sectional" (by default)  or "longitudinal" assigns the type of data (2 or 3 level nested structure)
#' @param time a logical (FALSE, if no time trends are expected, otherwise TRUE) is only relevant for evaluation of cross-sectional data
#' @return Aproximated power of two tailed test, although the design matrix is fractionated, then power is not valid
#' formula used for cross-sectional data provided by Michael A. Hussey and James P. Hughes,
#' Design and analysis of stepped wedge cluster randomized trials, 
#' Contemporary Clinical Trials(28),2007, 
#' and for longitudinal data by Heo M., Kim N., Rinke ML., Wylie-Rosett J., 
#' Sample size determinations for stepped-wedge clinical trials from a three-level data hierarchy perspective, 
#' Stat Methods Med Res., 2016
#' @examples
#' noCl<-10
#' noT<-6
#' switches<-2
#' DM<-designMatrix(noCl,noT,switches)
#' sigma.e <- 2
#' sigma.alpha <- 2   
#' #Power for cross-sectional SWD design by formula of Hussey&Hughes 
#' calcPower.SWD(ThetaEst=1,Design=DM, sigmaq=sigma.e^2, tauq=sigma.alpha^2, time=FALSE)
#' calcPower.SWD(ThetaEst=1,Design=DM, sigmaq=sigma.e^2, tauq=sigma.alpha^2, time=TRUE)
#' #Power for longitudinal SWD design by formula of Heo&Kim 
#' DM.new<-NULL
#' for(i in 1:dim(DM)[2]){
#' DM.new<-cbind(DM.new,DM[,i], DM[,i])
#' }
#' sigma.e <- sqrt(7/10)
#' sigma <- sqrt(2/10)
#' sigma.alpha <- sqrt(1/10 )
#' K<- 10 #number of participants within each 'cell'
#' calcPower.SWD(ThetaEst=1, Design=DM.new, tauq=sigma.alpha^2, sigmaq=sigma^2, sigmaq.error =sigma.e^2,  noSub=K, type="longitudinal")
#' @export

calcPower.SWD<-function(ThetaEst,alpha=0.05, Design, sigmaq,  tauq, sigmaq.error =NULL,  noSub=NULL, time=TRUE, type="cross-sectional"){
  
  #quantile of the standard normal distribution function
  Z<-qnorm(p=1-alpha/2)
  
  if(type=="cross-sectional"){#power for cross-sectional, formula by Hussey&Hughes
    
    #number of cluster
    I<-dim(Design)[1]
    #number of timepoints
    TP<-dim(Design)[2]
    #number of all clusters and timepoints, which get the treatment
    U<-sum(Design)
    #W:sum over all timepoints: (sum of all clusters with treatment)^2
    W<-sum(colSums(Design)^2)
    #V:sum over all clusters: (sum of all timepoints with treatment)^2
    V<-sum(rowSums(Design)^2)
    
    if(is.null(sigmaq)){  #if not given, than can be calculated by  sigmaq.error/noSub
      
      if(is.null(sigmaq.error)||is.null(noSub)) {
        stop("A value of sigmaq.error and/or noSub is NA. Both is required if sigmaq not given for cross-sectional data.")
      }else{
        sigmaq<-sigmaq.error/noSub
      }     
    }
    
    if(time==TRUE){#with time trend
      #Variance of estimated Theta by lmm estimated with WLS and Designmatrixed encoded with nulls and zeros
      varTheta<-(I*sigmaq*(sigmaq+TP*tauq))/((I*U-W)*sigmaq+(U^2+I*TP*U-TP*W-I*V)*tauq) 
    }
    else{#no time trends
      #Variance of estimated TPheta by lmm estimated with WLS and Designmatrixed encoded with nulls and zeros
      varTheta<-(I*sigmaq*(sigmaq+TP*tauq))/((I*TP*U-U^2)*sigmaq+I*TP*(TP*U-V)*tauq) 
      
    }
    
  }else if (type=="longitudinal") {#power for longitudinal, formula by Heo&Kim, als wenn time=TRUE by H&H
    
    #number of participants
    K<-noSub  
    #number of periods/measurments under experimental condition in ith cluster
    J_i_E<-rowSums(Design)
    #number of periods/measurments under control condition in ith cluster
    J_i_C<-rowSums(abs(Design-1))
    sumJs<-sum(sapply(1:dim(Design)[1],function(i){J_i_E[i]*J_i_C[i]}))
    J<-dim(Design)[2]
    sigmaq.1<-sigmaq.error
    sigmaq.2<-sigmaq
    sigmaq.3<-tauq
    sigmaq_ges<-(sigmaq.1+sigmaq.2+sigmaq.3)
    f.multiplied.sigma<-(sigmaq.1+sigmaq.2+sigmaq.3 + (K-1)*(sigmaq.2+sigmaq.3) - K*sigmaq.3)
    varTheta<-f.multiplied.sigma*J/(K*sumJs)
    #rho.1<-(sigmaq.2+sigmaq.3)/sigmaq_ges
    #sigma^2/(sigma.e^2+sigma^2+sigma.alpha^2)
    #rho.2<-sigmaq.3/sigmaq_ges
    #f<-(1 + (K-1)*rho.1 - K*rho.2)
    #varTheta<-(J*f*sigmaq_ges)/(K*sumJs)
  }
  
  #calculated power by Z distributed value, 
  #for both methods the same, only different Variances of treatment effect
  power<-pnorm(q=(ThetaEst/sqrt(varTheta)-Z))
  return(power)
  
  
}
