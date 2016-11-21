
#' #' Sampling Response of individuals within a SWD model
#'
#' @description Sample data (response) for given numbers of individuals by given a model (of a parallel, cross-sectional, stepped wedge design study) 
#' @param type of the design is either cross-sectional or longitudinal 
#' @param K number of timepoints or measurments (design parameter)
#' @param J number of subjects
#' @param I number of clusters (design parameter)
#' @param D a complete data design matrix corresponding to the assumed model
#' @param A a complete data design matrix corresponding to the true data, if A is null, then A is equal to D
#' @param V covariance matrix for the normal distribution 
#' @param parameters corresponding to the model (regression fixed effects coefficients)
#' @return Data of individuals intensities corresponds to the SWD model and full model parameter information
#' @examples
#' K<-6  #measurement (or timepoints)
#' I<-10 #Cluster
#' J<-2 #number of subjects
#' X<-designMatrix(nC=I, nT=K, nSw=2)
#' D<-completeDataDesignMatrix(J, X)
#' sigma.1<-0.1
#' sigma.3<-0.9
#' type<-"cross-sec"
#' V<-CovMat_Design(K, J, I, sigma.1=sigma.1, sigma.3=sigma.3)
#' mu.0<-0
#' theta<-1
#' betas<-rep(0, K-1)
#' parameters<-c(mu.0, betas, theta)
#' sample.data<-sampleData(type = type, K=K,J=J,I=I, D=D, V=V, parameters=parameters)
#' xtabs(~cluster+measurement, data=sample.data)
#' @export
sampleData<-function(type, K,J,I, D, A=NULL, V, parameters ){
  
  library("mvtnorm")
  if(is.null(A)){A<-D}
  #sample I cluster from distribution of cluster, each have the same mean vector
  mean.vec<-A%*%parameters
  sample<-rmvnorm(1, mean=mean.vec, sigma=V) #row = number of subjects. # col = number of rep.Meas
  
  if(type=="cross-sec"){
    
    sample.data<-data.frame(val=as.vector(sample),
                            measurement = as.factor(rep(rep(1:K, times=J), times=I)),
                            #subject=rep(rep(1:J, each=K), times=I),
                            subject=1:(I*J*K),
                            cluster=rep(1:I, each = J*K),
                            intervention=D[,dim(D)[2]]
    )
  }
  
  if(type=="long"){
    
    sample.data<-data.frame(val=as.vector(sample),
                            measurement = as.factor(rep(rep(1:K, times=J), times=I)),
                            #subject=rep(rep(1:J, each=K), times=I),
                            subject=rep(1:(I*J), each=K),
                            cluster=rep(1:I, each = J*K),
                            intervention=D[,dim(D)[2]]
    )
  }
  
  completeData<-sample.data
  return(completeData)
}
