#' covariance matrix for the multivriate normal distributed variables
#' 
#' @description covariance matrix of the normal distribution under cluster randomized study type given a design and a type
#' @param K number of timepoints or measurments (design parameter)
#' @param J number of subjects
#' @param I number of clusters (design parameter)
#' @param sigma.1.q variance of the lowest level (error variance or within subject variance)
#' @param sigma.2.q secound level variance (e.g. within cluster and between subject variance),  by default NULL and then a cross-sectional type
#' @param sigma.3.q third level variance (e.g. between cluster variance)
#' @return V covariance matrix
#' @examples
#' K<-6  #measurement (or timepoints)
#' I<-10 #Cluster
#' J<-2 #number of subjects
#' 
#' sigma.1<-0.1
#' sigma.3<-0.9
#' CovMat.Design(K, J, I,sigma.1.q=sigma.1, sigma.3.q=sigma.3)
#' 
#' sigma.1<-0.1
#' sigma.2<-0.4
#' sigma.3<-0.9
#' CovMat.Design(K, J, I,sigma.1.q=sigma.1, sigma.2.q=sigma.2, sigma.3.q=sigma.3)
#' @export
CovMat.Design<-function(K, J, I, sigma.1.q, sigma.2.q=NULL, sigma.3.q){

    
    if(is.null(sigma.2.q)){ #type=="cross-sec"
      #Covarianmatrix of one ...  (all measurements)
      W.j<-sigma.1.q*diag(1,nrow=K,ncol=K) #KxK  
    } else{#type=="long"
      #Covarianmatrix of one subject (all rep. measurements)
      W.j<-sigma.2.q*matrix(1,nrow=K,ncol=K)+sigma.1.q*diag(1,nrow=K,ncol=K) #KxK
    }
    
    x<-lapply(1:J, function(X) W.j)
    #blockMatrixDiagonal(sigma,sigma,sigma)
    W<-blockMatrixDiagonal(x) #J*K x J*K
    #Covarianmatrix of one cluster (all Masurements zu each timepoint of all/each subjects within one cluster)
    V.i<-sigma.3.q*matrix(1,nrow=J*K,ncol=J*K)+W#J*K x J*K
    #gemeinsame Verteilung aller Cluster
    x<-lapply(1:I, function(X) V.i)
    #blockMatrixDiagonal(sigma,sigma,sigma)
    V<-blockMatrixDiagonal(x)   
  
  
  return(V)
  
}
