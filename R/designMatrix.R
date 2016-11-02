#' Design matrix for SWD model 
#'
#' @description create design matrix for a given setup of a stepped wedge design
#' @param nC number of cluster
#' @param nT number of timepoints
#' @param nSw number of cluster switches per time point from control to intervention
#' @return design matrix for a given setup of a stepped wedge design
#' @examples
#' 
#' designMatrix.SWD(5,6,1)
#'
#' K<-6  #measurement (or timepoints)
#' I<-10 #Cluster
#' designMatrix.SWD(nC=I, nT=K, nSw=2)
#'
#' @export
designMatrix.SWD<-function(nC, nT, nSw){
  
  ma<-sapply(1:nT, function(i){
    
    noTr<-(i-1)*nSw
    noC<-nC-noTr
    return(c(rep(1,noTr), rep(0,noC)))
  })
  return(ma)
  
}


#' Design matrix for SWD model under a grade of intervention implementation pattern
#'
#' @description Creates a implementation matrix for a given stepped wedge design and grade of intervention implementation pattern
#' @param nC Number of clusters
#' @param nT Number of timepoint
#' @param nSw number of clusters switches from control to treatment at each timepoint
#' @param pattern a vector for grade of intervention implementation pattern, which gives the derivation from 100 percent effectiveness over time
#' @return Design matrix for SWD model under a grade of intervention implementation pattern
#' @examples
#' 
#' implemMatrix.SWD(5,6,1, c(seq(0.4,1,0.2),1))
#' 
#' K<-6  #measurement (or timepoints)
#' I<-10 #Cluster
#' implemMatrix.SWD(nC=I, nT=K, nSw=2, pattern=c(seq(0.4,1,0.2),1))
#' @export
implemMatrix.SWD<-function(nC, nT, nSw, pattern){
  
  if(length(pattern)+1==nT){
    ma<-sapply(1:nT, function(i){
      
      noTr<-(i-1)*nSw
      noC<-nC-noTr      
      vec<-c(if( i>1 ){rep(pattern[(i-1):1],each=noTr/(i-1))} , 
             rep(0,noC) )
      return(vec)
    })
    return(ma) 
  }else{
    stop("T he length of the pattern must be one less than the number of timepoints.")
  }
  
}


#' Design matrix for complete data within design
#'
#' @description create design matrix for complete data within design
#' @param J number of subjects
#' @param X  given design matrix
#' @return design matrix for complete data within design
#' @examples
#' 
#' designMatrix.SWD(5,6,1)
#'
#' K<-6  #measurement (or timepoints)
#' I<-10 #Cluster
#' J<-2 number of subjects
#' X<-designMatrix.SWD(nC=I, nT=K, nSw=2)
#' completeDataDesignMatrix(J, X)
#' @export
completeDataDesignMatrix<-function(J, X){
  
  I<-nrow(X)
  K<-ncol(X)
  D<-NULL
  for(i in 1:I){
    
    D.i.j<-cbind(rep(1, K), #erste Spalte f?r mu
                 rbind(rep(0,K-1), diag(1,nrow=K-1,ncol=K-1)),#for times      
                 X[i,]#ergibt Spalte f?r Theta  in gro?er Effektmarix
    )
    D.i<-do.call("rbind", rep(list(D.i.j), J))
    D<-rbind(D,D.i )
  }
  
  return(D)
  
}
