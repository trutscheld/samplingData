#' Design matrix for SWD model 
#'
#' @description create design matrix for a given setup of a stepped wedge design
#' @param nC number of cluster
#' @param nT number of timepoints
#' @param nSw number of cluster : within parallel recieve the control (nC-nSw receive the intervention), within cross-over recieve the pattern (0, 1) (nC-nSw receive the pattern (1,0)) for nearly the same number of time points, within SWD switches from control to intervention per time point 
#' @param swP is the time point the cluster cross over the condition in a cross over study, if not given then it is nearly half of the time past
#' @param design is the study type (parallel, cross-sectional, stepped wedge)
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
designMatrix<-function(nC, nT, nSw, swP=NULL, design="SWD"){
  
  
  if(design=="parallel"){
    
    ma<-rbind(t(replicate(nSw, rep(0, nT))), t(replicate(nC-nSw, rep(1, nT))))
  }
  
  if(design=="cross-over"){
    
    if(is.null(swP)){swP<-ceiling(nT/2)}
    swP<-
    ma<-rbind(t(replicate(nSw, c(rep(0, swP), rep(1,nT-swP)))), 
              t(replicate(nC-nSw, c(rep(1,swP), rep(0, nT-swP)))))
  }
  
  
  
  if(design=="SWD"){
    ma<-sapply(1:nT, function(i){
      
      noTr<-(i-1)*nSw
      noC<-nC-noTr
      return(c(rep(1,noTr), rep(0,noC)))
    })
  }
  
  
  
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
#' K<-6  #measurement (or timepoints)
#' I<-10 #Cluster
#' J<-2 #number of subjects
#' X<-designMatrix(nC=I, nT=K, nSw=2)
#' completeDataDesignMatrix(J, X)
#' @export
completeDataDesignMatrix<-function(J, X){
  
  if(!is.matrix(X)){
    stop("X must be a matrix")
  }
  if(is.na(J)|is.null(J)|(J==0)){
    stop("The number of subjects must be at least 1")
  }
  
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
