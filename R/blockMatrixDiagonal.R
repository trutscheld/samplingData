#' diagonal block matrix
#'
#' @description create a diagonal block matrix
#' @param ... a list of matrices
#' @return diagonal block matrix concatinated from this list of matrices
#' @examples
#' m1<-matrix(round(runif(4*4),1),nrow=4,ncol=4)
#' m2<-matrix(round(runif(4*4),1),nrow=4,ncol=4)
#' blockMatrixDiagonal(m1,m2,m2,m1)
#'
#' sigma.1<-0.1
#' sigma.2<-0.4
#' J<-10 #subjects
#' I<-3 #cluster
#' V.i<-sigma.2*matrix(1,nrow=J,ncol=J)+sigma.1*diag(1, nrow=J,ncol=J) #Covarianmatrix of one cluster
#' x<-lapply(1:I, function(X) V.i)
#' blockMatrixDiagonal(x) #Covarianmatrix of all cluster
#' @export
blockMatrixDiagonal<-function(...){  
  matrixList<-list(...)
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
  
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  finalMatrix
}

# 
# 
# 
# 