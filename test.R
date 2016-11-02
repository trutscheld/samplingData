
source("blockMatrixDiagonal.R")

m1<-matrix(round(runif(4*4),1),nrow=4,ncol=4)
m2<-matrix(round(runif(4*4),1),nrow=4,ncol=4)
blockMatrixDiagonal(m1,m2,m2,m1)

sigma.1<-0.1
sigma.2<-0.4
J<-10 #subjects
I<-3 #cluster
V.i<-sigma.2*matrix(1,nrow=J,ncol=J)+sigma.1*diag(1, nrow=J,ncol=J) #Covarianmatrix of one cluster
x<-lapply(1:I, function(X) V.i)
blockMatrixDiagonal(x) #Covarianmatrix of all cluster


rm(list=ls())
source("designMatrix.R")

I<-10 #Cluster
K<-6  #measurement (or timepoints)
designMatrix.SWD(nC=I, nT=K, nSw=2)
implemMatrix.SWD(nC=I, nT=K, nSw=2, pattern=c(seq(0.4,1,0.2),1))

rm(list=ls())
source("CovarianceMatrix.R")

K<-4#rep. measurement (unterste Ebene)
J<-2#Subjects
I<-3#Cluster

sigma.1<-1
sigma.3<-9
CovMat_Design(K=K, J=J, I=I,sigma.1=sigma.1, sigma.3=sigma.3, design="SWD", type="cross-sec")

sigma.2<-4
CovMat_Design(K, J, I,sigma.1=sigma.1, sigma.2=sigma.2, sigma.3=sigma.3, design="SWD", type="long")


rm(list=ls())
source("designMatrix.R")
source("blockMatrixDiagonal.R")
source("CovarianceMatrix.R")
source("SamplingDesign.R")
library(lme4)

K<-6  #measurement (or timepoints)
I<-10 #Cluster
sw<-2
X<-designMatrix.SWD(nC=I, nT=K, nSw=sw) #Designmatrix of SWD

J<-2 #number of subjects
D<-completeDataDesignMatrix(J, X)

sigma.1<-0.1
sigma.3<-0.9
type<-"cross-sec"
V<-CovMat_Design(K=K, J=J, I=I, sigma.1=sigma.1, sigma.3=sigma.3, design="SWD", type=type)

mu.0<-0
theta<-1
betas<-rep(0, K-1)
parameters<-c(mu.0, betas, theta)

sample.data<-sampleData(type = type, K=K,J=J,I=I, D=D, V=V, parameters=parameters )

xtabs(~cluster+measurement, data=sample.data) #how many observation per cluster and timepoint
xtabs(val~cluster+measurement, data=sample.data) #means of cluster and timepoint
lmer(val~intervention+measurement + (1|cluster), data=sample.data)


sigma.2<-0.4
type<-"long"
V<-CovMat_Design(K, J, I,sigma.1=sigma.1, sigma.2=sigma.2, sigma.3=sigma.3, design="SWD", type=type)
sample.data<-sampleData(type = type, K=K,J=J,I=I, D=D, V=V, parameters=parameters )

xtabs(~cluster+measurement, data=sample.data) #how many observation per cluster and timepoint
xtabs(val~cluster+measurement, data=sample.data) #means of cluster and timepoint
lmer(val~intervention+measurement + (1|cluster), data=sample.data)
lmer(val~intervention+measurement + (1|cluster) +(1|subject), data=sample.data)

rm(list=ls())
source("designMatrix.R")
source("PowerSWD.R")

noCl<-10
noT<-6
switches<-2
DM<-designMatrix.SWD(noCl,noT,switches)
sigma.e <- 2
sigma.alpha <- 2    
calcPower.SWD(ThetaEst=1,Design=DM, sigmaq=sigma.e^2, tauq=sigma.alpha^2, time=FALSE)
calcPower.SWD(ThetaEst=1,Design=DM, sigmaq=sigma.e^2, tauq=sigma.alpha^2, time=TRUE)
