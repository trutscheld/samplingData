
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
designMatrix(nC=I, nT=K, nSw=2)
implemMatrix.SWD(nC=I, nT=K, nSw=2, pattern=c(seq(0.4,1,0.2),1))

rm(list=ls())
source("CovarianceMatrix.R")

K<-4#rep. measurement (unterste Ebene)
J<-2#Subjects
I<-3#Cluster

sigma.1<-1
sigma.3<-9
CovMat_Design(K=K, J=J, I=I,sigma.1=sigma.1, sigma.3=sigma.3)

sigma.2<-4
CovMat_Design(K, J, I,sigma.1=sigma.1, sigma.2=sigma.2, sigma.3=sigma.3)


rm(list=ls())
source("designMatrix.R")
source("blockMatrixDiagonal.R")
source("CovarianceMatrix.R")
source("SamplingDesign.R")
library(lme4)

K<-6  #measurement (or timepoints)
I<-10 #Cluster
sw<-2
X<-designMatrix(nC=I, nT=K, nSw=sw) #Designmatrix of SWD

J<-2 #number of subjects
D<-completeDataDesignMatrix(J, X)

sigma.1<-0.1
sigma.3<-0.9
#type "cross-sec"
V<-CovMat_Design(K=K, J=J, I=I, sigma.1.q=sigma.1, sigma.3.q=sigma.3)

mu.0<-0
theta<-1
betas<-rep(0, K-1)
parameters<-c(mu.0, betas, theta)

sample.data<-sampleData(type = "cross-sec", K=K,J=J,I=I, D=D, V=V, parameters=parameters )

xtabs(~cluster+measurement, data=sample.data) #how many observation per cluster and timepoint
xtabs(val~cluster+measurement, data=sample.data) #means of cluster and timepoint
lmer(val~intervention+measurement + (1|cluster), data=sample.data)


sigma.2<-0.4
#type 
V<-CovMat_Design(K, J, I,sigma.1.q=sigma.1, sigma.2.q=sigma.2, sigma.3.q=sigma.3)
sample.data<-sampleData(type = "long", K=K,J=J,I=I, D=D, V=V, parameters=parameters )

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
DM<-designMatrix(noCl,noT,switches)
sigma.e <- 2
sigma.alpha <- 2 
clustersize<-10
Theta<-1
calcPower.SWD(ThetaEst=Theta,Design=DM, sigmaq=sigma.e^2/clustersize, tauq=sigma.alpha^2, time=FALSE)
calcPower.SWD(ThetaEst=Theta,Design=DM, sigmaq=NULL, tauq=sigma.alpha^2, sigmaq.error =sigma.e^2, noSub=clustersize, time=FALSE)
calcPower.SWD(ThetaEst=Theta,Design=DM, sigmaq=sigma.e^2/clustersize, tauq=sigma.alpha^2)
calcPower.SWD(ThetaEst=Theta,Design=DM.new, sigmaq=sigma.e^2/clustersize, tauq=sigma.alpha^2)


#Example Heo&Kim Table 1
sigma.e <- sqrt(7/10)
sigma <- sqrt(2/10)
sigma.alpha <- sqrt(1/10 )

###Table 1, 1 row
delta<- 0.3# treatment effect
rho.1<- 0.3#correlation among level 1 data
(sigma^2+sigma.alpha^2)/(sigma.e^2+sigma^2+sigma.alpha^2)
rho.2<- 0.1#correlation among level 2 data
sigma.alpha^2/(sigma.e^2+sigma^2+sigma.alpha^2)
K<- 5#number of participant within each period (each cell9)
b<- 0#number of baseline measurements
##fixed parameters
c<-2# number of cluster in each step
p<-2#number of periods within step
S<-5 # number of steps
sigma.quad<-1#gesamtvarianz
(sigma.e^2+sigma^2+sigma.alpha^2)
#calculate
(J<-b+p*S) # total number of periods per cluster
(f<-(1+(K-1)*rho.1-K*rho.2))
(sumJs<-c*p*S*(S+1)*(p*(S-1)+3*b)/6)
(Z<-qnorm(p=1-0.05/2))
(Var.delta<-f*J*sigma.quad/(K*sumJs))
pnorm(q=(delta/sqrt(Var.delta)-Z))#correct

DM.new<-NULL
for(i in 1:dim(DM)[2]){
  
  DM.new<-cbind(DM.new,DM[,i], DM[,i])
}
calcPower.SWD(ThetaEst=delta, Design=DM.new[,3:12], 
              tauq=sigma.alpha^2, sigmaq=sigma^2, sigmaq.error =sigma.e^2,  
              noSub=K, type="longitudinal")
#Table 1, row 5
calcPower.SWD(ThetaEst=delta, Design=DM.new, 
              tauq=sigma.alpha^2, sigmaq=sigma^2, sigmaq.error =sigma.e^2,  
              noSub=K, type="longitudinal")

