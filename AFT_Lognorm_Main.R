##
##Demonstrate the application of POI-SIMEX method using the simulated data
##In the real application, model statement needs to be modified accordingly
###

library(tidyverse);library(survival);library(survivalAnalysis)
library(purrr);library(dplyr);library(plyr);library(stringr)

setwd("C:/Users/workingDir")
source("./POI_SIMEX_AFT.R")

sim.lognormT<-function(a,b,n,unif.a0=0.5,unif.b0=9,seed=12345,rcensor=0){
  set.seed(seed)
  simdata<-c()
  #area in POI process--set to 1 for simplicity
  A<-rep(1,n)
  ##covariate without ME
  Z<-runif(n,unif.a0,unif.b0)
  
  ##a=gamma shape b=gamma scale
  
  lambda<-rgamma(n,a,b)
  
    lambda<-rgamma(n,shape=a,scale=b)  #rate=1/scale
    
    W<- rpois(n, lambda*A) #generate Poissons, mis-measured covariate
    lambda.naive<-W/A
    
    log.Y<-b0+b1*lambda+b2*Z+sig.e*rnorm(n)  ##set standard normal dist rnorm(0,1)
    Y<-exp(log.Y)
    delta<- rep(1, n) #event indicator
    delta[sample(n, rcensor*n)] <- 0 # random censorship
    
    simdata<-data.frame(Y=Y,Z=Z,lambda=lambda,lambda.naive=lambda.naive,A=A,delta=delta)
  return (simdata)
}


progStart=Sys.time()

##step 1 simulate one data and with sample n=50
n<-50;
a=1;b=2;b0=2;b1=1;b2=0.5;sig.e=2

mydata<-sim.lognormT(a,b,n,unif.a0=0.5,unif.b0=9,seed=12345,rcensor=0.2)

##step 2 run the analysis

naive.est<-c(); simex.est<-c()

    ##naive fit

    formula<-Surv(Y,delta) ~ lambda.naive + Z
    
    naive.fit = survreg(formula=formula,data=mydata,dist='lognormal',robust=T)
    naive.est0<-summary(naive.fit)$table
    para<-row.names(naive.est0)
    naive.est<-cbind(naive.est0,para)
    
    ##simex correction
    
    simexaft.est<-POI.simexaft(formula=formula,data=mydata,SIMEXvariable="lambda.naive",areaVariable="A",repind=list(),B=200,lambda=seq(0,2,0.1),extrapolation="quadratic",dist="lognormal")
    simex0<-cbind(simexaft.est[1]$coefficients,simexaft.est[2]$se,simexaft.est[3]$scalereg,
                  simexaft.est[4]$pvalue)
    colnames(simex0)<-c("coef","se","scale","pvalue")
    simex.est<-cbind(para=row.names(simex0),simex0)

##save the result
write.csv(naive.est,row.names=F,file="./output/naive.est.csv")
write.csv(simex.est,row.names=F,file="./output/simex.est.csv")

progEnd=Sys.time()-progStart

cat("run time \n")
print(progEnd)
