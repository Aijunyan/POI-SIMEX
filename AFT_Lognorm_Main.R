##
##Demonstrate the application of POI-SIMEX method using the simulated data
##In the real application, model statement needs to be modified accordingly
###
##install required packages

list.of.packages <- c("tidyverse", "survival","survivalAnalysis","SurvRegCensCov","purrr","dplyr","plyr","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

setwd("~")
source("./POI_SIMEX_AFT.R")

lapply(list.of.packages, require, character.only = TRUE)

sim.lognormT<-function(a,b,n,unif.a0=0.5,unif.b0=9,seed=12345,rcensor=0){
  # Arguments
  # @para a, b are gamma shape and scale parameters to simulate true covariate X
  # @para n, the number of observations
  # @para unif.a0, unif.b0 are two parameters of uniform distribution to simulate Z covariate
  # @para seed: random seed
  # @para rcensor: % of censoring, values from 0 to 1
  # output dataset: simdata, type is dataframe
  
  set.seed(seed)
  simdata<-c()
  #area in POI process--set to 1 for simplicity
  A<-rep(1,n)
  ##covariate without ME
  Z<-runif(n,unif.a0,unif.b0)
  
  ##a=gamma shape b=gamma scale
  
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

    ##for simulated data, fit true value
    true.est<-c()
    my.formula<-Surv(Y,delta) ~ lambda + Z
    true.fit = survreg(formula=my.formula,data=mydata,dist='lognormal',robust=T)
    true.est0<-summary(true.fit)$table
    para<-row.names(true.est0)
    true.est<-cbind(true.est0,para)
  
    ##naive fit
    naive.est<-c()
    my.formula<-Surv(Y,delta) ~ lambda.naive + Z
    
    naive.fit = survreg(formula=my.formula,data=mydata,dist='lognormal',robust=T)
    naive.est0<-summary(naive.fit)$table
    para<-row.names(naive.est0)
    naive.est<-cbind(naive.est0,para)

   
    ##simex correction
    simex.est<-c()
    simexaft.est<-POI.simexaft(formula=my.formula,data=mydata,SIMEXvariable="lambda.naive",areaVariable="A",B=200,
                               lambda=seq(0,2,0.1),extrapolation="quadratic",dist="lognormal")

    simex0<-cbind(simexaft.est[1]$coefficients,simexaft.est[2]$se,simexaft.est[3]$scalereg,
                  simexaft.est[4]$pvalue)

    colnames(simex0)<-c("coef","se","scale","pvalue")
    simex.est<-cbind(para=row.names(simex0),simex0)

##Step 3 save the result
##if true.est exist
write.csv(true.est, row.names=F,file="./true.est.csv")
write.csv(naive.est,row.names=F,file="./naive.est.csv")
write.csv(simex.est,row.names=F,file="./simex.est.csv")

progEnd=Sys.time()-progStart

cat("run time \n")
print(progEnd)
