##
##Demonstrate the application of POI-SIMEX method using the simulated data.
##In the real application, model statement needs to be modified accordingly.
###
# REFERENCE:  "POI-SIMEX for conditionally 
#    Poisson distributed biomarkers from tissue microarrays" (2024)
#    Aijun Yang (aijunyan@uvic.ca), Phineas T. Hamilton, Brad H. Nelson, Julian J. Lum, 
#    Mary Lesperance, Farouk S. Nathoo

##install required packages
#Package plyr must be loaded before dplyr
list.of.packages <- c("tidyverse", "survival","survivalAnalysis","SurvRegCensCov","purrr","plyr","dplyr","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE) #load the packages

# Get the path of the current R file
current_file_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
 # Set the working directory to the directory of the current file
setwd(current_file_path)

# if it is not working, simply hard code the path where the two R files are saved
# setwd("path_two_R_files_saved")

source("./POI_SIMEX_AFT.R")

progStart=Sys.time()

##Step 1 ------ simulate one data and with sample n=50
n<-50;
a=1;b=2;b0=2;b1=1;b2=0.5;sig.e=2

mydata<-sim.lognormT(a,b,n,unif.a0=0.5,unif.b0=9,seed=12345,rcensor=0.2)

##Step 2 ------ run the analysis

    ##for simulated data, fit true value
    true.est<-c()
    my.formula<-Surv(Y,delta) ~ lambda + Z ##lambda and Z are variables in mydata
    true.fit = survreg(formula=my.formula,data=mydata,dist='lognormal',robust=T)
    true.est0<-summary(true.fit)$table
    para<-row.names(true.est0)
    true.est<-cbind(true.est0,para)
  
    ##naive fit
    naive.est<-c()
    my.formula<-Surv(Y,delta) ~ lambda.naive + Z  ##lambda.naive and Z are variables in mydata
    naive.fit = survreg(formula=my.formula,data=mydata,dist='lognormal',robust=T)
    naive.est0<-summary(naive.fit)$table
    para<-row.names(naive.est0)
    naive.est<-cbind(naive.est0,para)

   
    ##simex correction
    simex.est<-c()
    simexaft.est<-POI.simexaft(formula=my.formula,data=mydata,SIMEXvariable="lambda.naive",
                               areaVariable="A",B=200,
                               lambda=seq(0,2,0.1),extrapolation="quadratic",dist="lognormal")

    simex0<-cbind(simexaft.est[1]$coefficients,simexaft.est[2]$se,simexaft.est[3]$scalereg,
                  simexaft.est[4]$pvalue)

    colnames(simex0)<-c("coef","se","scale","pvalue")
    simex.est<-cbind(para=row.names(simex0),simex0)

##Step 3 ----- save the result
##for simulated data, write true.est
write.csv(true.est, row.names=F,file="./true.est.csv")
write.csv(naive.est,row.names=F,file="./naive.est.csv")
write.csv(simex.est,row.names=F,file="./simex.est.csv")

progEnd=Sys.time()-progStart

cat("run time \n")
print(progEnd)
