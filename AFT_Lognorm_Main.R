##
##Demonstrate the application of POI-SIMEX method using the simulated data.
##In the real application, model statement needs to be modified accordingly.
###
# REFERENCE:  "POI-SIMEX for conditionally 
#    Poisson distributed biomarkers from tissue microarrays" (2024)
#    Aijun Yang (aijunyan@uvic.ca), Phineas T. Hamilton, Brad H. Nelson, Julian J. Lum, 
#    Mary Lesperance, Farouk S. Nathoo

##load required packages

library(survival)
library(survivalAnalysis)
library(SurvRegCensCov)

# Get the path of the current R file
current_file_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
# Set the working directory to the directory of the current file
setwd(current_file_path)

# if it is not working, simply hard code the path where the two R files are saved
# setwd("path_two_R_files_saved")
source("./POI_SIMEX_AFT.R")

progStart <- Sys.time()

# -------------------------------------------------------------- #
##Step 1 ------ simulate one data and with sample n=50
# paramter values:
n <- 50
a <- 1 
b <- 2
b0 <- 2 
b1 <- 1
b2 <- 0.5
sig.e <- 2

mydata <- sim.lognormT(a, b, n,
                       unif.a0 = 0.5,
                       unif.b0 = 9,
                       seed = 12345, 
                       rcensor = 0.2)

# -------------------------------------------------------------- #
##Step 2 ------ run the analysis

##for simulated data, fit true value

my.formula <- Surv(Y, delta) ~ lambda + Z ##lambda and Z are variables in mydata
true.fit <- survreg(formula = my.formula,
                    data = mydata,
                    dist='lognormal',
                    robust = TRUE)
true.est <- summary(true.fit)$table


## naive fit --------------#

my.formula <- Surv(Y, delta) ~ lambda.naive + Z  ##lambda.naive and Z are variables in mydata
naive.fit <- survreg(formula = my.formula,
                     data = mydata,
                     dist = 'lognormal',
                     robust = TRUE)
naive.est <- summary(naive.fit)$table


## simex correction --------------#

simexaft.est <- POI.simexaft(formula = my.formula,
                             data = mydata,
                             SIMEXvariable = "lambda.naive",
                             areaVariable = "A", 
                             B = 200,
                             lambda = seq(0, 2, 0.1),
                             extrapolation = "quadratic", 
                             dist="lognormal")

simex.est <- cbind(coef = simexaft.est[1]$coefficients, 
                   se = simexaft.est[2]$se,
                   scale = simexaft.est[3]$scalereg,
                   pvalue = simexaft.est[4]$pvalue)


# -------------------------------------------------------------- #
##Step 3 ----- save the result
##for simulated data, write true.est

write.csv(true.est, row.names = TRUE, file = "./true.est.csv")
write.csv(naive.est, row.names = TRUE, file = "./naive.est.csv")
write.csv(simex.est, row.names = TRUE, file = "./simex.est.csv")

progEnd <- Sys.time() - progStart

cat("run time \n")
print(progEnd)
