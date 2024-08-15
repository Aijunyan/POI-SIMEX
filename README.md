# Project Title

POI-SIMEX for correcting measurement error in the conditionally Poisson Distributed Biomarkers from tissue microarrays

## Description
Motivated by cancer studies where biomarkers obtained from a tissue microarray subsampled from a large toumer tissue, we need to correct the bias in the regression coefficent estimation when using observed cell densities as surrogates for the true tissue cell densities. We extend SIMEX to the conditional Poisson case (POI-SIMEX) where measurement errors are non-Gaussian with heteroscedastic variance.

## Getting Started

### Dependencies

* R libraries: library(tidyverse);library(survival);library(survivalAnalysis);library(SurvRegCensCov)
library(purrr);library(dplyr);library(plyr);library(stringr)

### Installing

*  Download or copy the R programs  POI-SIMEX_AFT.R and AFT_Lognorm_Main. R
*  Modification of the working directory dependoning where the R codes are saved in AFT_Lognorm_Main.R and POI_SIMEX_AFT.R
*  setwd("~")
*  source("./POI_SIMEX_AFT.R")
*  modify the output path if needed
* Write.csv(true.est,row.names=F,file="./true.est.csv")
* Write.csv(naive.est,row.names=F,file="./naive.est.csv")
* write.csv(simex.est,row.names=F,file="./simex.est.csv")

### Executing program

* After the modification, the program can run by block
* Block 1--Simulate the data (or use real data, real data needs to modify the model statment)
* Block 2--Run the models using true value for simulated data, naive method, and POI-SIMEX method
* Block 3--Save the outputs

## Authors

Aijun Yang (aijunyan@uvic.ca)
Phineas T. Hamilton, Brad H. Nelson, Julian J. Lum Mary Lesperance, Farouk S. Nathoo
## Reference
POI-SIMEX for Conditionally Poisson Distributed Biomarkers from Tissue Microarrays
## Version History

* 0.1
    * Initial Release

## Acknowledgments

* [AFT SIMEX]( https://CRAN.R-project.org/package=simexaft)

