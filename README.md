# Project Title

POI-SIMEX for conditionally Poisson distributed biomarkers from tissue microarrays

## Description
Motivated by cancer studies where biomarkers obtained from a tissue microarray subsampled from a large tumor tissue, we need to correct the bias in the regression coefficient estimation when using observed cell densities as surrogates for the true tissue cell densities. We extend SIMEX to the conditional Poisson case (POI-SIMEX) where measurement errors are non-Gaussian with heteroscedastic variance.

## Getting Started

### Dependencies

* R libraries: tidyverse, survival, survivalAnalysis, SurvRegCensCov, purrr, 
plyr, dplyr, stringr
### Installing

*  Download or copy the R programs  POI-SIMEX_AFT.R and AFT_Lognorm_Main. R
*  Modify the working directory depending on where the two R files are saved:
* current_file_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
* setwd(current_file_path)
* if it is not working, simply hard code the path where the two R files are saved
* setwd("path_two_R_files_saved")
* source("./POI_SIMEX_AFT.R")
* modify the output path if needed
* Write.csv(true.est,row.names=F,file="./true.est.csv")
* Write.csv(naive.est,row.names=F,file="./naive.est.csv")
* write.csv(simex.est,row.names=F,file="./simex.est.csv")

### Executing program

* After the modification, the program can run by section
* Block 1--Simulate the data (or use real data, real data needs to modify the model statment)
* Block 2--Run the models: using true value if data is simulated; using observed value for the naive method; and the POI-SIMEX method
* Block 3--Save the outputs

## Authors

Aijun Yang (aijunyan@uvic.ca)
Phineas T. Hamilton, Brad H. Nelson, Julian J. Lum, Mary Lesperance, Farouk S. Nathoo
## Reference
POI-SIMEX for Conditionally Poisson Distributed Biomarkers from Tissue Microarrays
## Version History

* 0.1
    * Initial Release

## Acknowledgments

* Xiong J, He W, Yi GY (2019). simexaft: R package version 1.0.7.1, [AFT SIMEX]( https://CRAN.R-project.org/package=simexaft)

