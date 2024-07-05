# POI-SIMEX
# POI-SIMEX_AFT.R is R code for correcting the measurement error for the conditionally Poisson distributed surrogates
# POI-SIMEX_AFT.R is modified based on SIMEX.AFT package source code
# Y--log survival time
# X--error-prone covariate(s)
# Z--covariate(s) without measurement error
# download the AFT_lognorm_Main.R and  POI_SIMEX_AFT.R to your preferred folder
# You might need to modify the working path and output path in  AFT_lognorm_Main.R 
# Check the required package
# The model statement needs to be modified according to the covariates included

# Project Title

POI-SIMEX for correcting measurement error in the conditionally Poisson Distributed Biomarkers from tissue microarrays

## Description
Motivated by cancer studies where biomarkers obtained from a tissue microarray subsampled from a large toumer tissue, we need to correct the bias in the regressio n coefficent estimation when using observed cell densities as surrogates for the true tissue cell densities. We extend SIMEX to the conditional Poisson case (POI-SIMEX) where measurement errors are non-Gaussian with heteroscedastic variance.

## Getting Started

### Dependencies

* R libraries: library(tidyverse);library(survival);library(survivalAnalysis);library(SurvRegCensCov)
library(purrr);library(dplyr);library(plyr);library(stringr)

### Installing

* How/where to download your program
* Any modifications needed to be made to files/folders

### Executing program

* How to run the program
* Step-by-step bullets
```
code blocks for commands
```

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors

Contributors names and contact info

ex. Dominique Pizzie  
ex. [@DomPizzie](https://twitter.com/dompizzie)

## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
