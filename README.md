# Experiments for "Pairwise Stochastic Approximation for Confirmatory Factor Analysis of Categorical Data"

This repository collects all the code needed to reproduce results reported in "Pairwise Stochastic Approximation for Confirmatory Factor Analysis of Categorical Data".

- In the `customPackage` folder, it is provided the `.tar.gz` source file to install the `plFA` package, which contains the custom estimation routines described in the paper.


- The `scripts` folder collects the three `.R` files used to produce the results in the paper:
    - `sims.R`
    - `big5.R`
    - `capabilities.R`
    
    
- The `data` folder contains the dataset used for the learning capabilities application. The one related to the Big Five survey, instead, will be automatically downloaded by the `big5.R` script.

Running any of the three `.R` files will create an `output` folder where the results will be saved.
