# Experiments for "Pairwise Stochastic Approximation for Confirmatory Factor Analysis of Categorical Data"

This repository collects all the code needed to reproduce results reported in "Pairwise Stochastic Approximation for Confirmatory Factor Analysis of Categorical Data".

- The estimation requires the installation of the development package `plFA` via
``` r
# install.packages("devtools")
devtools::install_github("giuseppealfonzetti/plFA@3b4e05b")
```
- The `scripts` folder collects the three `.R` main files used to reproduce the results in the paper:
    - `sims.R`
    - `big5.R`
    - `capabilities.R`
  Here you can also find additional experiments related to the stability of the algorithm:
    - `stepsizes.R`
    - `starting_points.R`
    
- The `data` folder contains the dataset used for the learning capabilities application. The one related to the Big Five survey, instead, will be automatically downloaded by the `big5.R` script.

Running any of the `.R` files will create an `output` folder where the results will be saved.
