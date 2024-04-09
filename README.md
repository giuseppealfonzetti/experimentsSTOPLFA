# Experiments for "Pairwise Stochastic Approximation for Confirmatory Factor Analysis of Categorical Data"

This repository collects all the code needed to reproduce results reported in "Pairwise Stochastic Approximation for Confirmatory Factor Analysis of Categorical Data".

- The estimation requires the installation of the development package `plFA` via
``` r
# install.packages("devtools")
devtools::install_github("giuseppealfonzetti/plFA@3b4e05b")
```
- The `scripts` folder collects the `.R` files used to reproduce the results in the paper:
    - `sims.R`: main simulations
    - `big5.R`: big five real data experiment
    - `capabilities.R`: capabilities real data experiment
    - `stepsizes.R`: additional simulation experiments
    - `starting_points.R` : additional simulation experiments
    
- The `data` folder contains the dataset used for the learning capabilities application. The one related to the Big Five survey, instead, will be automatically downloaded by the `big5.R` script.

Running any of the `.R` files will create an `output` folder where the results will be saved.
