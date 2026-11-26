# An Integrated Population Model to Incorporate Spatio-Temporal Heterogeneity in Demographic Rates

Code and data for the spatially explicit integrated population model (sIPM) applied to Gray Catbird data, together with the simulation study.

## Repository contents

| Item | Description |
|---|---|
| [`packages_needed.R`](packages_needed.R) | Installs and loads all R packages used in the simulation study and case study. |
| [`Case Study/`](Case%20Study) | Scripts to reproduce the sIPM analysis of Gray Catbird data: model fitting, posterior summaries, spatial predictions and figures. Includes the cleaned data. See [`Case Study/README.md`](Case%20Study/README.md). |
| [`Simulation study/`](Simulation%20study) | Scripts to generate simulated datasets, fit the sIPM, evaluate model performance and reproduce all simulation results in the manuscript. |

## Installation

1. Run `packages_needed.R` to install the required CRAN and GitHub packages.
2. The BayesNSGP scripts need the development version of **BayesNSGP**. The release on CRAN does not yet include all the changes they use. `packages_needed.R` installs it, or you can install it directly:

   ```r
   remotes::install_github("danielturek/BayesNSGP", subdir = "BayesNSGP")
   ```

All scripts assume that the working directory is the root of this repository (the folder containing `Case Study` and `Simulation study`).

## Software

The nearest-neighbour Gaussian process (NNGP) used in the manuscript was first implemented as standalone R functions ([`Case Study/NNGP_functions.R`](Case%20Study/NNGP_functions.R)), and the manuscript results were produced with them.

These functions have since been developed, documented and tested as the R package [**BayesNSGP**](https://github.com/danielturek/BayesNSGP). It provides:

- neighbor construction,
- the NNGP density and simulation,
- a local-update MCMC sampler (`RW_NN_GP`),
- posterior prediction (`NNGP.pred`).

The package is described in Ketwaroo and Turek (submitted to the *Journal of Open Source Software).

The Case Study folder contains, in addition to the original scripts, a re-run of the case study with BayesNSGP and a comparison of run time and posterior estimates with the original implementation (see [Re-running the case study with BayesNSGP](Case%20Study/README.md#re-running-the-case-study-with-bayesnsgp)). The original scripts and results are unchanged.
