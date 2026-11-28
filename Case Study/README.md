# Case study: Gray Catbird

This folder contains all scripts and processed data required to reproduce the case study analysis presented in the manuscript.

We apply the proposed spatially explicit integrated population model (sIPM) to jointly analyse:

- count data from the North American Breeding Bird Survey (BBS),
- capture–recapture (CR) data from the Monitoring Avian Productivity and Survivorship (MAPS) program.

We consider BBS and MAPS sampling locations along the eastern coast of North America from 2004–2014.

## Contents

- [Data](#data)
- [Original analysis (manuscript)](#original-analysis-manuscript)
- [Re-running the case study with BayesNSGP](#re-running-the-case-study-with-bayesnsgp)
- [Recommended execution order](#recommended-execution-order)
- [Computational note](#computational-note)

## Data

| File / folder | Description |
|---|---|
| `GRCA_data_2004_2014.Rdata` | Processed data required to reproduce the analysis (objects listed below). |
| `Data_csv/` | The same data as `GRCA_data_2004_2014.Rdata` in csv form. csv cannot store 3-dimensional arrays, so `Fm` and `Rm` are saved per sampling location. |
| `bcr_terrestrial_shape/` | Shapefiles of Bird Conservation Regions (BCRs) used for mapping and visualization. |

Objects in `GRCA_data_2004_2014.Rdata`:

| Object | Description |
|---|---|
| `C` | BBS count spatio-temporal data |
| `obs` | observer identity |
| `oby1` | first-time observer covariate: 1 = first-time observer |
| `Fm` | spatio-temporal m-array for a mixture of transient and residual individuals |
| `Rm` | spatio-temporal m-array for residual individuals |
| `visits` | binary covariate: 1 = sampling location was sampled, 0 = not sampled |
| `GRCA_splocations` | spatial coordinates of BBS and MAPS sampling locations |
| `Tmax.resid` | spatio-temporal values of residual maximum temperature |

## Original analysis (manuscript)

These scripts produced the results in the manuscript.

### Scripts

| File | Description |
|---|---|
| `NNGP_functions.R` | Functions required to implement the nearest-neighbour Gaussian process (NNGP). |
| `GRCA_SV_NNGP.R` | Analysis with an NNGP and a **spatially varying** covariate effect on residual maximum temperature. **This is the primary model used in the manuscript.** Also includes code to perform model checking. |
| `GRCA_SV_No_NNGP.R` | Spatially varying covariate effect on residual maximum temperature, **without** the NNGP. |
| `GRCA_SC_NNGP.R` | NNGP and a **constant** covariate effect on residual maximum temperature. |
| `GRCA_SC_No_NNGP.R` | Spatially constant covariate effect on residual maximum temperature, **without** the NNGP. |
| `GOF.R` | Goodness of fit for the selected model (`GRCA_SV_NNGP.R`). |
| `GRCA_predictions_phi_f.R` | Spatially continuous predictions of survival and recruitment averaged over 2004–2014, from `GRCA_SV_NNGP.R`. Needs `GRCA_main_samples_NNGP_SV.Rdata` and `GRCA_samples2_NNGP_SV.Rdata`. |
| `GRCA_predictions_N_gr.R` | Spatially continuous predictions of population growth rate and relative abundance averaged over 2004–2014, from `GRCA_SV_NNGP.R`. Needs `GRCA_samples2_NNGP_SV.Rdata`. |
| `Case study workflow.pdf` | Visual roadmap of the execution steps and key researcher decisions. It mirrors the sequence in the manuscript's Methods section and guides navigation and modification of the R scripts. |

### MCMC output

| File | Parameters |
|---|---|
| `GRCA_main_samples_NNGP_SV.Rdata` | MCMC samples from `GRCA_SV_NNGP.R`: `rf`, `w`, `alpha.p`, `eta`, `wst`, `wft`, `rho.as`, `alpha2.as`, `rho.af`, `alpha2.af`, `beta.f.tmax`, `beta.s.tmax`, `beta.f.long`, `beta.s.long`, `beta.f.lat`, `beta.s.lat`, `beta.f.tmaxlong`, `beta.s.tmaxlong`, `beta.f.tmaxlat`, `beta.s.tmaxlat` |
| `GRCA_samples2_NNGP_SV.Rdata` | MCMC samples from `GRCA_SV_NNGP.R`: `C.pred`, `Fm.pred`, `Rm.pred`, `N`, `sa`, `f`, `growth.rate`, `w.af`, `w.as`. Available at <https://doi.org/10.5281/zenodo.22029221>. |

## Re-running the case study with BayesNSGP

The manuscript results were produced with the standalone NNGP functions in `NNGP_functions.R` and scalar random-walk (`"RW"`) samplers on the spatial effects (`GRCA_SV_NNGP.R`). The scripts below repeat the analysis with the [BayesNSGP](https://github.com/danielturek/BayesNSGP) package, using its local-update sampler (`RW_NN_GP`) and its prediction function (`NNGP.pred`), and compare the two. They do not replace the original scripts.

**Requirements:** the development version of BayesNSGP (see the [main README](../README.md#installation)) and the packages in `packages_needed.R`. Run all scripts from the repository root.

### Files

| File | Description |
|---|---|
| `GRCA_SV_BayesNSGP.R` | Same model and data as `GRCA_SV_NNGP.R`, with the NNGP functions from BayesNSGP. Set `sampler_type` to `"RW_NN_GP"` (BayesNSGP sampler) or `"RW"` (NIMBLE scalar random walk, as in `GRCA_SV_NNGP.R`) and run the script once for each, in a fresh R session. Saves `sipm_<sampler_type>.rds` with the samples, run time and MCMC settings. |
| `GRCA_compare_samplers.R` | Reads the two results files and reports run time, effective sample size per second, the agreement of the identifiable quantities (survival `sa`, recruitment `f`, hyperparameters, fixed effects) and the agreement of the survival spatial field `w.as`, which is judged against the agreement between two chains of the same run. Writes `comparison_summary.csv` and two figures. Stops if the two runs used different MCMC settings. |
| `GRCA_predictions_phi_f_BayesNSGP.R` | Spatially continuous predictions of survival and recruitment, as `GRCA_predictions_phi_f.R`, but with `NNGP.pred()` from BayesNSGP. Uses `sipm_RW_NN_GP.rds`. |
| `sipm_RW.rds`, `sipm_RW_NN_GP.rds` | Output of `GRCA_SV_BayesNSGP.R` for the two sampler types. [Add: the spatial effects (`samples2`) are large; state whether they are included or where they are archived, e.g. Zenodo.] |

### Settings

Both runs use identical settings, initial values and model:

```r
runMCMC(FullModelmod,
        niter = 350000, nburnin = 200000, thin = 30, thin2 = 100,
        samplesAsCodaMCMC = TRUE, nchains = 2,
        inits = list(inits1, inits2), setSeed = c(1, 2), summary = TRUE)
```

Software and hardware (from `sessionInfo()`, which is stored inside each `.rds` file as `$session`):
 
| | |
|---|---|
| R | 4.4.2 |
| nimble | 1.4.1 |
| BayesNSGP | development version from GitHub (package version 0.3.1; commit [https://github.com/danielturek/BayesNSGP/pull/14] [https://github.com/danielturek/BayesNSGP/pull/15] |
| sf | 1.0.19 |
| Platform | aarch64-apple-darwin20 (macOS), single core (Apple M4 Max) |
 
Each run used one core (NIMBLE runs the two chains one after the other), and nothing else was running on the machine during the timing.

### Differences from the original scripts

1. **Distances.** `BayesNSGP::computeNeighbors()` computes distances on the WGS84 ellipsoid, as before. The script converts them to kilometres, so `rho` and its prior are on the same scale as in the manuscript.
2. **Initial values.** Initial values of the spatial effects are simulated with `BayesNSGP::rmnorm_NN_GP()` instead of a dense multivariate normal draw.
3. **Prediction.** The neighbor search is done once per prediction location and all posterior draws are predicted together. The spatial effects and hyperparameters used for one prediction draw come from the same MCMC iteration (rows of `samples` and `samples2` are matched by iteration number, because they are thinned differently).
4. **Locations west of all sampling locations.** These have no preceding neighbors and are drawn from the marginal distribution of the spatial effect.

[Optional: state whether and how 3 and 4 change the maps in the manuscript, once checked.]

### Results of the comparison

Output of `GRCA_compare_samplers.R` (28 November 2026). Both runs use the settings above.

**Run time and efficiency**

| | RW | RW_NN_GP |
|---|---|---|
| Run time (hours) | 20.0 | 9.0 |
| Speed-up | | 2.2× |
| Effective sample size per second, ratio RW_NN_GP / RW | | median 2.4 (range 1.2 to 4.8 over monitored parameters) |

**Agreement of the identifiable quantities.** Both samplers give the same posterior for the quantities that the data identify:

- Posterior means of survival (`sa`) and recruitment (`f`) are almost identical: correlation 0.983 for `sa` (7590 nodes) and 0.984 for `f` (6970 nodes).
- The hyperparameters agree between the two runs (for example `rho.as` is about 900 km in both runs, `alpha2.as` about 0.08).

**Survival spatial field (`w.as`).** The posterior means of `w.as` agree only up to Monte Carlo error. This is a limit of the model, not of the sampler:

| Correlation of posterior means of `w.as` | |
|---|---|
| Between the two runs | 0.73 |
| Between chain 1 and chain 2 of RW | 0.62 |
| Between chain 1 and chain 2 of RW_NN_GP | 0.84 |

The between-run value lies between the two within-run values, so the two samplers differ no more than two chains of the same sampler. The spatial field is only weakly identified. Survival at a site is the sum of a temporal intercept, the covariate terms (including linear longitude and latitude) and `w.as`, and with very few capture-recapture locations the smooth part of the field trades off with the intercept and the longitude/latitude terms. This is visible in the per-chain posterior means of `wst[1]`, which differ more between chains than between samplers. The posterior of `rho.as` is also close to the upper end of its prior (`dmax/3`, about 990 km), that is, the data favour a very smooth survival field.

**Mixing.** The convergence diagnostics are similar for the two runs and are limited by the model:

| | RW | RW_NN_GP |
|---|---|---|
| Largest R-hat, monitored parameters | 1.35 | 1.19 |
| `w.as` nodes: median effective sample size | 166 | 153 |
| `w.as` nodes: share with R-hat above 1.1 | 0.18 | 0.14 |
| `w.as` nodes: largest R-hat | 2.82 | 1.38 |

Posterior means of the monitored parameters of the two runs differ by less than 2 Monte Carlo SE for 75% of the parameters. Because the effective sample size ignores differences between chains, the Monte Carlo SEs are too small where R-hat is large, so this share is below the 95% expected under identical posteriors even though the two samplers target the same distribution.

**Summary.** The BayesNSGP sampler gives the same estimates of survival, recruitment and the hyperparameters as the scalar random-walk sampler, 2.2 times faster. The survival spatial effects themselves are only weakly identified and, with this chain length, not precisely estimated by either sampler; R-hat is above 1.1 for some nodes in both runs.

## Recommended execution order

**Original analysis**

1. Run `GRCA_SV_NNGP.R` (or an alternative model script).
2. Assess model fit with `GOF.R`.
3. Generate spatial predictions with `GRCA_predictions_phi_f.R` and `GRCA_predictions_N_gr.R`.

**BayesNSGP re-run**

1. Run `GRCA_SV_BayesNSGP.R` with `sampler_type = "RW_NN_GP"`.
2. Run `GRCA_SV_BayesNSGP.R` with `sampler_type = "RW"`, in a fresh R session on the same machine with nothing else running, so that run times are comparable.
3. Run `GRCA_compare_samplers.R`.
4. Run `GRCA_predictions_phi_f_BayesNSGP.R`.

## Computational note

All scripts assume the dependencies in `packages_needed.R` have been installed.

To skip long model runs, the pre-computed MCMC samples from `GRCA_SV_NNGP.R` are provided in `GRCA_main_samples_NNGP_SV.Rdata` and `GRCA_samples2_NNGP_SV.Rdata` (the latter on Zenodo, see above). You can skip step 1 and go straight to step 3 to generate the spatial predictions.

