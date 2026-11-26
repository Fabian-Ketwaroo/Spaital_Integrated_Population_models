# Simulation study

This folder contains the scripts used to conduct the simulation study presented in the manuscript. The simulations explore the performance of the spatially explicit integrated population model (sIPM).

## Main simulation study

| Script | Scenario |
|---|---|
| `Perfect_alignment.R` | Count and capture-recapture data are perfectly aligned in space and time. |
| `Spatial_misalignment.R` | Count and capture-recapture data are spatially misaligned, without temporal misalignment. |
| `SC_Spatio_temporal_misalignment_high.R` | Spatial misalignment and **high** temporal misalignment, with a spatially **constant** covariate effect. |
| `SC_Spatio_temporal_misalignment_low.R` | Spatial misalignment and **low** temporal misalignment, with a spatially **constant** covariate effect. |
| `SV_Spatio_temporal_misalignment_high.R` | Spatial misalignment and **high** temporal misalignment, with a spatially **varying** covariate effect. |
| `SV_Spatio_temporal_misalignment_low.R` | Spatial misalignment and **low** temporal misalignment, with a spatially **varying** covariate effect. |

## Sampling intensity

[`Sampling_intensity/`](Sampling_intensity) contains the additional experiment reported in Section S1.6 of the Supporting Information. It repeats the spatial misalignment scenario at four levels of sampling intensity (M = 80, 200, 500 and 800 locations used for model fitting). See the README in that folder for design details, run instructions and computation times.

## Recommended execution order

1. `Perfect_alignment.R`: establishes baseline performance.
2. `Spatial_misalignment.R`: explores purely spatial misalignment.
3. SC / SV scripts: explore combinations of spatial and temporal misalignment with different covariate specifications.
4. `Sampling_intensity/`: explores the effect of the number of sampling locations.

## Notes

All scripts assume the dependencies in `packages_needed.R` have been installed. Expected computation times for each script are reported in Section S1.3 and Table S8 of the Supporting Information.

## Note on MCMC settings

Each script uses the final 1000 posterior samples per chain for prediction, selected by explicit row indices. Two MCMC configurations are used across the simulation study, and the indices differ accordingly.

| Scripts | `niter` | `nburnin` | `thin` | `thin2` | Primary indices | Secondary indices |
|---|---|---|---|---|---|---|
| `Perfect_alignment.R`, `Spatial_misalignment.R`, `SC_Spatio_temporal_misalignment_high.R`, `SC_Spatio_temporal_misalignment_low.R`, `Sampling_intensity/M80.R`, `M200.R`, `M500.R`, `M800.R` | 85000 | 60000 | 10 | 20 | `1501:2500` | `251:1250` |
| `SV_Spatio_temporal_misalignment_high.R`, `SV_Spatio_temporal_misalignment_low.R` | 125000 | 80000 | 15 | 25 | `2001:3000` | `801:1800` |

The first configuration retains (85000 − 60000)/10 = 2500 primary and (85000 − 60000)/20 = 1250 secondary samples per chain. The second retains 45000/15 = 3000 and 45000/25 = 1800. In each case the indices select the last 1000 rows.

The spatially varying scripts use longer chains because the additional regression coefficients mix more slowly.

> [!IMPORTANT]
> If you change `niter`, `nburnin`, `thin` or `thin2`, these row indices will no longer be valid. R returns `NA` rather than an error for out-of-range row indices, so the scripts will run to completion but produce `NA` or invalid predictions with no warning. If you alter the MCMC settings, update both index ranges to select the final 1000 rows of each chain.
