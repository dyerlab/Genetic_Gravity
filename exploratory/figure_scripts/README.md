# Figure and data scripts

Run from the manuscript repository root unless noted. `GG_TAG=g05` selects the γ = ½ inputs written by `exploratory/half_bandwidth_results.R`.

| Script | Output |
|---|---|
| `scenarios_combined_v2.R` | Fig_ScenarioComparison → `media/fig-scenarios-combined-v5.png` (with `GG_TAG=g05`) |
| `ibgd_fit_v2.R` | Fig_IBGDFit → `media/fig-ibgd-fit-v7.png` (with `GG_TAG=g05`; v5 = pGD/Isotropic, v6 = dcGD/Isotropic) |
| `gravity_w_vs_divmigrate_3x3_v3.R` | Fig_GravityVsDivMigrate → `media/fig-gravity-w-vs-divmigrate-3x3-v8.png` (with `GG_TAG=g05`) |
| `sampling_intensity_v3.R` | Figure S6 (Sampling intensity) → `media/fig-sampling-intensity-v3.png` (v2 = sample-size x breaks and x grid; v1 = duplicated legend); input `data/derived/density_power_tables.rda` from `exploratory/density_power_results.R` |
| `export_truth_compare.R` | `data/divmigrate_truth_compare.csv` (input of `R/divmigrate_truth.R`) |
| `build_timecourse_nm.R` | Run from the private research repo: the Nm time course copied to `data/divmigrate_timecourse_timecourse.csv` |

Earlier figure versions and dropped analyses are in `exploratory/archive/figure_scripts/` (see its README).
