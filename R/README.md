# Analysis code

Every script here is run from the repository root, for example `Rscript R/figures/fig_ibgd_fit.R`. They need R (≥ 4.4) and [gstudio](https://github.com/dyerlab/gstudio) ≥ 1.16.0, plus `tidyverse`, `igraph`, `patchwork` and `ggrepel`. The case-study map also needs `sf` and `rnaturalearth`/`rnaturalearthhires`.

## What can be rebuilt from this repository

- **Every figure and table:** from the cached results in `data/derived/` and the summary tables in `data/`. These run in seconds to minutes and don't need the raw simulations.
- **The cached results themselves:** these are recomputed from the simulated Population Graphs and genotypes (about 50 lineages × 3,000 generations × several scenarios) and from the simulator that produced them. Neither is in this repository ==[simulation output and simulator: deposit location to be added]==. Point the scripts at a local copy with:

  ```sh
  export GG_SIM_DIR=/path/to/simulations/data   # simulated graphs and genotypes (replicate<k>/...)
  export GG_SIM_REPO=/path/to/simulator         # contains R/simulate.R, for simulation/fresh_replicates.R
  ```

  Several of these runs take hours and cache their work under `data/derived/` (ignored).

The *Araptus attenuatus* case study is fully reproducible from `gstudio::arapat`.

Reports that scripts write as they run go to `output/` (created on demand, not tracked).

## Layout

| Folder | Contents |
|---|---|
| `common/` | Shared I/O and definitions: loading simulated graphs, axis positions, lineage phases (`lineage_phases.R`), allele-frequency helpers. |
| `simulation/` | Scripts that read the raw simulations (`GG_SIM_DIR`) and write the cached results in `data/` and `data/derived/`. |
| `analysis/` | Summaries and tables computed from the cached results. |
| `case_study/` | *Araptus attenuatus*: the analysis, design-matched calibration and all robustness checks. |
| `figures/` | One script per manuscript figure. |
| `legacy/` | The original code behind five cached inputs that the current scripts still read (see `legacy/README.md`). It needs helpers retired at gstudio 1.15, so it documents how those files were made rather than rerunning them. |

## Run order (from the raw simulations)

1. `simulation/extract_forward_scenarios.R` and `simulation/extract_information_null.R` produce `data/forward_scenarios*.csv` and `data/information_null*.csv`. Run each twice, once as is (γ = 1) and once with `GAMMA=0.5`.
2. `simulation/he_gradient_boundary.R`, then `analysis/lineage_phase_windows.R` and `analysis/phase_boundaries.R`, produce the lineage phases (`data/derived/lineage_phases.csv`, read through `common/lineage_phases.R`).
3. `simulation/half_bandwidth_compute.R` computes the per-census statistics at γ = ½ and γ = 1. `analysis/half_bandwidth_results.R` then summarizes them and writes the γ = ½ inputs of the figures.
4. Run these independently:
   - `analysis/phase_core_rerun.R`
   - `analysis/core_tables_auc.R`
   - `analysis/fpr_cluster_ci.R`
   - `analysis/bandwidth_gamma_sweep.R`
   - `simulation/fresh_replicates.R`, then `analysis/fresh_replicates_report.R`
   - `simulation/permutation_nulls_1145.R`, then `analysis/permutation_nulls_1145_report.R`
   - `simulation/density_power_full.R` and `simulation/density_power_hs.R`, then `analysis/density_power_results.R`
5. Case study (no simulations needed): start with `case_study/case_study_arapat.R`, then run the remaining `case_study/` scripts.
6. Figures: run the scripts in `figures/` with `GG_TAG=g05`, which selects the γ = ½ inputs.

## Figures

| Figure | File | Script |
|---|---|---|
| 1 | `media/fig-gravity-field-bandwidth-v4.png` | `figures/fig_gravity_field_bandwidth.R` |
| 2 | `media/fig-scenarios-combined-v5.png` | `GG_TAG=g05 figures/fig_scenarios_combined.R` |
| 3 | `media/fig-source-sink-score-v7.png` | `GG_TAG=g05 figures/fig_source_sink_score.R` |
| 4 | `media/fig-ibgd-fit-v7.png` | `GG_TAG=g05 figures/fig_ibgd_fit.R` |
| 5 | `media/fig-gravity-w-vs-divmigrate-3x3-v8.png` | `GG_TAG=g05 figures/fig_gravity_vs_divmigrate.R` |
| 6 | `media/fig-case-study-arapat-v6.png` | `figures/fig_case_study_arapat.R` |
| S1 | `media/fig-information-null-v3.png` | `GG_TAG=g05 figures/fig_information_null.R` |
| S2 | `media/fig-isotropic-baseline-v2.png` | `GG_TAG=g05 figures/fig_isotropic_baseline.R` |
| S3 | `media/fig-signed-asymmetry-v6.png` | `GG_TAG=g05 figures/fig_signed_asymmetry.R` |
| S4 | `media/fig-source-sink-examples-v3.png` | `figures/fig_source_sink_examples.R` |
| S5 | `media/fig-bandwidth-gamma-sweep-v4.png` | `analysis/bandwidth_gamma_sweep.R` |
| S6 | `media/fig-sampling-intensity-v3.png` | `figures/fig_sampling_intensity.R` |

## Tables

| Table | Script |
|---|---|
| 1 | Simulation parameters (no computation) |
| 2, S5 (source–sink test, false positives) | `analysis/half_bandwidth_results.R` |
| 3, S7 (IBGD fits, directional isolation) | `analysis/half_bandwidth_results.R`, `analysis/phase_core_rerun.R` |
| S1 | Notation (no computation) |
| S2 (identities) | `legacy/node_potential_step0.R` |
| S3, S4 (decomposition of the late sign change) | `analysis/half_bandwidth_results.R` |
| S6 (source–sink detection by phase) | `analysis/phase_core_rerun.R` |
| S8 (lineage-clustered intervals) | `analysis/fpr_cluster_ci.R` |
| S9, S10 (test agreement; bandwidth sensitivity of Δ) | `analysis/core_tables_auc.R` |
| S11 (source–sink test by bandwidth) | `analysis/phase_core_rerun.R` |
| S12 (second set of lineages) | `simulation/fresh_replicates.R`, `analysis/fresh_replicates_report.R` |
| S13 (sampling intensity) | `analysis/density_power_results.R`; S10.4 numbers from `simulation/density_power_hs.R` |
| S14 (*Araptus* calibration) | `case_study/case_study_arapat_calibration.R`; generator check `analysis/popgraph_null_validation.R`; simulated design `case_study/araptus_design_calibration.R` |
| S15 (*Araptus* bandwidth sweep) | `case_study/araptus_gamma_sweep.R` |
| S16 (*Araptus* calibration at γ = 1) | `case_study/araptus_calibration_gamma1.R` |
| S17–S20 (information content, locus jackknife, northern calibration) | `case_study/araptus_locus_jackknife.R`, `case_study/araptus_northern_calibration.R` |
| S21 (nulls with diversity differences) | `case_study/araptus_heteroNe_null.R`, `case_study/araptus_degree_structure.R` |
| S22 | Summary of the case-study checks (no computation) |

Permutation nulls for Δ (Section S5.4): `simulation/permutation_nulls_1145.R` and `analysis/permutation_nulls_1145_report.R`. Slope-ratio AUC (Section S6.4): `analysis/core_tables_auc.R`.

`common/gravity_convention.R`, run on its own (`Rscript R/common/gravity_convention.R`), checks that the stock gstudio functions reproduce the archived edge values. It uses `data/derived/node_potential_followup.rda` and `node_potential_nodes.rds`.
