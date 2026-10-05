# Figure and data scripts (recovered from session scratchpad)

Saved 2026-09-27. 'Run from' is the working directory each script expects; scratchpad caches now go to `data/derived/` (manuscript repo) or `data/divMigrate/scratch_outputs/` (original repo).

| Script | Run from | Purpose / output |
|---|---|---|
| `scenarios_combined.R` | manuscript repo | Fig_ScenarioComparison: 3x2 facet of C_D-bar, diameter, mean |Delta| for Redistributed | Obstructed -> media/fig-scenarios-combined.png |
| `divmig_correlation_all.R` | manuscript repo | Fig_DivMigrateCorrelation: per-census Spearman rho(S, Nm), all replicates + median -> media/fig-divmig-correlation-all.png |
| `lag1_compute.R` | original repo | Compute lag-1 series for every replicate with saved matrices (step 1 of Fig_DivMigrateLag1) -> divMig.{rep}.lag1.rda |
| `lag1_plot.R` | manuscript repo | Lag-1 median traces, Nm vs S (step 2; also caches lag1_all.rds used by lag1_isoband.R) -> media/fig-divmig-lag1-all.png |
| `lag1_isoband.R` | manuscript repo | Fig_DivMigrateLag1: lag-1 medians with isotropic mean +/- 1 SD bands (needs lag1_all.rds from lag1_plot.R) -> media/fig-divmig-lag1-isoband.png |
| `gravity_w_vs_divmigrate_3x3.R` | manuscript repo | Manuscript 3x3: Spearman rho, RMSE, replicates present; w_ij vs divMigrate Nm -> media/fig-gravity-w-vs-divmigrate-3x3.png |
| `build_timecourse_nm.R` | original repo | Builds the Nm time-course table later copied to data/divmigrate_timecourse_timecourse.csv (provenance of the manuscript CSV) |
| `export_truth_compare.R` | manuscript repo | Exports per-census truth metrics (pGD similarity + Nm) -> data/divmigrate_truth_compare.csv (input of R/divmigrate_truth.R) |
| `n_edges.R` | manuscript repo | Mean edge count through time by scenario -> media/fig-n-edges.png |
| `truth_raw.R` | manuscript repo | Truth comparison, raw p_ij vs Nm, five metrics (caches truth_all.rds) -> media/fig-divmig-truth.png |
| `truth_sim_vs_raw.R` | manuscript repo | Truth comparison with pGD similarity; raw vs similarity vs Nm -> media/fig-divmig-truth-sim.png, fig-divmig-truth-raw-vs-sim.png |
| `truth_sim_overall.R` | manuscript repo | Overall rho + RMSE, pGD similarity vs Nm -> media/fig-divmig-truth-sim-overall.png |
| `truth_sim_rho_mae.R` | manuscript repo | Overall rho + MAE, pGD similarity vs Nm -> media/fig-divmig-truth-sim-rho-mae.png |
| `gravity_w_truth.R` | manuscript repo | Gravity scores 1/(1+p), w_ij, 1/(1+w) vs truth -> media/fig-gravity-w-truth.png |
| `gravity_w_vs_divmigrate.R` | manuscript repo | w_ij vs Nm, rho/RMSE/MAE -> media/fig-gravity-w-vs-divmigrate.png |
| `gravity_w_vs_divmigrate_present.R` | manuscript repo | As above + replicates-present row -> media/fig-gravity-w-vs-divmigrate-present.png |
| `divmigrate_burnin_significance.R` | manuscript repo | [dropped thread] rep-39 burn-in divMigrate significance profile (OUR BH z-test) -> media/fig-divmigrate-burnin-significance.png |
| `divmigrate_significance_timeline.R` | manuscript repo | [dropped thread] rep-39 full-timeline divMigrate significance (OUR BH z-test) -> media/fig-divmigrate-significance-timeline.png |
| `divmigrate_first_significance_fig.R` | manuscript repo | [dropped thread] time to first significance, 50 reps (OUR BH z-test) -> media/fig-divmigrate-first-significance.png |

## Verification (2026-09-27)
Re-running from this folder reproduced these manuscript figures **byte-identically**:
- `fig-scenarios-combined.png`
- `fig-divmig-correlation-all.png`
- `fig-divmig-lag1-all.png`
- `fig-divmig-lag1-isoband.png`

**Exception:** the manuscript's `media/fig-gravity-w-vs-divmigrate-3x3.png` (1500 × 1087 px, commit `3139ea7`) is Claude's render **trimmed by hand** to remove the caption text under the plot, which the figure legend already covers. Claude's untrimmed original is kept as `media/fig-gravity-w-vs-divmigrate-3x3-claude-original.png`. `gravity_w_vs_divmigrate_3x3.R` now omits the caption, so a fresh render needs no trimming. Its dimensions won't match the hand-trimmed file byte for byte, so don't overwrite that file without checking first.
