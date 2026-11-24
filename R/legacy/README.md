# Legacy producers

These scripts wrote cached inputs in `data/derived/` that the current code still reads. They came before the analysis moved to stock gstudio functions (gstudio 1.15.0, October 2026). Several rely on project helpers that have since been retired: `gc_census()`, `node_fields()` and `ibd_one()` (the last in `ibd_cgd_vs_pgd.R`). Several also read divMigrate output held in the simulation repository. They are kept as the record of how each file was made. `R/common/gravity_convention.R` checks that the current gstudio functions reproduce the archived values.

| Script | Cached output | Used by |
|---|---|---|
| `node_potential_reanalysis.R`, `node_potential_reanalysis_followup.R` | `node_potential_reanalysis.rda` (divMigrate *Nm* and its comparison with the true migration matrix; γ = 1 test results) | `analysis/half_bandwidth_results.R`, `analysis/phase_boundaries.R`, `analysis/phase_core_rerun.R`, `figures/fig_gravity_vs_divmigrate.R`, `figures/fig_source_sink_score.R` |
| `node_potential_test.R` (with `node_potential_helpers.R`) | `node_potential_test.rda` (first-round γ = 1 source–sink test) | `analysis/half_bandwidth_results.R` |
| `node_potential_followup.R` | `node_potential_followup.rda`, `node_potential_nodes.rds` (archived edge and node values) | the check in `common/gravity_convention.R` |
| `node_potential_step0.R` | Table S2 (identities), printed | — |
| `phase_partitioned_results.R`, `ss_bandwidth_by_phase.R` | `phase_partitioned_results.rda`, `ss_bandwidth_by_phase.rda` (γ = 1 bandwidth tables) | `analysis/half_bandwidth_results.R` |
| `ibgd_trajectory.R` (with `ibd_cgd_vs_pgd.R`) | `ibgd_trajectory_fit.rda` (γ = 1 IBGD fits) | `analysis/half_bandwidth_results.R`, `analysis/lineage_phase_windows.R`, `analysis/phase_boundaries.R`, `figures/fig_ibgd_fit.R` |
