# export_truth_compare.R
# Exports per-census truth metrics (pGD similarity + Nm) -> data/divmigrate_truth_compare.csv (input of R/divmigrate_truth.R)
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.

dm <- "~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
fs <- list.files(dm, "^divMig[.][0-9]+[.]truth_sim[.]rda$", full.names = TRUE)
d <- do.call(rbind, lapply(fs, function(f) { e <- new.env(); load(f, e); e$divmig_truth }))
d$scenario <- c(iso = "Isotropic", flux = "Redistributed", rate = "Obstructed")[as.character(d$Treatment)]
d$method   <- c(pGD = "gravity", Nm = "divMigrate")[as.character(d$Method)]
out <- d[order(d$scenario, d$Replicate, d$Generation, d$method),
         c("Replicate", "scenario", "Generation", "method", "Excluded", "Exclusion",
           "rho_all", "rmse_all", "mae_all", "rho_corr", "rmse_corr", "mae_corr", "auc_corr")]
names(out)[1:3] <- c("replicate", "scenario", "generation")
write.csv(out, "data/divmigrate_truth_compare.csv", row.names = FALSE)
cat(nrow(out), "rows;", length(unique(out$replicate)), "replicates\n")
