# lag1_compute.R
# Compute lag-1 series for every replicate with saved matrices (step 1 of Fig_DivMigrateLag1) -> divMig.{rep}.lag1.rda
# Run from the original repo root (~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs).
# Recovered from the 2026-09-24..27 working session scratchpad.

source("R/divmigrate_lag1.R")
f <- list.files("data/divMigrate", "^divMig[.][0-9]+[.]matrices[.]rda$")
done <- sort(as.integer(sub("^divMig[.]([0-9]+)[.]matrices[.]rda$", "\\1", f)))
for (r in done) tryCatch(suppressMessages(divmigrate_lag1(r)), error = function(e) message("skip ", r, ": ", conditionMessage(e)))
cat("lag-1 files for", length(done), "replicates\n")
