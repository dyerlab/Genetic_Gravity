# build_timecourse_nm.R
# Builds the Nm time-course table later copied to data/divmigrate_timecourse_timecourse.csv (provenance of the manuscript CSV)
# Run from the original repo root (~/Documents/Research/PopulationGenetics/AsymmetricPopGraphs).
# Recovered from the 2026-09-24..27 working session scratchpad.
S_DIR <- "data/divMigrate/scratch_outputs"; dir.create(S_DIR, showWarnings = FALSE, recursive = TRUE)

suppressMessages({library(gstudio); library(igraph); library(dplyr)})
gens <- c(seq(2004, 2954, by = 50), 2999)
scen <- c(iso = 1L, flux = 2L, rate = 3L); lab <- c(iso = "Isotropic", flux = "Redistributed", rate = "Obstructed")
one_rep <- function(r) {
  e <- new.env(); load(sprintf("data/divMigrate/divMig.%d.matrices.rda", r), e)
  bind_rows(lapply(names(scen), function(s) bind_rows(lapply(gens, function(gen) {
    ge <- new.env(); load(sprintf("data/replicate%d/rep%d-graph-scenario%d-%04d.rda", r, r, scen[[s]], gen), ge)
    g  <- graph_asymmetries(ge$graph); el <- as_edgelist(g, names = TRUE)
    M  <- e$divmig_M[[s]][, , as.character(gen)]
    A  <- M[el] - M[el[, 2:1]]
    data.frame(replicate = r, scenario = lab[[s]], generation = gen,
               gravity_signed = mean(E(g)$delta), gravity_abs = mean(abs(E(g)$delta)),
               divm_signed = mean(A), divm_abs = mean(abs(A)))
  }))))
}
tc <- bind_rows(parallel::mclapply(1:50, one_rep, mc.cores = 13))
write.csv(tc, file.path(S_DIR, "divmigrate_timecourse_nm.csv"), row.names = FALSE)
old <- read.csv("~/Documents/Manuscripts/Genetic Gravity/data/divmigrate_timecourse_timecourse.csv")
m <- merge(tc, old, by = c("replicate", "scenario", "generation"))
cat("rows:", nrow(tc), " matched to D file:", nrow(m), "\n")
cat("gravity_signed max|diff| vs cached:", max(abs(m$gravity_signed.x - m$gravity_signed.y)), "\n")
cat("non-finite Nm divm_signed:", sum(!is.finite(tc$divm_signed)), "  D-based:", sum(!is.finite(old$divm_signed)), "\n")
cat("Nm divm_abs range:", signif(range(tc$divm_abs, na.rm = TRUE), 3), "\n")
cat("cor(D, Nm) divm_signed:", round(cor(m$divm_signed.x, m$divm_signed.y, use = "complete.obs", method = "spearman"), 2), "\n")
