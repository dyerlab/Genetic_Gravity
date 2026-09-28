# exploratory/node_potential_reanalysis_followup.R
#
# Two checks prompted by exploratory/node_potential_reanalysis.R (results file, Tasks 3a and the
# symmetric-gradient Delta R^2):
#   corridor  On the 24 true stepping-stone links, is the estimate larger in the forward direction?
#             Gravity: w_{i|i+1} (W[i+1, i]) vs w_{i+1|i} (W[i, i+1]), i.e. sign of Delta_{i->i+1},
#             links retained in the graph only. divMigrate: M[i, i+1] vs M[i+1, i] (row = source),
#             censuses where Nm is defined. Direction-specific, unlike the all-pairs rho.
#   symdir    Delta R^2 > 0 on the symmetric gradient: slope ratio b_fwd/b_rev among detections,
#             and how detections spread over lineages.
# Output: added to data/derived/node_potential_reanalysis.rda as rean$corridor and rean$symdir.

suppressPackageStartupMessages({ library(dplyr); library(parallel) })
e <- new.env(); load("data/derived/node_potential_reanalysis.rda", envir = e); rean <- e$rean
DM <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
POPS <- sprintf("Pop%02d", 1:25); LAB <- c(iso = "Isotropic", flux = "Redistributed", rate = "Obstructed")
i <- 1:24
rows <- mclapply(1:50, function(r) {
  ew <- new.env(); load(file.path(DM, sprintf("divMig.%d.w.rda", r)), envir = ew)
  em <- new.env(); load(file.path(DM, sprintf("divMig.%d.matrices.rda", r)), envir = em)
  out <- list()
  for (lab in names(LAB)) {
    Wa <- ew$gravity_W[[lab]]; Ma <- em$divmig_M[[lab]]; gens <- dimnames(Wa)[[3]]
    for (j in seq_along(gens)) {
      W <- Wa[POPS, POPS, j]; M <- Ma[POPS, POPS, gens[j]]
      fw <- W[cbind(i + 1, i)]; rv <- W[cbind(i, i + 1)]; ret <- fw > 0 & rv > 0
      off <- row(M) != col(M); ok <- all(is.finite(M[off])) && all(M[off] >= 0)
      mf <- M[cbind(i, i + 1)]; mr <- M[cbind(i + 1, i)]
      out[[length(out) + 1]] <- data.frame(Replicate = r, Treatment = LAB[[lab]], Generation = as.integer(gens[j]),
        links = sum(ret), w_fwd_gt = mean(fw[ret] > rv[ret]),
        nm_ok = ok, nm_fwd_gt = if (ok) mean(mf > mr) else NA)
    }
  }
  do.call(rbind, out)
}, mc.cores = max(1L, detectCores() - 1L))
rean$corridor <- list(d = do.call(rbind, rows))
tr <- rean$ibgd_sym$traj |> filter(status == "ok") |> mutate(det = r2_pgd_dir > r2_cgd)
rean$symdir <- list(
  by_scen = tr |> group_by(Treatment) |>
    summarise(n = n(), fpr = mean(det), det_ratio_lt1 = mean(rel_fwd[det] < 1), det_med_ratio = median(rel_fwd[det]),
              all_med_ratio = median(rel_fwd), .groups = "drop"),
  by_rep = tr |> group_by(Treatment, Replicate) |> summarise(fpr = mean(det), .groups = "drop"),
  by_block = tr |> mutate(block = cut(generation, c(2000, 2250, 2500, 2750, 3000))) |>
    group_by(Treatment, block) |> summarise(fpr = mean(det), .groups = "drop"))
save(rean, file = "data/derived/node_potential_reanalysis.rda")

if (sys.nframe() == 0L) {
  options(width = 150)
  cat("=== Corridor links: fraction where the forward-direction estimate is larger ===\n")
  print(rean$corridor$d |> mutate(win = cut(Generation, c(2000, 2250, 2500, 2750, 3000))) |>
          group_by(Treatment, win) |>
          summarise(links = median(links), w_fwd = mean(w_fwd_gt, na.rm = TRUE),
                    nm_fwd = mean(nm_fwd_gt, na.rm = TRUE), nm_censuses = sum(nm_ok), .groups = "drop") |>
          as.data.frame(), digits = 3)
  cat("\n=== Symmetric gradient, Delta R^2 > 0 ===\n"); print(as.data.frame(rean$symdir$by_scen), digits = 3)
  print(as.data.frame(rean$symdir$by_block), digits = 3)
  cat("\nPer-lineage FPR quantiles:\n")
  print(rean$symdir$by_rep |> group_by(Treatment) |>
          summarise(med = median(fpr), q90 = quantile(fpr, .9), max = max(fpr), zero = sum(fpr == 0)) |> as.data.frame(),
        digits = 3)
  # isotropic comparison from the stored trajectory
  ef <- new.env(); load("data/derived/ibgd_trajectory_fit.rda", envir = ef)
  iso <- ef$traj |> filter(Treatment %in% c("Isotropic", "Burn-in"), status == "ok") |> mutate(det = r2_pgd_dir > r2_cgd)
  cat("\nIsotropic / burn-in detections: fraction with ratio < 1\n")
  print(iso |> group_by(Treatment) |> summarise(fpr = mean(det), det_ratio_lt1 = mean(rel_fwd[det] < 1),
                                                med_cgd = median(r2_cgd)) |> as.data.frame(), digits = 3)
  cat("\nMedian R^2_cGD on the gradient:\n")
  print(rean$ibgd_sym$traj |> group_by(Treatment) |> summarise(r2_cgd = median(r2_cgd, na.rm = TRUE),
        r2_dir = median(r2_pgd_dir, na.rm = TRUE), edges = median(edges, na.rm = TRUE)) |> as.data.frame(), digits = 3)
}
