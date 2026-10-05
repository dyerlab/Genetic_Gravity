# exploratory/ibd_naive.R
#
# Naive isolation by distance: the analyst knows only where the populations
# were sampled. Geographic distance is the Euclidean distance between sampling
# locations (1-D landscape, unit spacing: x = deme index). Each method's
# genetic distance matrix is correlated with it, with no direction term:
#   cGD         undirected shortest-path distance over Population Graph edge
#               weights (symmetric; unordered pairs);
#   pGD         directed shortest-path distance over the p_ij arc costs,
#               D = 1/S - 1 (asymmetric; ordered pairs i != j);
#   divMigrate  -log(M), M = relative migration, stat = "Nm" (asymmetric;
#               ordered pairs). Non-finite entries are dropped.
# Statistics per census and method: Pearson r (the Mantel r) and Spearman rho
# between genetic and geographic distance, over pairs with finite values.
#
# Sources: forward-phase graphs (rep{r}-graph-scenario{s}-{gen}.rda) and the
# saved divMigrate / pGD matrices (divMig.{rep}.matrices.rda), every census
# 2004-2999, 50 replicates, three treatments.
#
# Run from the repository root:  Rscript exploratory/ibd_naive.R
# Outputs: data/derived/ibd_naive.rda; media/fig-ibd-naive.png

suppressPackageStartupMessages({ library(igraph); library(dplyr); library(tidyr); library(parallel)
  library(ggplot2) })

SIM_DIR <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"
OUT <- "data/derived/ibd_naive.rda"
FIG <- "media/fig-ibd-naive.png"
trt <- c(iso = "Isotropic", flux = "Redistributed", rate = "Obstructed")
scen <- c(iso = 1L, flux = 2L, rate = 3L)
pops <- sprintf("Pop%02d", 1:25)
GEO <- abs(outer(1:25, 1:25, "-")); dimnames(GEO) <- list(pops, pops)
UP <- upper.tri(GEO); OFF <- row(GEO) != col(GEO)
methods <- c(cgd = "cGD (Population Graph)", pgd = "pGD (genetic gravity)", dm = "divMigrate (Nm)")
cols <- c("cGD (Population Graph)" = "#2a78d6", "pGD (genetic gravity)" = "#eb6834",
          "divMigrate (Nm)" = "#1baf7a")

cor2 <- function(y, mask) { ok <- mask & is.finite(y)
  c(pairs = sum(ok), pearson = cor(y[ok], GEO[ok]), spearman = cor(y[ok], GEO[ok], method = "spearman")) }

if (file.exists(OUT)) load(OUT) else {
  res <- mclapply(1:50, function(rep) {
    e <- new.env(); load(file.path(SIM_DIR, "divMigrate", sprintf("divMig.%d.matrices.rda", rep)), envir = e)
    out <- list()
    for (t in names(trt)) {
      gens <- as.integer(dimnames(e$divmig_M[[t]])[[3]])
      for (k in seq_along(gens)) {
        f <- sprintf("%s/replicate%d/rep%d-graph-scenario%d-%d.rda", SIM_DIR, rep, rep, scen[[t]], gens[k])
        C <- matrix(NA_real_, 25, 25, dimnames = list(pops, pops))
        if (file.exists(f)) {
          g <- new.env(); load(f, envir = g)
          d <- distances(g$graph, weights = E(g$graph)$weight)
          C[rownames(d), colnames(d)] <- d
        }
        P  <- 1 / e$divmig_S[[t]][pops, pops, k] - 1
        DM <- suppressWarnings(-log(e$divmig_M[[t]][pops, pops, k]))
        s <- rbind(cgd = cor2(C, UP), pgd = cor2(P, OFF), dm = cor2(DM, OFF))
        out[[length(out) + 1]] <- data.frame(Replicate = rep, Treatment = trt[[t]], generation = gens[k],
                                             method = rownames(s), s, row.names = NULL)
      }
    }
    do.call(rbind, out)
  }, mc.cores = max(1L, detectCores() - 1L))
  ibd <- do.call(rbind, res)
  dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
  save(ibd, file = OUT)
}
ibd <- ibd |> mutate(Treatment = factor(Treatment, trt), method = factor(methods[method], methods),
                     block = cut(generation, seq(2000, 3000, 100), labels = paste0(seq(2000, 2900, 100), "s")))

cat("Median Pearson r / Spearman rho with geographic distance, by 100-generation block:\n")
print(ibd |> group_by(Treatment, block, method) |>
        summarise(r = median(pearson, na.rm = TRUE), rho = median(spearman, na.rm = TRUE), .groups = "drop") |>
        mutate(v = sprintf("%.2f / %.2f", r, rho)) |> select(-r, -rho) |>
        pivot_wider(names_from = method, values_from = v) |> as.data.frame())

cat("\nPaired within census: share where each graph measure beats divMigrate (Pearson | Spearman)\n")
pw <- ibd |> select(Replicate, Treatment, generation, block, method, pearson, spearman) |>
  pivot_wider(names_from = method, values_from = c(pearson, spearman))
nm <- function(stat, m) paste0(stat, "_", methods[[m]])
print(pw |> group_by(Treatment, block) |>
        summarise(cgd_P = mean(.data[[nm("pearson", "cgd")]] > .data[[nm("pearson", "dm")]], na.rm = TRUE),
                  cgd_S = mean(.data[[nm("spearman", "cgd")]] > .data[[nm("spearman", "dm")]], na.rm = TRUE),
                  pgd_P = mean(.data[[nm("pearson", "pgd")]] > .data[[nm("pearson", "dm")]], na.rm = TRUE),
                  pgd_S = mean(.data[[nm("spearman", "pgd")]] > .data[[nm("spearman", "dm")]], na.rm = TRUE),
                  .groups = "drop") |> as.data.frame(), digits = 2)

## ---- figure: Pearson (Mantel) r over time -------------------------------------------
fl <- ibd |> group_by(Treatment, method, generation) |>
  summarise(med = median(pearson, na.rm = TRUE), lo = quantile(pearson, .25, na.rm = TRUE),
            hi = quantile(pearson, .75, na.rm = TRUE), .groups = "drop")
p <- ggplot(fl, aes(generation, med, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  facet_wrap(~Treatment, nrow = 1) +
  labs(x = "Generation", y = "IBD, r(genetic, geographic distance)", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(FIG, p, width = 10, height = 4, dpi = 150)
cat("\nFigure:", FIG, "\n")
