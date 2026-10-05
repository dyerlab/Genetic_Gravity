# exploratory/ibgd_gravity_vs_divmigrate.R
#
# Does genetic gravity's closer fit to the true migration matrix carry through
# to isolation by distance? The same directional IBD models are fitted to each
# method's all-pairs directional matrix at every forward census:
#   gravity     D = 1/S - 1, the directed pGD shortest-path distance
#               (S = 1/(1 + D), saved as divmig_S);
#   divMigrate  -log(M), M = relative migration (stat = "Nm"), [i, j] = i -> j
#               (saved as divmig_M).
# Both are "distances" (small = well connected), oriented from row i to column j.
# Landscape: 1-D, unit spacing, x = deme index.
#
# Per census and method, over ordered pairs i != j with a finite value:
#   blind   adj. R^2 of  y ~ |dx|                      (direction ignored)
#   aware   adj. R^2 of  y ~ |dx| + |dx|:reverse       (forward/reverse slopes)
#   gain    aware - blind                              (fit from modelling direction)
#   ratio   forward slope / reverse slope              (< 1: forward is cheaper,
#                                                        the imposed direction)
# Two response scales: "rank" (y ranked within the census; transform-free, the
# primary comparison) and "raw" (the distances as defined above).
#
# Source matrices: divMig.{rep}.matrices.rda (private repo,
# R/divmigrate_matrix_census.R), every census 2004-2999, 50 replicates.
#
# Run from the repository root:  Rscript exploratory/ibgd_gravity_vs_divmigrate.R
# Outputs: data/derived/ibgd_gravity_vs_divmigrate.rda;
#          media/fig-ibgd-gravity-vs-divmigrate.png

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(parallel); library(ggplot2) })

DM_DIR <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
OUT    <- "data/derived/ibgd_gravity_vs_divmigrate.rda"
FIG    <- "media/fig-ibgd-gravity-vs-divmigrate.png"
trt    <- c(iso = "Isotropic", flux = "Redistributed", rate = "Obstructed")
cols   <- c("Genetic gravity (pGD)" = "#2a78d6", "divMigrate (Nm)" = "#eb6834")

x  <- 1:25
X  <- abs(outer(x, x, "-")); REV <- outer(x, x, function(a, b) b < a)   # j upstream of i
OFF <- row(X) != col(X)

#' Directional IBD fit for one K x K distance matrix
fit_one <- function(Y) {
  ok <- OFF & is.finite(Y)
  if (sum(ok) < 50) return(c(n = sum(ok), blind = NA, aware = NA, ratio = NA))
  d <- data.frame(y = Y[ok], X = X[ok], rev = as.numeric(REV[ok]))
  a <- lm(y ~ X + X:rev, d); cf <- coef(a)
  c(n = sum(ok), blind = summary(lm(y ~ X, d))$adj.r.squared, aware = summary(a)$adj.r.squared,
    ratio = unname(cf["X"] / (cf["X"] + cf["X:rev"])))
}
rank_mat <- function(Y) { r <- Y; r[] <- NA; ok <- OFF & is.finite(Y); r[ok] <- rank(Y[ok]); r }

if (file.exists(OUT)) load(OUT) else {
  res <- mclapply(1:50, function(rep) {
    e <- new.env(); load(file.path(DM_DIR, sprintf("divMig.%d.matrices.rda", rep)), envir = e)
    out <- list()
    for (t in names(trt)) {
      M <- e$divmig_M[[t]]; S <- e$divmig_S[[t]]; gens <- as.integer(dimnames(M)[[3]])
      for (k in seq_along(gens)) {
        Dg <- 1 / S[, , k] - 1
        Dm <- suppressWarnings(-log(M[, , k])); Dm[!is.finite(Dm)] <- NA
        for (m in c("gravity", "divmigrate")) {
          Y <- if (m == "gravity") Dg else Dm
          for (sc in c("rank", "raw")) {
            f <- fit_one(if (sc == "rank") rank_mat(Y) else Y)
            out[[length(out) + 1]] <- data.frame(Replicate = rep, Treatment = trt[[t]], generation = gens[k],
                                                 method = m, scale = sc, t(f))
          }
        }
      }
    }
    do.call(rbind, out)
  }, mc.cores = max(1L, detectCores() - 1L))
  fits <- do.call(rbind, res)
  fits$gain <- fits$aware - fits$blind
  dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
  save(fits, file = OUT)
}
fits <- fits |> mutate(Treatment = factor(Treatment, trt),
                       method = factor(method, c("gravity", "divmigrate"), names(cols)),
                       block = cut(generation, seq(2000, 3000, 100), labels = paste0(seq(2000, 2900, 100), "s")))

for (sc in c("rank", "raw")) {
  cat(sprintf("\n==== Response scale: %s ====\n", sc))
  print(fits |> filter(scale == sc) |> group_by(Treatment, block, method) |>
          summarise(pairs = median(n), unfit = mean(is.na(aware)), blind = median(blind, na.rm = TRUE),
                    aware = median(aware, na.rm = TRUE), gain = median(gain, na.rm = TRUE),
                    fwd_cheaper = mean(ratio < 1, na.rm = TRUE), ratio = median(ratio, na.rm = TRUE),
                    .groups = "drop") |>
          pivot_wider(names_from = method, values_from = c(pairs, unfit, blind, aware, gain, ratio, fwd_cheaper),
                      names_glue = "{.value}_{substr(method, 1, 3)}") |>
          select(Treatment, block, starts_with("aware"), starts_with("gain"), starts_with("ratio"),
                 starts_with("fwd"), starts_with("unfit")) |> as.data.frame(), digits = 3)
}

cat("\nPaired within census (rank scale): share where gravity > divMigrate\n")
pw <- fits |> filter(scale == "rank") |>
  select(Replicate, Treatment, generation, block, method, aware, gain) |>
  pivot_wider(names_from = method, values_from = c(aware, gain))
names(pw) <- sub("Genetic gravity \\(pGD\\)", "g", sub("divMigrate \\(Nm\\)", "d", names(pw)))
print(pw |> group_by(Treatment, block) |>
        summarise(aware_g_better = mean(aware_g > aware_d, na.rm = TRUE),
                  gain_g_better = mean(gain_g > gain_d, na.rm = TRUE), .groups = "drop") |>
        as.data.frame(), digits = 3)

## ---- figure: rank scale -------------------------------------------------------------
lab <- c(aware = "Direction-aware IBD fit\nadj. R²", gain = "Fit gained by modelling direction\nΔR²")
fl <- fits |> filter(scale == "rank") |> select(Treatment, method, generation, aware, gain) |>
  pivot_longer(c(aware, gain), names_to = "metric") |>
  mutate(metric = factor(metric, names(lab), lab)) |>
  group_by(Treatment, method, metric, generation) |>
  summarise(med = median(value, na.rm = TRUE), lo = quantile(value, .25, na.rm = TRUE),
            hi = quantile(value, .75, na.rm = TRUE), .groups = "drop")
p <- ggplot(fl, aes(generation, med, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  facet_grid(metric ~ Treatment, scales = "free_y", switch = "y") +
  labs(x = "Generation", y = NULL, colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", strip.placement = "outside", panel.grid.minor = element_blank())
ggsave(FIG, p, width = 10, height = 6, dpi = 150)
cat("\nFigure:", FIG, "\n")
