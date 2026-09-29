# exploratory/he_gradient_boundary.R
#
# A population-genetic (gravity-free) marker for the end of reorganization.
#
# Why: the current boundary (R_gain, t90 of the within-pGD directional gain) is
# itself a gravity statistic, so "the plateau is where gravity inference is
# strongest" is partly true by construction. The falling-apart boundary (the
# informative horizon, polymorphism) is already gravity-free; this looks for a
# matching reorganization marker read directly from allele frequencies.
#
# Marker: the within-deme expected heterozygosity gradient along the chain.
# Under directional flow, upstream demes receive less immigration and lose He
# faster, so He rises with deme index. Per census: He_d = mean over loci of
# 2 p (1 - p) in deme d (biallelic), and the gradient is the OLS slope of He_d on
# deme index over the interior demes 3-23 (terminal demes carry their own
# boundary deficit; Supplementary Materials). Also stored: the slope over all
# demes, and the Spearman correlation.
# Per lineage, as for R_gain: penalized-spline smooth (k = 10) of the slope over
# 2004-2999, rise from the lineage's own burn-in mean (1904-1999) to its
# smoothed maximum, R_he = first generation the smooth completes 90% of the rise.
# Specificity: a lineage counts as reorganized only if its rise exceeds the 95th
# percentile of isotropic rises (the same rule used for R_gain).
#
# Run from the repository root:  Rscript exploratory/he_gradient_boundary.R
# Outputs: data/derived/he_gradient.rds (per census: slope_int, slope_all, rho),
#          data/derived/he_gradient_boundary.rda (per lineage: R_he and comparison)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(mgcv); library(parallel) })
SIM <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"
CACHE <- "data/derived/he_gradient.rds"
SC <- c(Isotropic = 1L, Redistributed = 2L, Obstructed = 3L)
grid <- seq(2004, 2999, 5)

he_census <- function(f) {
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e); g <- e$genotypes
  L <- grep("^L[0-9]", names(g))
  M <- vapply(L, function(j) { a <- as.character(g[[j]]); (substr(a, 1, 2) == "01") + (substr(a, 4, 5) == "01") },
              numeric(nrow(g)))
  n <- as.vector(table(g$Population)[sort(unique(g$Population))])
  p <- rowsum(M, g$Population)[sort(unique(g$Population)), , drop = FALSE] / (2 * n)
  he <- rowMeans(2 * p * (1 - p))
  x <- as.integer(sub("\\D+", "", names(he))); int <- x >= 3 & x <= 23
  c(slope_int = unname(coef(lm(he[int] ~ x[int]))[2]), slope_all = unname(coef(lm(he ~ x))[2]),
    rho = suppressWarnings(cor(he, x, method = "spearman")), he_mean = mean(he))
}

if (file.exists(CACHE)) hg <- readRDS(CACHE) else {
  jobs <- rbind(expand.grid(Replicate = 1:50, scenario = "Burn-in", generation = seq(1904L, 1999L, 5L), stringsAsFactors = FALSE),
                expand.grid(Replicate = 1:50, scenario = names(SC), generation = seq(2004L, 2999L, 5L), stringsAsFactors = FALSE))
  jobs$file <- with(jobs, ifelse(scenario == "Burn-in",
    sprintf("%s/replicate%d/rep%d-genotypes-%d.rda", SIM, Replicate, Replicate, generation),
    sprintf("%s/replicate%d/rep%d-genotypes-scenario%d-%d.rda", SIM, Replicate, Replicate, SC[scenario], generation)))
  res <- mclapply(seq_len(nrow(jobs)), function(i) { v <- he_census(jobs$file[i]); if (is.null(v)) NULL else cbind(jobs[i, 1:3], t(v)) },
                  mc.cores = max(1L, detectCores() - 1L))
  hg <- do.call(rbind, res); saveRDS(hg, CACHE)
}

sm <- function(g, y) as.numeric(predict(gam(y ~ s(g, k = 10), data = data.frame(g = g, y = y), method = "REML"), data.frame(g = grid)))
base <- hg |> filter(scenario == "Burn-in") |> group_by(Replicate) |> summarise(b = mean(slope_int), .groups = "drop")
out <- list()
for (s in names(SC)) for (r in 1:50) {
  x <- hg |> filter(scenario == s, Replicate == r) |> arrange(generation)
  f <- sm(x$generation, x$slope_int); b <- base$b[base$Replicate == r]; amp <- max(f) - b
  i90 <- which(f >= b + 0.9 * amp)[1]
  out[[length(out) + 1]] <- data.frame(Scenario = s, Replicate = r, he_amp = amp,
                                       R_he = if (is.na(i90)) NA else grid[i90], he_peak_gen = grid[which.max(f)])
}
hb <- do.call(rbind, out)
thr <- quantile(hb$he_amp[hb$Scenario == "Isotropic"], 0.95)
hb$reorg_he <- hb$he_amp > thr
lp <- read.csv("data/derived/lineage_phases.csv", stringsAsFactors = FALSE)
cmp <- left_join(hb, lp, by = c("Scenario", "Replicate"))
save(hb, thr, cmp, file = "data/derived/he_gradient_boundary.rda")

if (sys.nframe() == 0L) {
  options(width = 160)
  cat("Censuses:", nrow(hg), "\n\nHe gradient (interior OLS slope), median by scenario and 100-generation block:\n")
  print(hg |> filter(scenario != "Burn-in") |> mutate(block = cut(generation, seq(2000, 3000, 200), dig.lab = 4)) |>
          group_by(scenario, block) |> summarise(slope = signif(median(slope_int), 3), rho = round(median(rho), 2), .groups = "drop") |>
          pivot_wider(names_from = scenario, values_from = c(slope, rho)) |> as.data.frame())
  cat(sprintf("\nBurn-in slope: median %.2e; isotropic 95th-pct rise threshold %.2e\n", median(base$b), thr))
  cat("\nLineages reorganized (rise > isotropic 95th pct) and R_he median (IQR):\n")
  print(hb |> group_by(Scenario) |> summarise(reorg = sum(reorg_he), q25 = quantile(R_he[reorg_he], .25, na.rm = TRUE),
          q75 = quantile(R_he[reorg_he], .75, na.rm = TRUE), R_he_med = median(R_he[reorg_he], na.rm = TRUE),
          peak = median(he_peak_gen[reorg_he])) |> as.data.frame())
  cat("\nAgreement with R_gain (asymmetric lineages reorganized under both):\n")
  print(cmp |> filter(Scenario != "Isotropic") |> group_by(Scenario) |>
          summarise(both = sum(reorg_he & reorganized), he_only = sum(reorg_he & !reorganized), gain_only = sum(!reorg_he & reorganized),
                    med_diff = median((R_he - R_gen)[reorg_he & reorganized], na.rm = TRUE),
                    rho = cor(R_he, R_gen, use = "complete.obs", method = "spearman"),
                    R_he_after_F = sum(R_he >= F_gen, na.rm = TRUE)) |> as.data.frame(), digits = 3)
}
