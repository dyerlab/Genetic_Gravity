# R/simulation/density_power_hs.R
# Within-population diversity per census of the density/power panel, and whether test power
# collapses onto a within-population information index across the n and L axes.
#   Hs = multilocus Hs (genetic_structure(mode = "Gst"), all 25 demes x 100 individuals)
#   info = n * L_eff * Hs_poly, Hs_poly = Hs rescaled to the polymorphic loci
#          (expected heterozygous genotypes per population summed over loci).
# Output: data/derived/density_power_hs.rda; tables printed (local/density_power_full_results.md, section 8)
suppressPackageStartupMessages({ library(dplyr); library(gstudio); library(parallel) })
source("R/common/gravity_convention.R"); options(width = 200)
load("data/derived/density_power_full.rda")
SCN <- c(Symmetric = 1L, Redistributed = 2L, Obstructed = 3L, `sym-low` = 5L)
f <- "data/derived/density_power_hs.rda"
if (!file.exists(f)) {
  hs <- bind_rows(mclapply(seq_len(nrow(cens)), function(i) { cz <- cens[i, ]
    e <- new.env(); load(sprintf("%s/replicate%d/rep%d-genotypes-scenario%d-%d.rda", GC_SIM_DIR, cz$Replicate, cz$Replicate, SCN[[cz$scenario]], cz$generation), envir = e)
    s <- genetic_structure(e$genotypes, mode = "Gst"); m <- s[s$Locus != "Multilocus", ]
    data.frame(cz, Hs = s$Hs[s$Locus == "Multilocus"], Hs_poly = mean(m$Hs[m$Ht > 0]), n_poly = sum(m$Ht > 0)) },
    mc.cores = max(1L, detectCores() - 1L)))
  save(hs, file = f)
} else load(f)
cat("Within-population Hs (all loci; polymorphic loci only):\n")
print(as.data.frame(hs |> group_by(scenario) |> summarise(lo = quantile(Hs_poly, .1), hi = quantile(Hs_poly, .9), Hs_all = median(Hs), Hs_poly_med = median(Hs_poly))), digits = 3, row.names = FALSE)
d <- res |> filter(!is.na(n_edges), scheme %in% c("spread", "loci"), K %in% c(25, 18, 12, 8)) |>
  left_join(hs |> select(scenario, Replicate, generation, Hs_poly), by = c("scenario", "Replicate", "generation")) |>
  mutate(L_eff = ifelse(scheme == "spread", n_poly, L_eff), info = n * L_eff * Hs_poly,
         axis = ifelse(scheme == "loci", "L varied (n = 30)", "n varied (all loci)"),
         ibin = cut(info, c(0, 15, 30, 60, 120, 250, 1e4)))
cat("\nSource-sink rejection and fragmentation by information bin, K and axis (plateau, both asymmetric scenarios):\n")
print(as.data.frame(d |> filter(scenario %in% c("Redistributed", "Obstructed")) |> group_by(K, ibin, axis) |>
  summarise(G = n(), rejS = mean(p_S < .05, na.rm = TRUE), frag = mean(n_comp > 1), .groups = "drop") |> filter(G >= 40) |>
  arrange(desc(K), ibin, axis)), digits = 2, row.names = FALSE)
cat("\nSymmetric: fragmentation and cor(S,S0) by information bin (K = 25, 12):\n")
print(as.data.frame(d |> filter(scenario == "Symmetric", K %in% c(25, 12)) |> group_by(K, ibin) |>
  summarise(G = n(), frag = mean(n_comp > 1), corSS0 = median(cor_S_S0, na.rm = TRUE), dshare = median(deg_share, na.rm = TRUE), .groups = "drop")), digits = 2, row.names = FALSE)
cat("\nLogistic (plateau): source-sink rejection ~ log K + log n + log L_eff + log Hs_poly\n")
a <- d |> filter(scenario %in% c("Redistributed", "Obstructed"))
print(round(summary(glm(I(p_S < .05) ~ log(K) + log(n) + log(L_eff) + log(Hs_poly), binomial, a))$coef, 3))
cat("\nLogistic (all): fragmented ~ log K + log n + log L_eff + log Hs_poly\n")
print(round(summary(glm(I(n_comp > 1) ~ log(K) + log(n) + log(L_eff) + log(Hs_poly), binomial, d))$coef, 3))
