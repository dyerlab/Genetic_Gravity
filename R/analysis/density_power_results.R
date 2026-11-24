# R/analysis/density_power_results.R
#
# Tables from data/derived/density_power_full.rda (R/simulation/density_power_full.R).
# alpha = 0.05 throughout; intervals are 95% lineage-bootstrap intervals (B = 1000, lineages
# resampled within scenario).
#
# Delta R^2 cutoff, cross-fitted by lineage: lineages split into folds A (1-25) and B (26-50).
# For each (K, n) cell the cutoff c is the 95th percentile of Delta R^2 over the Symmetric
# graphs of one fold, and every graph (any scenario) of lineage r is judged against the cutoff
# from the fold NOT containing r. Asymmetric lineages branch from the same burn-in state as the
# Symmetric lineage of the same number, so this keeps each lineage's own history out of its
# cutoff. A K-only cutoff (pooled over n) is computed the same way, to test dependence on n.
# The sym-low false-positive rate is read at the Symmetric cutoffs.
#
# Run: Rscript R/analysis/density_power_results.R   (env IN to override the input)
# Output: printed tables (paste into local/density_power_full_results.md) and
#         data/derived/density_power_tables.rda

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(gstudio) })
options(width = 220, dplyr.summarise.inform = FALSE)
IN <- Sys.getenv("IN", "data/derived/density_power_full.rda"); load(IN)
set.seed(20261124L)
ok <- res |> filter(!is.na(n_edges)) |> mutate(fold = ifelse(Replicate <= 25, "A", "B"),
  sample = factor(scenario, c("Symmetric", "sym-low", "Redistributed", "Obstructed")))
cat(sprintf("%d graphs, %d censuses, %d errors; NA p_S %d, NA dR2 %d\n\n", nrow(ok), nrow(cens),
            if ("err" %in% names(res)) sum(!is.na(res$err)) else 0L, sum(is.na(ok$p_S)), sum(is.na(ok$dR2))))

#' Lineage-bootstrap mean of a logical vector (NA dropped), as "p (lo-hi)" pieces
boot_rate <- function(y, lin, B = 1000L) {
  m <- !is.na(y); y <- y[m]; lin <- lin[m]; if (!length(y)) return(c(est = NA, lo = NA, hi = NA))
  s <- tapply(y, lin, sum); n <- tapply(y, lin, length); L <- length(s)
  bs <- replicate(B, { i <- sample.int(L, L, TRUE); sum(s[i]) / sum(n[i]) })
  c(est = mean(y), lo = unname(quantile(bs, .025)), hi = unname(quantile(bs, .975)))
}
fmt <- function(v) sprintf("%.3f (%.3f-%.3f)", v[1], v[2], v[3])
rate_tab <- function(d, yexpr, by) d |> group_by(across(all_of(by))) |>
  group_modify(~ { v <- boot_rate(eval(yexpr, .x), .x$Replicate); tibble(G = sum(!is.na(eval(yexpr, .x))), est = v[1], lo = v[2], hi = v[3]) }) |> ungroup()

sp <- ok |> filter(scheme == "spread")

## ---- 1. Delta R^2 cutoffs, cross-fitted -----------------------------------------------------
cut_cell <- sp |> filter(scenario == "Symmetric") |> group_by(fold, K, n) |>
  summarise(c_cell = quantile(dR2, .95, na.rm = TRUE), .groups = "drop")
cut_K <- sp |> filter(scenario == "Symmetric") |> group_by(fold, K) |>
  summarise(c_K = quantile(dR2, .95, na.rm = TRUE), .groups = "drop")
other <- function(f) ifelse(f == "A", "B", "A")
apply_cut <- function(d) d |> mutate(cf = other(fold)) |>
  left_join(cut_cell |> rename(cf = fold), by = c("cf", "K", "n")) |>
  left_join(cut_K |> rename(cf = fold), by = c("cf", "K")) |>
  mutate(det0 = dR2 > 0, det_cell = dR2 > pmax(0, c_cell), det_K = dR2 > pmax(0, c_K))
spc <- apply_cut(sp)
cat("== 1a. Symmetric Delta R^2 95th percentile by (K, n), mean of the two folds ==\n")
print(as.data.frame(cut_cell |> group_by(K, n) |> summarise(cutoff = mean(c_cell), fold_diff = diff(range(c_cell))) |>
  pivot_wider(id_cols = K, names_from = n, values_from = cutoff, names_prefix = "n=")), digits = 3, row.names = FALSE)
cat("\n== 1b. Delta R^2 rejection: rule > 0, held-out cell cutoff, held-out K-only cutoff ==\n")
dr <- bind_rows(lapply(c(det0 = "det0", det_cell = "det_cell", det_K = "det_K"), function(v)
  rate_tab(spc, as.name(v), c("sample", "K", "n")) |> mutate(rule = v)))
print(as.data.frame(dr |> mutate(v = sprintf("%.2f", est)) |> select(sample, K, n, rule, v) |>
  pivot_wider(names_from = c(rule), values_from = v) |> arrange(sample, desc(K), desc(n))), row.names = FALSE)

## ---- 2. Source-sink test: whole graph vs largest component ------------------------------------
cat("\n== 2. Source-sink test rejection (p < 0.05), whole graph; [largest component where fragmented] ==\n")
ss <- rate_tab(sp, quote(p_S < .05), c("sample", "K", "n"))
print(as.data.frame(ss |> mutate(v = sprintf("%.2f", est)) |> select(sample, K, n, v) |>
  pivot_wider(names_from = n, values_from = v, names_prefix = "n=") |> arrange(sample, desc(K))), row.names = FALSE)
fr <- sp |> filter(n_comp > 1)
cat(sprintf("\nFragmented graphs: %d of %d (%.0f%%)\n", nrow(fr), nrow(sp), 100 * nrow(fr) / nrow(sp)))
print(as.data.frame(fr |> group_by(sample) |> summarise(G = n(), lc_frac = median(lc_size / K),
  rejS_whole = mean(p_S < .05, na.rm = TRUE), rejS_lc = mean(p_S_lc < .05, na.rm = TRUE),
  dR2pos_whole = mean(dR2 > 0, na.rm = TRUE), dR2pos_lc = mean(dR2_lc > 0, na.rm = TRUE))), digits = 2, row.names = FALSE)
print(as.data.frame(sp |> group_by(sample, n) |> summarise(frag = mean(n_comp > 1), comp = mean(n_comp)) |>
  pivot_wider(id_cols = sample, names_from = n, values_from = frag, names_prefix = "frag n=")), digits = 2, row.names = FALSE)

## ---- 3. Descriptions: sign of r(S, x) and slope ratio -----------------------------------------
cat("\n== 3. Sources upstream (r(S,x) < 0) and slope ratio < 1, all graphs, by K (pooled n >= 30 | n = 10) ==\n")
desc <- sp |> mutate(nband = ifelse(n >= 30, "n>=30", "n<=20")) |> group_by(sample, K, nband) |>
  summarise(up = mean(r_S < 0, na.rm = TRUE), ratio_lt1 = mean(ratio < 1, na.rm = TRUE), .groups = "drop")
print(as.data.frame(desc |> pivot_wider(id_cols = c(sample, K), names_from = nband, values_from = c(up, ratio_lt1)) |>
  arrange(sample, desc(K))), digits = 2, row.names = FALSE)
lin_up <- sp |> filter(n >= 30) |> group_by(sample, K, Replicate) |> summarise(mr = mean(r_S, na.rm = TRUE), .groups = "drop") |>
  group_by(sample, K) |> summarise(lineages_up = mean(mr < 0, na.rm = TRUE), .groups = "drop")
cat("\nLineage-averaged r(S,x) < 0 (n >= 30):\n")
print(as.data.frame(lin_up |> pivot_wider(names_from = sample, values_from = lineages_up) |> arrange(desc(K))), digits = 2, row.names = FALSE)

## ---- 4. Contiguous vs spread ----------------------------------------------------------------------
cat("\n== 4. Contiguous block vs spread, same K and n ==\n")
cs <- apply_cut(ok |> filter(scheme %in% c("spread", "contiguous"), K %in% unique(res$K[res$scheme == "contiguous"]),
                             n %in% unique(res$n[res$scheme == "contiguous"]), scenario != "sym-low"))
print(as.data.frame(cs |> group_by(sample, K, n, scheme) |> summarise(G = n(), span = median(span), deg = median(mean_deg),
  rejS = mean(p_S < .05, na.rm = TRUE), up = mean(r_S < 0, na.rm = TRUE), dR2_0 = mean(det0, na.rm = TRUE),
  dR2_cut = mean(det_cell, na.rm = TRUE), .groups = "drop") |> arrange(sample, desc(K), desc(n), scheme)), digits = 2, row.names = FALSE)

## ---- 5. Loci ----------------------------------------------------------------------------------------
cat("\n== 5. Loci at n = 30 (L = 20 from the spread scheme; L = 10, 8 polymorphic loci) ==\n")
lo <- apply_cut(ok |> filter(n == 30, K %in% unique(res$K[res$scheme == "loci"]), scenario != "sym-low",
                             scheme %in% c("spread", "loci")) |> mutate(L = ifelse(scheme == "spread", 20L, L)))
print(as.data.frame(lo |> group_by(sample, K, L) |> summarise(G = n(), deg = median(mean_deg), dshare = median(deg_share, na.rm = TRUE),
  frag = mean(n_comp > 1), rejS = mean(p_S < .05, na.rm = TRUE), dR2_0 = mean(det0, na.rm = TRUE),
  dR2_cut = mean(det_cell, na.rm = TRUE), .groups = "drop") |> arrange(sample, desc(K), desc(L))), digits = 2, row.names = FALSE)

## ---- 6. Degree share, mean degree and gamma = 1/2 neutrality ------------------------------------
cat("\n== 6. Median degree share, mean degree, cor(S,S0) (Symmetric) by (K, n) ==\n")
print(as.data.frame(sp |> group_by(K, n) |> summarise(dshare = median(deg_share, na.rm = TRUE), deg = median(mean_deg),
  corSS0_sym = median(cor_S_S0[scenario == "Symmetric"], na.rm = TRUE), .groups = "drop") |> arrange(desc(K), desc(n))),
  digits = 2, row.names = FALSE)
asy <- sp |> filter(scenario %in% c("Redistributed", "Obstructed"))
cat("\nLogistic, plateau: source-sink rejection ~ degree share + log K + log n\n")
print(round(summary(glm(I(p_S < .05) ~ deg_share + log(K) + log(n), binomial, asy))$coef, 3))
cat("\nLogistic, plateau: Delta R^2 > held-out cell cutoff ~ degree share + log K + log n\n")
print(round(summary(glm(det_cell ~ deg_share + log(K) + log(n), binomial, spc |> filter(scenario %in% c("Redistributed", "Obstructed"))))$coef, 3))

## ---- 7. Phi and the Falconer sampling optimum -------------------------------------------------------
cat("\n== 7. Multilocus Fst per census and the Falconer optimum ==\n")
cf <- res |> filter(!is.na(fst)) |> distinct(scenario, Replicate, generation, fst, n_poly)
print(as.data.frame(cf |> group_by(scenario) |> summarise(censuses = n(), fst_med = median(fst), fst_lo = quantile(fst, .1),
  fst_hi = quantile(fst, .9), poly = median(n_poly))), digits = 3, row.names = FALSE)
phi <- median(cf$fst[cf$scenario %in% c("Redistributed", "Obstructed")])
fal <- bind_rows(lapply(c(250, 500, 1000, 2500), function(N) {
  o <- optimal_sampling(N, phi); best <- o[which.min(o$var_phi), ]
  data.frame(N = N, phi = phi, strata_at_min = best$Stratum, per_stratum = N / best$Stratum, var_min = best$var_phi) }))
cat("\nFalconer optimum at the plateau median Phi (optimal_sampling(N, phi)):\n"); print(fal, digits = 3, row.names = FALSE)

save(cut_cell, cut_K, dr, ss, desc, lin_up, cf, fal, file = "data/derived/density_power_tables.rda")
cat("\nWrote data/derived/density_power_tables.rda\n")
