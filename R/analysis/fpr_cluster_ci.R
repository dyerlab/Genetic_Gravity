# R/analysis/fpr_cluster_ci.R
#
# Lineage-clustered 95% intervals for the false-positive rates. Censuses within a lineage
# are autocorrelated (shared history), so the unit of replication is the lineage: the
# interval is a cluster bootstrap that resamples the 50 lineages with replacement
# (B = 10,000, percentile interval; the rate is pooled over the resampled lineages'
# censuses). The naive binomial (Clopper-Pearson) interval over censuses and the design
# effect (cluster-bootstrap variance / binomial variance) are given for comparison.
#
# Rates: source-sink test (graph spectral null) at alpha = 0.01 / 0.05 / 0.10, 50-generation
# grid (burn-in every 5 generations); Delta R^2 > 0 and the post hoc rule (Delta R^2 > 0 and
# slope ratio < 1), 5-generation grid. gamma = 1/2 and gamma = 1 (original lineages,
# R/simulation/half_bandwidth_compute.R) and the 50 fresh lineages (R/simulation/fresh_replicates.R).
#
# Run from the repository root:  Rscript R/analysis/fpr_cluster_ci.R
# Output: local/fpr_cluster_ci.md, data/derived/fpr_cluster_ci.rda

dir.create("output", showWarnings = FALSE)   # run reports (not tracked)
suppressPackageStartupMessages(library(dplyr))
set.seed(20261002); NB <- 10000L
ld <- function(f) { e <- new.env(); load(f, envir = e); e }
SYMS <- c("Burn-in", "Symmetric", "sym-mid", "sym-low", "sym-verylow")

clus <- function(lin, x) {                         # x logical per census, lin lineage id
  s <- tapply(x, lin, sum); n <- tapply(x, lin, length); k <- length(s)
  bt <- replicate(NB, { i <- sample.int(k, k, replace = TRUE); sum(s[i]) / sum(n[i]) })
  p <- sum(s) / sum(n); bn <- binom.test(sum(s), sum(n))$conf.int
  data.frame(censuses = sum(n), lineages = k, rate = p, lo = unname(quantile(bt, .025)), hi = unname(quantile(bt, .975)),
             bin_lo = bn[1], bin_hi = bn[2], deff = var(bt) / (p * (1 - p) / sum(n)),
             lin_with_any = sum(s > 0), lin_max = max(s / n))
}
rows <- list()
add <- function(set, test, sample, d) rows[[length(rows) + 1]] <<- cbind(data.frame(set = set, test = test, sample = sample), d)

for (g in c("05", "1")) {
  lab <- if (g == "05") "Original lineages, γ = ½" else "Original lineages, γ = 1"
  T <- ld(sprintf("data/derived/half_bandwidth_tests_g%s.rda", g))$tests |> mutate(sample = recode(scenario, Isotropic = "Symmetric"))
  G <- ld(sprintf("data/derived/half_bandwidth_grid_g%s.rda", g))$grid |> filter(generation >= 1904, !is.na(r2_pgd_dir))
  for (s in SYMS) {
    t <- T |> filter(sample == s); gg <- G |> filter(sample == s)
    for (a in c(0.01, 0.05, 0.10)) add(lab, sprintf("Source–sink, α = %.2f", a), s, clus(t$Replicate, t$p < a))
    add(lab, "ΔR² > 0", s, clus(gg$Replicate, gg$r2_pgd_dir > gg$r2_cgd))
    add(lab, "ΔR² > 0 and slope ratio < 1 (post hoc)", s, clus(gg$Replicate, gg$r2_pgd_dir > gg$r2_cgd & gg$rel_fwd < 1))
  }
}
fr <- ld("data/derived/fresh_replicates.rda")$fr_tests
for (st in c("1", "0.5")) for (s in c("Symmetric", "sym-low")) {
  t <- fr |> filter(setting == st, sample == s); lab <- sprintf("Fresh lineages, γ = %s", if (st == "0.5") "½" else "1")
  for (a in c(0.01, 0.05, 0.10)) add(lab, sprintf("Source–sink, α = %.2f", a), s, clus(t$Replicate, t$p < a))
  add(lab, "ΔR² > 0 (50-generation grid)", s, clus(t$Replicate[!is.na(t$dR2)], t$dR2[!is.na(t$dR2)] > 0))
}
ci <- bind_rows(rows)
save(ci, file = "data/derived/fpr_cluster_ci.rda")

pct <- function(x) sprintf("%.1f", 100 * x)
md <- c("# False-positive rates with lineage-clustered intervals", "",
  "Written by `R/analysis/fpr_cluster_ci.R`. 95% intervals from a cluster bootstrap over lineages (50 lineages resampled with replacement, 10,000 times; percentile interval). For comparison: the naive binomial interval over censuses, which treats censuses as independent, and the design effect (bootstrap variance ÷ binomial variance; 1 = no clustering). \"Lineages with ≥ 1\" counts lineages contributing any false positive; \"max\" is the highest single-lineage rate.", "",
  "Grids: source–sink, every 50 generations (1,000 censuses per sample; *burn-in* 1904–1999 every 5). ΔR², every 5 generations (1,000 *burn-in*, 10,000 per forward sample). Fresh lineages: every 50 generations for both tests.", "")
for (set in unique(ci$set)) {
  x <- ci[ci$set == set, ]
  md <- c(md, sprintf("## %s", set), "", "| Test | Sample | Censuses | Rate (%) | Lineage-clustered 95% CI | Binomial 95% CI | Design effect | Lineages with ≥ 1 / max |",
          "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
          sprintf("| %s | %s | %d | %s | %s–%s | %s–%s | %.1f | %d / %s%% |", x$test, x$sample, x$censuses, pct(x$rate), pct(x$lo), pct(x$hi), pct(x$bin_lo), pct(x$bin_hi),
                  x$deff, x$lin_with_any, pct(x$lin_max)), "")
}
writeLines(md, "output/fpr_cluster_ci.md"); cat("Wrote local/fpr_cluster_ci.md\n")
