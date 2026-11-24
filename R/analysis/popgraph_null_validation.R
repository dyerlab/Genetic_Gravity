# R/analysis/popgraph_null_validation.R
#
# Validation of gstudio::popgraph_null()'s Wright-Fisher generator against the
# individual-based simulations behind the density/power panel. The generator is run on
# the manuscript's Symmetric model exactly: 25 demes, symmetric stepping stone at 0.025 per
# direction (one neighbour at each end), N = Ne = 100, 20 biallelic loci started at
# p ~ N(0.5, 0.1) clipped to [0.01, 0.99] in every deme, sampled after 2000 + (4..954)
# generations (the Symmetric census window). Each simulated state is subsampled as in the
# panel, and the null quantiles of Delta R^2 (x = deme index) are compared with the panel's:
#   araptus  K = 14 spread, n = the Araptus per-site sizes in axis order, L = 8 polymorphic
#            (panel: R/case_study/araptus_design_calibration.R)
#   K25n30   K = 25, n = 30, all loci (panel: density_power_full.R)
#   K12n30   K = 12 spread, n = 30, all loci
# Run: Rscript R/analysis/popgraph_null_validation.R  -> data/derived/popgraph_null_validation.rda
suppressPackageStartupMessages({ library(gstudio); library(igraph); library(dplyr); library(parallel) })
options(width = 200)
B <- as.integer(Sys.getenv("B", "300")); OUT <- "data/derived/popgraph_null_validation.rda"
K <- 25L; pops <- sprintf("Pop%02d", 1:K)
M <- matrix(0, K, K, dimnames = list(pops, pops)); i <- 1:(K - 1)
M[cbind(i, i + 1)] <- 0.025; M[cbind(i + 1, i)] <- 0.025; diag(M) <- 1 - rowSums(M)
load("data/derived/araptus_design_calibration.rda")      # n_axis, and `cal` (the panel's Araptus-design null)
draw_spread <- function(k) vapply(split(1:25, cut(1:25, k, labels = FALSE)), function(b) b[sample.int(length(b), 1L)], 1L)
dR2_of <- function(P, demes, n, L) {
  Pl <- lapply(P[L], function(p) p[demes, , drop = FALSE])
  n_i <- setNames(rep_len(n, length(demes)), pops[demes])
  X <- gstudio:::.pgn_genotypes(Pl, n_i)
  g <- tryCatch(suppressWarnings(popgraph(X$mv, groups = X$groups)), error = function(e) NULL)
  if (is.null(g)) return(c(dR2 = NA, n_comp = NA))
  x <- setNames(as.numeric(sub("\\D+", "", V(g)$name)), V(g)$name)
  ib <- tryCatch(ibgd(g, x, mode = "dcGD", gamma = 0.5), error = function(e) NULL)
  c(dR2 = if (is.null(ib)) NA else unname(ib$statistic), n_comp = components(g)$no)
}
one <- function(b) {
  set.seed(20261127L + b)
  p <- pmin(pmax(rnorm(20, 0.5, 0.1), 0.01), 0.99); p0 <- lapply(p, function(q) c(q, 1 - q)); names(p0) <- sprintf("L%02d", 1:20)
  ngen <- 2000L + sample(seq(4L, 954L, 50L), 1L)
  P <- gstudio:::.pgn_wf(p0, M, 100L, ngen)
  poly <- names(P)[vapply(P, function(Pl) any(colMeans(Pl) > 0 & colMeans(Pl) < 1), TRUE)]
  rbind(data.frame(design = "araptus", t(dR2_of(P, sort(draw_spread(14L)), n_axis, sample(poly, min(8L, length(poly)))))),
        data.frame(design = "K25n30", t(dR2_of(P, 1:25, 30L, names(P)))),
        data.frame(design = "K12n30", t(dR2_of(P, sort(draw_spread(12L)), 30L, names(P))))) |>
    mutate(b = b, ngen = ngen, n_poly = length(poly))
}
if (!file.exists(OUT)) {
  t0 <- Sys.time(); val <- bind_rows(mclapply(seq_len(B), one, mc.cores = max(1L, detectCores() - 1L)))
  save(val, file = OUT); cat(sprintf("Wrote %s in %.1f min\n", OUT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
} else load(OUT)
load("data/derived/density_power_full.rda")
ref <- bind_rows(
  cal |> filter(scenario == "Symmetric") |> transmute(design = "araptus", dR2, n_comp, Replicate),
  res |> filter(scenario == "Symmetric", scheme == "spread", n == 30, K == 25) |> transmute(design = "K25n30", dR2, n_comp, Replicate),
  res |> filter(scenario == "Symmetric", scheme == "spread", n == 30, K == 12) |> transmute(design = "K12n30", dR2, n_comp, Replicate))
set.seed(20261128L)
bq <- function(z, g, pr) { sp <- split(z, g); ids <- names(sp); bs <- replicate(1000, quantile(unlist(sp[sample(ids, length(ids), TRUE)]), pr, na.rm = TRUE))
  sprintf("%.4f (%.4f-%.4f)", quantile(z, pr, na.rm = TRUE), quantile(bs, .025), quantile(bs, .975)) }
cat("\nDelta R^2 under the symmetric model: Wright-Fisher generator vs the individual-based panel\n")
for (d in c("araptus", "K25n30", "K12n30")) {
  w <- val |> filter(design == d); r <- ref |> filter(design == d)
  cat(sprintf("\n%s\n  q95  WF %s | panel %s\n  q90  WF %s | panel %s\n  median WF %.4f | panel %.4f\n  P(dR2>0) WF %.3f | panel %.3f;  fragmented WF %.2f | panel %.2f;  KS p = %.3f\n",
    d, bq(w$dR2, w$b, .95), bq(r$dR2, r$Replicate, .95), bq(w$dR2, w$b, .90), bq(r$dR2, r$Replicate, .90),
    median(w$dR2, na.rm = TRUE), median(r$dR2, na.rm = TRUE), mean(w$dR2 > 0, na.rm = TRUE), mean(r$dR2 > 0, na.rm = TRUE),
    mean(w$n_comp > 1, na.rm = TRUE), mean(r$n_comp > 1, na.rm = TRUE), suppressWarnings(ks.test(w$dR2, r$dR2)$p.value)))
}
