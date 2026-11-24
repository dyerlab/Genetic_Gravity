# R/case_study/araptus_design_calibration.R
#
# Worked example of design-specific calibration of the directional isolation comparison,
# using the simulated stepping stone (illustrative: the simulation is not Baja). The
# Araptus central component (14 sites; per-site n 5-10; 8 loci) is mimicked by:
#   K = 14 demes, stratified spread draw along the 25-deme chain (as in density_power_full.R);
#   n per deme = the 14 Araptus sample sizes in axis order (distance from the refuge, which
#     is position along the hypothesised expansion), assigned to the drawn demes in index order;
#   L = 8 loci drawn from those polymorphic in the census.
# Null: Symmetric forward censuses (50 lineages x 6, the density-power census list) and, as
# the strong-drift check, sym-low (50 x 4); N_DRAW draws per census. The cutoff is the
# 95th percentile of the Symmetric Delta R^2, floored at 0; one-sided empirical p of the
# observed Delta R^2 = (1 + #{null >= obs}) / (1 + B). Power at the cutoff from the plateau
# cores of Redistributed and Obstructed (the same census list).
#
# Run: Rscript R/case_study/araptus_design_calibration.R   (env N_DRAW, default 5)
# Output: data/derived/araptus_design_calibration.rda; printed summary.

suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(parallel) })
source("R/common/gravity_convention.R")
N_DRAW <- as.integer(Sys.getenv("N_DRAW", "5")); GAMMA <- 0.5
SCN <- c(Symmetric = 1L, Redistributed = 2L, Obstructed = 3L, `sym-low` = 5L)
OUT <- "data/derived/araptus_design_calibration.rda"

## ---- the Araptus design ----------------------------------------------------------------------
ar <- new.env(); load("data/case_study_arapat.rda", envir = ar)
sub <- ar$results$South$graph; nm <- V(sub)$name
n_axis <- as.integer(table(as.character(ar$peninsula$Population))[nm[order(ar$d[nm])]])   # refuge outward
obs <- ibgd(sub, mode = "dcGD", distance = abs(outer(ar$d[nm], ar$d[nm], "-")),
            forward = outer(ar$d[nm], ar$d[nm], function(di, dj) dj > di), gamma = GAMMA)
dR2_obs <- unname(obs$statistic)
cat(sprintf("Araptus central component: K = %d, n = %s (mean %.1f), L = 8; observed Delta R^2 = %.4f\n\n",
            length(nm), paste(n_axis, collapse = ","), mean(n_axis), dR2_obs))

## ---- null and alternative graphs -------------------------------------------------------------
load("data/derived/density_power_full.rda")          # `cens`: the panel's census list (seed 20261123)
draw_spread <- function(K) vapply(split(1:25, cut(1:25, K, labels = FALSE)), function(b) b[sample.int(length(b), 1L)], 1L)
one_census <- function(i) {
  cz <- cens[i, ]; set.seed(20261125L + i)
  e <- new.env(); load(sprintf("%s/replicate%d/rep%d-genotypes-scenario%d-%d.rda", GC_SIM_DIR, cz$Replicate, cz$Replicate,
                               SCN[[cz$scenario]], cz$generation), envir = e); gt <- e$genotypes
  loci <- setdiff(names(gt), c("ID", "Population"))
  poly <- loci[vapply(loci, function(l) length(unique(as.character(gt[[l]]))) > 1, TRUE)]
  bind_rows(lapply(seq_len(N_DRAW), function(d) {
    demes <- sort(draw_spread(14L)); pops <- sprintf("Pop%02d", demes)
    keep <- unlist(mapply(function(p, n) sample(which(gt$Population == p), n), pops, n_axis, SIMPLIFY = FALSE))
    L <- sample(poly, min(8L, length(poly))); s <- gt[keep, ]
    g <- tryCatch(popgraph(to_mv(s[, L, drop = FALSE]), groups = as.factor(as.character(s$Population))), error = function(e) NULL)
    if (is.null(g)) return(NULL)
    x <- gc_xg(g); ib <- tryCatch(ibgd(g, x, mode = "dcGD", gamma = GAMMA), error = function(e) NULL)
    t_s <- tryCatch(source_sink_test(g, x, GAMMA), error = function(e) NULL)
    data.frame(cz, draw = d, n_edges = ecount(g), mean_deg = mean(degree(g)), n_comp = components(g)$no,
               dR2 = if (is.null(ib)) NA else unname(ib$statistic), ratio = if (is.null(ib)) NA else ib$estimate[["slope ratio"]],
               p_S = if (is.null(t_s)) NA else t_s$p.value, r_S = if (is.null(t_s)) NA else unname(t_s$statistic),
               deg_share = if (is.null(t_s)) NA else attr(t_s$scores, "degree_share"))
  }))
}
if (!file.exists(OUT)) {
  t0 <- Sys.time()
  cal <- bind_rows(mclapply(seq_len(nrow(cens)), function(i) tryCatch(one_census(i), error = function(e) NULL),
                            mc.cores = max(1L, detectCores() - 1L), mc.preschedule = FALSE))
  save(cal, n_axis, dR2_obs, file = OUT)
  cat(sprintf("Wrote %s: %d graphs in %.1f min\n\n", OUT, nrow(cal), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
} else load(OUT)

## ---- summary ---------------------------------------------------------------------------------
boot_q <- function(z, lin, prob = .95, B = 1000L) { m <- !is.na(z); z <- z[m]; lin <- lin[m]; ids <- unique(lin); sp <- split(z, lin)
  bs <- replicate(B, quantile(unlist(sp[as.character(sample(ids, length(ids), TRUE))]), prob)); c(est = unname(quantile(z, prob)), lo = unname(quantile(bs, .025)), hi = unname(quantile(bs, .975))) }
boot_rate <- function(y, lin, B = 1000L) { m <- !is.na(y); y <- y[m]; lin <- lin[m]; s <- tapply(y, lin, sum); n <- tapply(y, lin, length); L <- length(s)
  bs <- replicate(B, { i <- sample.int(L, L, TRUE); sum(s[i]) / sum(n[i]) }); c(est = mean(y), lo = unname(quantile(bs, .025)), hi = unname(quantile(bs, .975))) }
set.seed(20261126L)
sy <- cal |> filter(scenario == "Symmetric"); sl <- cal |> filter(scenario == "sym-low")
q <- boot_q(sy$dR2, sy$Replicate); cutoff <- max(0, q[["est"]])
p_emp <- (1 + sum(sy$dR2 >= dR2_obs, na.rm = TRUE)) / (1 + sum(!is.na(sy$dR2)))
p_low <- (1 + sum(sl$dR2 >= dR2_obs, na.rm = TRUE)) / (1 + sum(!is.na(sl$dR2)))
cat(sprintf("Null graphs: Symmetric %d (from %d lineages), sym-low %d; fragmented: %.0f%% / %.0f%%\n", sum(!is.na(sy$dR2)),
            n_distinct(sy$Replicate), sum(!is.na(sl$dR2)), 100 * mean(sy$n_comp > 1), 100 * mean(sl$n_comp > 1)))
cat(sprintf("Uncalibrated rule Delta R^2 > 0 under Symmetric: %.3f (%.3f-%.3f); sym-low: %.3f\n",
            boot_rate(sy$dR2 > 0, sy$Replicate)[1], boot_rate(sy$dR2 > 0, sy$Replicate)[2], boot_rate(sy$dR2 > 0, sy$Replicate)[3], mean(sl$dR2 > 0, na.rm = TRUE)))
cat(sprintf("Symmetric 95th percentile of Delta R^2: %.4f (lineage bootstrap %.4f-%.4f); cutoff = %.4f\n", q[1], q[2], q[3], cutoff))
cat(sprintf("sym-low 95th percentile: %.4f; FPR at the Symmetric cutoff: %.3f\n", quantile(sl$dR2, .95, na.rm = TRUE), mean(sl$dR2 > cutoff, na.rm = TRUE)))
cat(sprintf("Observed Araptus Delta R^2 = %.4f: %s the cutoff; empirical one-sided p = %.4f (Symmetric), %.4f (sym-low)\n",
            dR2_obs, if (dR2_obs > cutoff) "EXCEEDS" else "does not exceed", p_emp, p_low))
cat(sprintf("Percentile of the observed value in the Symmetric null: %.1f%%\n\n", 100 * mean(sy$dR2 < dR2_obs, na.rm = TRUE)))
pw <- cal |> filter(scenario %in% c("Redistributed", "Obstructed")) |> group_by(scenario) |>
  group_modify(~ { v <- boot_rate(.x$dR2 > cutoff, .x$Replicate); s <- boot_rate(.x$p_S < .05, .x$Replicate)
    tibble(G = nrow(.x), power_dR2 = v[1], lo = v[2], hi = v[3], power_S = s[1], ratio_lt1 = mean(.x$ratio[.x$dR2 > cutoff] < 1, na.rm = TRUE),
           dR2_median = median(.x$dR2, na.rm = TRUE)) }) |> ungroup()
cat("Plateau power at the cutoff (and the source-sink test at alpha = 0.05):\n"); print(as.data.frame(pw), digits = 3, row.names = FALSE)
cat(sprintf("\nSource-sink test, Symmetric FPR: %.3f; Araptus r(S,d) = -0.55, p = 0.44\n", mean(sy$p_S < .05, na.rm = TRUE)))
cat(sprintf("Null degree share (median): %.2f; mean degree %.1f (Araptus central: 0.85, 2.7)\n", median(sy$deg_share, na.rm = TRUE), median(sy$mean_deg)))
