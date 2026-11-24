# R/analysis/bandwidth_gamma_sweep.R
#
# Bandwidth gamma sweep: sensitivity analysis (Fig_BandwidthSweep).
#
# The bandwidth b_i = gamma * s_i / k_i is varied on a geometric sqrt(2) grid
# gamma = 2^(-2, -1.5, ..., 1) (labelled 0.25, 0.35, 0.50, 0.71, 1.00, 1.41, 2.00;
# 0.71 = 1/sqrt(2), where the kernel exponent becomes e^2 / b^2), plus two limits:
#   NN    w_{j|i} = 1 for i's nearest retained neighbour (smallest e_ij), ties split
#         equally: neighbourhood_weights() at gamma = 1e-6, where the log-sum-exp
#         kernel has reached the limit (checked below). Edges whose e_ij agree to
#         floating-point precision (relative gap < 1e-13; 86 of 7,000 censuses) are
#         treated as near-ties and share the weight, rather than going to one of them.
#   flat  w_{j|i} = 1 / k_i (gamma = Inf), so Delta = Delta0 = 1/k_j - 1/k_i and
#         S = -1/k centred within each component (checked).
# Every quantity is computed with stock gstudio (>= 1.15.0) at each setting:
# gravity_edges(), source_sink_scores(), source_sink_test() (Spearman r(S, x), Moran
# spectral randomization on the graph, B = 999, two-sided), ibgd(mode = "dcGD") and
# neighbourhood_weights(). The graphs, cGD edge weights and the test's null adjacency
# exp(-e^2 / 2 bbar^2) do not depend on gamma. The same seed is set before every test
# in a census, so every setting sees the same random stream.
#
# Samples (50 replicates each): Burn-in 1904-1999 every 5 generations (1,000);
# Symmetric, Redistributed, Obstructed, sym-mid, sym-low, sym-verylow 2004-2954
# every 50 (1,000 each). Phases from lp_phase() with the policy column
# phase_policy_R_gen (unflagged lineages stay in ascent).
#
# Run from the repository root:  Rscript R/analysis/bandwidth_gamma_sweep.R
#   PILOT=1 Rscript R/analysis/bandwidth_gamma_sweep.R   (unit checks + timing only)
# Outputs: data/derived/bandwidth_gamma_sweep.rda (`sweep`: per census x setting;
#          `censuses`; `checks`; `meta`), local/bandwidth_gamma_sweep_results.md
#          (tables and checks; the summary is written by hand above them),
#          media/fig-bandwidth-gamma-sweep-v4.png.

dir.create("output", showWarnings = FALSE)   # run reports (not tracked)
suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(tidyr)
  library(parallel); library(ggplot2) })
source("R/common/gravity_convention.R")
source("R/common/lineage_phases.R")

OUT  <-"data/derived/bandwidth_gamma_sweep.rda"
RES  <- "output/bandwidth_gamma_sweep_results.md"
FIG  <- "media/fig-bandwidth-gamma-sweep-v4.png"
SEED <- 20260930L; B <- 999L
GAMMA <- setNames(2^seq(-2, 1, 0.5), c("0.25", "0.35", "0.50", "0.71", "1.00", "1.41", "2.00"))
SETTINGS <- c("NN", names(GAMMA), "flat")
GSET <- c(NN = 1e-6, GAMMA, flat = Inf)                         # gamma passed to gstudio for each setting
SCN <- c(Symmetric = 1L, Redistributed = 2L, Obstructed = 3L, `sym-mid` = 4L, `sym-low` = 5L, `sym-verylow` = 6L)
SYM <- c("Burn-in", "Symmetric", "sym-mid", "sym-low", "sym-verylow")
RATES <- list(Symmetric = c(0.025, 0.025), Redistributed = c(0.040, 0.010), Obstructed = c(0.025, 0.010))
POPS <- sprintf("Pop%02d", 1:25)
for (f in c(OUT, FIG)) if (file.exists(f) && !nzchar(Sys.getenv("PILOT"))) cat("note:", f, "exists (cached / will not be overwritten)\n")

## ---- census jobs -----------------------------------------------------------------------
bw_jobs <- expand.grid(Replicate = 1:50, scenario = c("Isotropic", "Redistributed", "Obstructed"),
                       generation = seq(2004L, 2954L, 50L), stringsAsFactors = FALSE)
jobs <- rbind(
  expand.grid(Replicate = 1:50, sample = "Burn-in", generation = seq(1904L, 1999L, 5L), stringsAsFactors = FALSE),
  expand.grid(Replicate = 1:50, sample = names(SCN), generation = seq(2004L, 2954L, 50L), stringsAsFactors = FALSE))
old <- match(paste(jobs$Replicate, sub("^Symmetric$", "Isotropic", jobs$sample), jobs$generation),
             paste(bw_jobs$Replicate, bw_jobs$scenario, bw_jobs$generation))
jobs$seed <- ifelse(is.na(old), SEED + 100000L + seq_len(nrow(jobs)), SEED + old)

load_census <- function(j) gc_load_graph(j$Replicate, j$generation, if (j$sample == "Burn-in") NULL else SCN[[j$sample]])

## ---- one census, all settings -----------------------------------------------------------
one_census <- function(j) {
  g <- load_census(j); if (is.null(g)) return(NULL)
  x <- gc_xg(g); k <- degree(g)[names(x)]
  truth <- j$sample %in% names(RATES)
  if (truth) { off <- row(diag(25)) != col(diag(25)); r <- RATES[[j$sample]]
    T0 <- matrix(0, 25, 25, dimnames = list(POPS, POPS)); i <- 1:24
    T0[cbind(i, i + 1L)] <- r[1]; T0[cbind(i + 1L, i)] <- r[2]; TrT <- t(T0 / max(T0[off])) }
  d1 <- NULL; out <- list()
  for (s in c("1.00", setdiff(SETTINGS, "1.00"))) {
    row <- tryCatch({
      gam <- GSET[[s]]
      set.seed(j$seed); tst <- source_sink_test(g, x, gam)
      sc <- tst$scores; S <- signif(sc$S, 12); ed <- gravity_edges(g, gam, x)
      if (s == "1.00") d1 <- ed$Delta
      ib <- tryCatch(ibgd(g, x, mode = "dcGD", gamma = gam)$estimate, error = function(e) NULL)
      v <- c(r_S = unname(tst$statistic), p = tst$p.value,
             cor_Sk = suppressWarnings(cor(S, sc$degree, method = "spearman")),
             cor_SS0 = suppressWarnings(cor(S, sc$S0, method = "spearman", use = "complete.obs")),
             r2_d0 = attr(sc, "degree_share"), grad_frac = attr(sc, "gradient_share"), dbar = mean(ed$Delta),
             sign_agree1 = mean(sign(ed$Delta) == sign(d1)),
             rho_d1 = suppressWarnings(cor(ed$Delta, d1, method = "spearman")),
             dR2 = if (is.null(ib)) NA else unname(ib["adj R2 dcGD"] - ib["R2 cGD"]),
             rel_fwd = if (is.null(ib)) NA else unname(ib["slope ratio"]),
             r2_cgd = if (is.null(ib)) NA else unname(ib["R2 cGD"]), r2_pgd_dir = if (is.null(ib)) NA else unname(ib["adj R2 dcGD"]),
             w_rho_tr = NA, corr_fwd = NA,
             flat_err = if (s == "flat") max(abs(S - sc$S0), na.rm = TRUE) else NA,
             flat_cor = if (s == "flat" && max(sc$component) == 1) abs(cor(S, -1 / sc$degree)) else NA)
      if (truth) { W <- neighbourhood_weights(g, gam)[POPS, POPS]; Ws <- W / max(W[off])
        v["w_rho_tr"] <- suppressWarnings(cor(Ws[off], TrT[off], method = "spearman"))
        link <- which(as_adjacency_matrix(g, sparse = FALSE)[POPS, POPS][cbind(1:24, 2:25)] > 0)   # true links retained as edges
        if (length(link)) v["corr_fwd"] <- mean(W[cbind(link + 1L, link)] > W[cbind(link, link + 1L)]) }
      data.frame(setting = s, t(v), failed = FALSE)
    }, error = function(e) data.frame(setting = s, failed = TRUE, err = conditionMessage(e)))
    out[[s]] <- row
  }
  bind_rows(out) |> mutate(Replicate = j$Replicate, sample = j$sample, generation = j$generation,
                           edges = ecount(g), n_deg1 = sum(k == 1), n_iso = sum(k == 0), n_comp_graph = components(g)$no)
}

## ---- unit checks (fast; always run) -----------------------------------------------------
#' NN limit built directly (nearest retained neighbour, ties split) against
#' neighbourhood_weights(g, 1e-6); kernel rows sum to one; S invariant to rescaling e
unit_checks <- function(n = 12) {
  set.seed(1); pick <- jobs[sample(nrow(jobs), n), ]
  res <- lapply(seq_len(n), function(r) {
    j <- pick[r, ]; g <- load_census(j); if (is.null(g)) return(NULL)
    A <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE); adj <- A > 0
    Wnn <- matrix(0, nrow(A), ncol(A), dimnames = dimnames(A))
    for (i in which(rowSums(adj) > 0)) { e <- A[i, adj[i, ]]; m <- e == min(e); Wnn[i, adj[i, ]] <- m / sum(m) }
    diag(Wnn) <- NA
    g10 <- g; E(g10)$weight <- 10 * E(g)$weight
    data.frame(sample = j$sample, Replicate = j$Replicate, generation = j$generation,
               nn_vs_gamma_1e6 = max(abs(neighbourhood_weights(g, GSET[["NN"]]) - Wnn), na.rm = TRUE),
               rowsum_err = max(abs(rowSums(neighbourhood_weights(g, 1), na.rm = TRUE)[rowSums(adj) > 0] - 1)),
               scale10_g050 = max(abs(source_sink_scores(g10, 0.5)$S - source_sink_scores(g, 0.5)$S)),
               scale10_g200 = max(abs(source_sink_scores(g10, 2)$S - source_sink_scores(g, 2)$S)))
  })
  do.call(rbind, res)
}
checks <- unit_checks()

if (nzchar(Sys.getenv("PILOT"))) {
  options(width = 200); print(checks, digits = 3)
  t0 <- Sys.time(); pj <- jobs[c(1, 1500, 2600, 3600, 5000, 6900), ]
  pr <- bind_rows(lapply(seq_len(nrow(pj)), function(i) one_census(pj[i, ])))
  cat(sprintf("\n%.1f s per census (serial)\n", as.numeric(difftime(Sys.time(), t0, units = "secs")) / nrow(pj)))
  print(pr |> select(sample, Replicate, generation, setting, r_S, p, cor_Sk, grad_frac, dR2, w_rho_tr, corr_fwd, flat_err, failed), digits = 3)
  quit(save = "no")
}

## ---- run (cached) -----------------------------------------------------------------------
if (file.exists(OUT)) load(OUT) else {
  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(jobs)), function(i) one_census(jobs[i, ]),
                  mc.cores = max(1L, detectCores() - 1L), mc.preschedule = FALSE)
  bad <- vapply(res, function(r) inherits(r, "try-error"), TRUE)
  sweep <- bind_rows(res[!bad])
  censuses <- jobs |> mutate(loaded = paste(Replicate, sample, generation) %in% paste(sweep$Replicate, sweep$sample, sweep$generation),
                             worker_error = bad)
  meta <- list(seed = SEED, B = B, gamma = GAMMA, wall_min = as.numeric(difftime(Sys.time(), t0, units = "mins")),
               gstudio = as.character(packageVersion("gstudio")), date = as.character(Sys.Date()))
  save(sweep, censuses, checks, meta, file = OUT)
}

## ---- phases and labels ----------------------------------------------------------------
PH_LAB <- c(Reorganize = "ascent", Plateau = "plateau", `Fall apart` = "decay")
d <- sweep |> filter(!failed) |>
  mutate(setting = factor(setting, SETTINGS),
         phase = ifelse(sample %in% c("Redistributed", "Obstructed"),
                        PH_LAB[as.character(lp_phase(sample, Replicate, generation))], "all"),
         col = ifelse(phase == "all", sample, paste(sample, phase)))
PCOLS <- c("Symmetric", "Redistributed ascent", "Redistributed plateau", "Redistributed decay",
           "Obstructed ascent", "Obstructed plateau")
d_ph <- d |> filter(col %in% PCOLS) |> mutate(col = factor(col, PCOLS))
d_sym <- d |> filter(sample %in% SYM) |> mutate(sample = factor(sample, SYM))

## ---- tables ---------------------------------------------------------------------------------
md <- function(df) {
  df <- as.data.frame(df); hdr <- paste("|", paste(names(df), collapse = " | "), "|")
  sep <- paste("|", paste(c(":----", rep(":---:", ncol(df) - 1)), collapse = " | "), "|")
  body <- apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))
  c(hdr, sep, body, "")
}
pct <- function(v, dig = 1) sprintf(paste0("%.", dig, "f%%"), 100 * v)
num <- function(v, dig = 2) ifelse(is.na(v), "—", formatC(v, format = "f", digits = dig))
wide <- function(x, colvar, val) x |> select(setting, all_of(colvar), all_of(val)) |>
  pivot_wider(names_from = all_of(colvar), values_from = all_of(val)) |> arrange(setting) |>
  mutate(setting = recode(as.character(setting), `0.71` = "0.71 (1/√2)", `1.00` = "1.00 (canonical)",
                          NN = "NN limit", flat = "flat limit")) |> rename(γ = setting)

n_cens <- d_ph |> filter(setting == "1.00") |> count(col) |> mutate(lab = sprintf("%s (n=%d)", col, n))
ph_names <- setNames(n_cens$lab, n_cens$col)

# A1
fpr <- d_sym |> group_by(setting, sample) |>
  summarise(v = sprintf("%s / %s / %s", pct(mean(p < .01)), pct(mean(p < .05)), pct(mean(p < .10))),
            f05 = mean(p < .05), .groups = "drop")
flag <- fpr |> filter(f05 > 0.075)
dec <- d_sym |> filter(sample == "Symmetric") |> mutate(bin = cut(p, seq(0, 1, .1), include.lowest = TRUE)) |>
  count(setting, bin) |> group_by(setting) |> mutate(v = pct(n / sum(n), 0)) |> ungroup()
# A2-A4
ss_ph <- d_ph |> group_by(setting, col) |>
  summarise(rej = pct(mean(p < .05)), med = num(median(r_S)), up = pct(mean(r_S < 0), 0), .groups = "drop") |>
  left_join(d_ph |> group_by(setting, col, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |>
              group_by(setting, col) |> summarise(lin = pct(mean(m < 0), 0), nl = n(), .groups = "drop"),
            by = c("setting", "col")) |> mutate(col = ph_names[as.character(col)])
# A5
lean <- d_sym |> filter(sample == "Symmetric") |> group_by(setting, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |>
  group_by(setting) |> summarise(`lineages with mean r > 0` = sprintf("%d / %d", sum(m > 0), n()),
                                 `sign-test p` = signif(binom.test(sum(m > 0), sum(m != 0))$p.value, 2), .groups = "drop")
# B6-B9
fl_ph <- d_ph |> group_by(setting, col) |>
  summarise(Sk = num(median(abs(cor_Sk), na.rm = TRUE)), SS0 = num(median(cor_SS0, na.rm = TRUE)),
            r2d0 = num(median(r2_d0, na.rm = TRUE)), grad = num(median(grad_frac, na.rm = TRUE), 3), .groups = "drop") |>
  mutate(col = ph_names[as.character(col)])
# C10-C11
c10 <- d_ph |> group_by(setting, col) |>
  summarise(v = sprintf("%s / %s", pct(mean(sign_agree1), 0), num(median(rho_d1, na.rm = TRUE))), .groups = "drop") |>
  mutate(col = ph_names[as.character(col)])
c11 <- d_ph |> filter(grepl("plateau", col)) |> group_by(setting, col, Replicate) |> summarise(m = mean(dbar), .groups = "drop") |>
  group_by(setting, col) |> summarise(v = sprintf("%s (%d)", pct(mean(m > 0), 0), n()), .groups = "drop")
# D12-D14
d12 <- d_sym |> group_by(setting, sample) |> summarise(v = pct(mean(dR2 > 0, na.rm = TRUE)), .groups = "drop")
d13 <- d_ph |> group_by(setting, col) |>
  summarise(det = pct(mean(dR2 > 0, na.rm = TRUE)),
            lt1 = pct(mean(rel_fwd[dR2 > 0 & !is.na(dR2)] < 1), 0), .groups = "drop") |> mutate(col = ph_names[as.character(col)])
# E15-E16
e15 <- d_ph |> filter(col %in% c("Symmetric", "Redistributed plateau", "Obstructed plateau")) |> group_by(setting, col) |>
  summarise(v = sprintf("%s / %s", num(median(w_rho_tr, na.rm = TRUE)), pct(mean(corr_fwd, na.rm = TRUE), 0)), .groups = "drop")

## ---- consistency checks ----------------------------------------------------------------------
# Against the run with the project's own implementation (before gstudio 1.15), when its
# output is present: the statistics must agree to machine precision; p-values come from
# new null draws and agree only within Monte Carlo error.
PRE <- "data/derived/pre_gstudio115/bandwidth_gamma_sweep.rda"
pre <- if (file.exists(PRE)) { e <- new.env(); load(PRE, envir = e)
  inner_join(e$sweep |> filter(!failed), sweep |> filter(!failed), by = c("sample", "Replicate", "generation", "setting"), suffix = c(".pre", ""))
} else NULL
STATS <- c("r_S", "cor_Sk", "cor_SS0", "r2_d0", "grad_frac", "dbar", "rho_d1", "dR2", "rel_fwd", "w_rho_tr", "corr_fwd")
pre_stat <- if (!is.null(pre)) sapply(STATS, function(v) max(abs(pre[[v]] - pre[[paste0(v, ".pre")]]), na.rm = TRUE))
pre_rej <- if (!is.null(pre)) pre |> mutate(col = ifelse(sample %in% c("Redistributed", "Obstructed"),
                                                         paste(sample, PH_LAB[as.character(lp_phase(sample, Replicate, generation))]), sample)) |>
  filter(col %in% c("Symmetric", "Redistributed plateau", "Obstructed plateau"), setting %in% c("0.50", "1.00")) |>
  group_by(setting, col) |> summarise(n = n(), pre = pct(mean(p.pre < .05)), now = pct(mean(p < .05)), .groups = "drop")
flat <- d |> filter(setting == "flat")

## ---- results markdown ---------------------------------------------------------------------
L <- c("<!-- TABLES BELOW ARE WRITTEN BY R/analysis/bandwidth_gamma_sweep.R -->", "",
  sprintf("Run %s; %d censuses loaded of %d; %d census x setting rows failed; wall time %.0f min; B = %d; gstudio %s.",
          meta$date, sum(censuses$loaded), nrow(censuses), sum(sweep$failed), meta$wall_min, meta$B, meta$gstudio), "",
  "Grid: γ = 2^(−2, −1.5, …, 1), i.e. 0.25, 0.354, 0.5, 0.707, 1, 1.414, 2 (labelled to two decimals). Columns give censuses (n) at γ = 1.", "",
  "## Consistency checks", "",
  "**Unit checks** (12 random censuses; max absolute difference):", "",
  md(checks |> mutate(across(where(is.numeric) & !c(Replicate, generation), ~ signif(.x, 2))) |> rename_with(~ gsub("_", " ", .x))),
  if (!is.null(pre)) c(
    sprintf("**Against the pre-gstudio-1.15 run** (%d census × setting rows): largest absolute difference in each statistic:", nrow(pre)), "",
    md(data.frame(statistic = names(pre_stat), `max difference` = signif(pre_stat, 2), check.names = FALSE)),
    "Rejection at α = 0.05 (new null draws):", "", md(pre_rej |> rename(γ = setting, `sample / phase` = col, `before` = pre, `gstudio 1.15` = now))),
  sprintf("**Flat limit:** max |S − S⁰| (S⁰ = −1/k centred within component) = %.2g over %d censuses; |cor(S, −1/k)| = 1 in %d of %d single-component censuses (min %.6f); %d censuses have more than one component, where only the within-component identity applies.",
          max(flat$flat_err, na.rm = TRUE), nrow(flat), sum(flat$flat_cor > 1 - 1e-10, na.rm = TRUE), sum(!is.na(flat$flat_cor)),
          min(flat$flat_cor, na.rm = TRUE), sum(is.na(flat$flat_cor))), "",
  sprintf("**Degree-one populations** (kept at every γ): median %g per census (range %d–%d); isolated nodes: %d censuses.",
          median(flat$n_deg1), min(flat$n_deg1), max(flat$n_deg1), sum(flat$n_iso > 0)), "",
  if (nrow(flag)) c("**FLAG: false-positive rate above 7.5% at α = 0.05:**", "", md(flag |> transmute(γ = setting, sample, FPR = pct(f05))))
  else "**No γ or limit has a *Symmetric* or sym-* false-positive rate above 7.5% at α = 0.05.**", "",
  "## A. Source–sink gradient test", "",
  "**A1. False-positive rate, α = 0.01 / 0.05 / 0.10, symmetric samples**", "", md(wide(fpr, "sample", "v")),
  "**A1b. *Symmetric* p-value deciles** (share of 1,000 censuses; 10% each under a uniform null)", "", md(wide(dec, "bin", "v")),
  "**A2. Rejection rate at α = 0.05** (*Symmetric* column = false-positive rate)", "", md(wide(ss_ph, "col", "rej")),
  "**A3. Median r(S,x) / censuses with r < 0 (sources upstream)**", "",
  md(wide(ss_ph |> mutate(v = paste(med, "/", up)), "col", "v")),
  "**A4. Lineages whose mean r(S,x) < 0**", "", md(wide(ss_ph, "col", "lin")),
  "**A5. Lineage lean, *Symmetric***", "", md(lean |> mutate(setting = as.character(setting)) |> rename(γ = setting)),
  "## B. Topological floor", "",
  "**B6. Median |Spearman cor(S, k)|**", "", md(wide(fl_ph, "col", "Sk")),
  "**B7. Median Spearman cor(S, S⁰)**", "", md(wide(fl_ph, "col", "SS0")),
  "**B8. Median edge-level R² of Δ on Δ⁰**", "", md(wide(fl_ph, "col", "r2d0")),
  "**B9. Median share of Δ variance reproduced by the gradient S_i − S_j** (source_sink_scores() `gradient_share`)", "", md(wide(fl_ph, "col", "grad")),
  "## C. Edge and graph level", "",
  "**C10. Sign agreement of Δ with γ = 1 / median Spearman ρ of signed Δ with γ = 1**", "", md(wide(c10, "col", "v")),
  "**C11. Lineages with mean Δ̄ > 0 on the *plateau* (lineages)**", "", md(wide(c11, "col", "v")),
  "## D. Directional isolation (ΔR², pGD from the γ-specific w)", "",
  "**D12. False-positive rate, ΔR² > 0, symmetric samples**", "", md(wide(d12, "sample", "v")),
  "**D13. Detection (ΔR² > 0) by scenario and phase**", "", md(wide(d13, "col", "det")),
  "**D14. Share of detections with b_fwd/b_rev < 1**", "", md(wide(d13, "col", "lt1")),
  "## E. Truth comparison (computed directly from w and the analytic migration matrix, not from divMigrate output)", "",
  "**E15–E16. Median Spearman ρ of w_{j|i} with the transposed truth / share of retained true links on which w is larger in the imposed direction**", "",
  md(wide(e15, "col", "v")))
summ <- if (file.exists(RES)) { x <- readLines(RES); i <- grep("^<!-- TABLES BELOW", x); if (length(i)) x[seq_len(i[1] - 1)] else x } else c("# Bandwidth γ sweep: results", "", "_Summary to be written._", "")
writeLines(c(summ, L), RES)
cat("Wrote", RES, "\n")

## ---- figure -------------------------------------------------------------------------------------
xpos <- c(NN = -3, setNames(log2(GAMMA), names(GAMMA)), flat = 2)
COL <- c(`Burn-in` = "#eda100", Symmetric = "#2a78d6", `sym-mid` = "#e87ba4", `sym-low` = "#008300",
         `sym-verylow` = "#4a3aa7", Redistributed = "#eb6834", Obstructed = "#1baf7a")
SHP <- c(`Burn-in` = 15, Symmetric = 16, `sym-mid` = 17, `sym-low` = 18, `sym-verylow` = 4, Redistributed = 16, Obstructed = 17)
fs <- function(x, panel, series, stat) x |> mutate(panel = panel, series = series, x = xpos[as.character(setting)]) |> rename(y = all_of(stat))
pd <- bind_rows(
  d_sym |> group_by(setting, series = as.character(sample)) |> summarise(y = mean(p < .05), .groups = "drop") |> mutate(panel = "a", what = "FPR"),
  d_ph |> filter(grepl("plateau", col)) |> group_by(setting, series = sub(" plateau", "", col)) |>
    summarise(y = mean(p < .05), .groups = "drop") |> mutate(panel = "b", what = "plateau"),
  d_ph |> filter(col %in% c("Symmetric", "Redistributed plateau", "Obstructed plateau")) |>
    group_by(setting, series = sub(" plateau", "", col)) |> summarise(y = median(abs(cor_Sk), na.rm = TRUE), .groups = "drop") |>
    mutate(panel = "c", what = ifelse(series == "Symmetric", "whole forward", "plateau")),
  d_sym |> group_by(setting, series = as.character(sample)) |> summarise(y = mean(dR2 > 0, na.rm = TRUE), .groups = "drop") |> mutate(panel = "d", what = "FPR"),
  d_ph |> filter(grepl("plateau", col)) |> group_by(setting, series = sub(" plateau", "", col)) |>
    summarise(y = mean(dR2 > 0, na.rm = TRUE), .groups = "drop") |> mutate(panel = "d", what = "plateau"),
  d_ph |> filter(col %in% c("Redistributed plateau", "Obstructed plateau", "Redistributed decay")) |>
    group_by(setting, col, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |>
    group_by(setting, col) |> summarise(y = mean(m < 0), .groups = "drop") |>
    mutate(panel = "e", series = sub(" (plateau|decay)", "", col), what = ifelse(grepl("decay", col), "decay", "plateau")) |> select(-col)) |>
  mutate(x = xpos[as.character(setting)], lim = setting %in% c("NN", "flat"),
         lt = ifelse(what == "decay", "dashed", "solid"))
PANEL <- c(a = "(a) Source–sink false-positive rate (α = 0.05)", b = "(b) Source–sink plateau power",
           c = "(c) Median |cor(S, k)|", d = "(d) ΔR² > 0: FPR (symmetric) and plateau detection",
           e = "(e) Lineages placing sources upstream")
pd$panel <- factor(PANEL[pd$panel], PANEL)
ref <- data.frame(panel = factor(PANEL["a"], PANEL), y = 0.05)
p <- ggplot(pd, aes(x, y, colour = series, shape = series)) +
  geom_vline(data = data.frame(x = c(0, -1), lt = c("solid", "dashed")), aes(xintercept = x, linetype = lt),
             colour = "grey70", linewidth = 0.4, inherit.aes = FALSE) +
  geom_hline(data = ref, aes(yintercept = y), colour = "grey50", linewidth = 0.4, linetype = "dotted") +
  geom_line(data = filter(pd, !lim), aes(group = interaction(series, what), linetype = lt), linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = COL, breaks = names(COL), name = NULL) +
  scale_shape_manual(values = SHP, breaks = names(COL), name = NULL) +
  scale_linetype_identity() +
  scale_x_continuous(breaks = xpos, labels = c("NN", names(GAMMA), "flat"), minor_breaks = NULL) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  labs(x = "γ (log scale; NN and flat limits at the edges; solid line γ = 1, dashed γ = 0.5; dashed series = Redistributed decay)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), strip.text = element_text(hjust = 0, face = "bold"),
        axis.title.x = element_text(size = 9, colour = "grey30"))
if (!file.exists(FIG)) { ggsave(FIG, p, width = 8, height = 13, dpi = 150); cat("Wrote", FIG, "\n") } else cat(FIG, "exists; not overwritten\n")
