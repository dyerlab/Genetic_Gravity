# exploratory/phase_boundaries.R
#
# Per-lineage phase boundaries (reorganize -> plateau -> fall apart) and the
# inferences that survive in each phase.
#
# Boundaries, per replicate x scenario (Isotropic run through the same rules as
# the control for how often a rule fires on noise):
#   End of reorganization (candidates)
#     R_gain   t90 of the smoothed within-pGD directional gain (direction-aware
#              minus direction-blind pGD R^2), from its burn-in baseline to its
#              maximum (exploratory/lineage_phase_windows.R; lineage_phase_windows.csv).
#     R_dR2    t90 of the smoothed pGD - cGD Delta R^2 (same file).
#   Start of falling apart (candidates)
#     F_hzn    informative horizon: first census with fewer than 4.5 loci
#              polymorphic per deme (R/information_loss.R; below it, information
#              loss alone inflates |Delta| by > 10%). Independent of every
#              statistic being evaluated.
#     F_hzn_s  the same on a penalized-spline smooth of n_poly_pop.
#     F_gain   first census after the gain's smoothed peak at which the smoothed
#              gain drops out of the top 10% of its rise (the mirror of t90).
#     F_diam   smoothed diameter 20% below its maximum, not before 2350
#              (exploratory/phase_corridor_groups.R).
#     F_abs    smoothed peak of mean |Delta| by 2950 (Fig_PeakAlignment).
#     F_nm     divMigrate breakdown: first of three consecutive censuses with a
#              non-finite or negative Nm entry.
#
# Primary scheme (chosen below, see results): reorganize = 2004 .. R_gain - 5,
# plateau = R_gain .. F_hzn - 5, fall apart = F_hzn onward. Lineages that never
# cross the horizon have no fall-apart phase.
#
# Run from the repository root:  Rscript exploratory/phase_boundaries.R
# Output: data/derived/phase_boundaries.rda (bounds, by_phase tables)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(mgcv) })
source("exploratory/gravity_convention.R")
grid <- seq(2004, 2999, 5); SC <- c("Isotropic", "Redistributed", "Obstructed")
HZN <- 4.5

fwd <- gc_read_forward() |> filter(Scenario %in% SC)
lpw <- read.csv("data/derived/lineage_phase_windows.csv", stringsAsFactors = FALSE)
ef <- new.env(); load("data/derived/ibgd_trajectory_fit.rda", envir = ef)
e  <- new.env(); load("data/derived/node_potential_reanalysis.rda", envir = e); rean <- e$rean

sm <- function(g, y, k = 10) as.numeric(predict(gam(y ~ s(g, k = k), data = data.frame(g = g, y = y), method = "REML"),
                                                 data.frame(g = grid)))
first_run <- function(b, gens, k = 3) { r <- rle(b); ends <- cumsum(r$lengths)
  i <- which(r$values & r$lengths >= k); if (length(i)) gens[ends[i[1]] - r$lengths[i[1]] + 1] else NA_real_ }

gain <- ef$traj |> filter(status == "ok") |> mutate(gain = r2_pgd_dir - r2_pgd_blind)
gbase <- gain |> filter(Treatment == "Burn-in") |> group_by(Replicate) |> summarise(b = mean(gain))
nmx <- rean$truth$d |> mutate(bad = Exclusion != "")

bounds <- list()
for (s in SC) for (r in 1:50) {
  a <- fwd |> filter(Scenario == s, replicate == r) |> arrange(generation)
  hz <- a$generation[a$n_poly_pop < HZN][1]
  hzs <- { f <- sm(a$generation, a$n_poly_pop); grid[f < HZN][1] }
  fd <- { f <- sm(a$generation, a$diameter); im <- which.max(f)
          i <- which(seq_along(f) > im & f <= 0.8 * f[im] & grid >= 2350)[1]; if (is.na(i)) NA else grid[i] }
  x <- gain |> filter(Treatment == s, Replicate == r) |> arrange(generation)
  fg <- { f <- sm(x$generation, x$gain); b <- gbase$b[gbase$Replicate == r]; im <- which.max(f)
          thr <- b + 0.9 * (f[im] - b); i <- which(seq_along(f) > im & f < thr)[1]; if (is.na(i)) NA else grid[i] }
  n <- nmx |> filter(Treatment == s, Replicate == r) |> arrange(Generation)
  fnm <- first_run(n$bad, n$Generation)
  l <- lpw[lpw$Scenario == s & lpw$Replicate == r, ]
  bounds[[length(bounds) + 1]] <- data.frame(Scenario = s, Replicate = r,
    R_gain = l$gain_t90, R_dR2 = l$dR2_t90, reorg = l$reorg_gain,
    F_hzn = hz, F_hzn_s = hzs, F_gain = fg, F_diam = fd,
    F_abs = if (isTRUE(l$turned) && l$absdelta_peak_gen <= 2950) l$absdelta_peak_gen else NA, F_nm = fnm)
}
bounds <- do.call(rbind, bounds)

## ---- phase labels per census (primary scheme) -------------------------------------
phase_of <- function(scn, rep, gen, R = "R_gain", F = "F_hzn") {
  b <- bounds[match(paste(scn, rep), paste(bounds$Scenario, bounds$Replicate)), ]
  Rg <- ifelse(is.na(b[[R]]), Inf, b[[R]]); Fg <- ifelse(is.na(b[[F]]), Inf, b[[F]])
  factor(ifelse(gen >= Fg, "3 fall apart", ifelse(gen >= Rg, "2 plateau", "1 reorganize")),
         c("1 reorganize", "2 plateau", "3 fall apart"))
}
by_phase <- function(d, scn, rep, gen, ..., R = "R_gain", F = "F_hzn")
  d |> mutate(phase = phase_of({{ scn }}, {{ rep }}, {{ gen }}, R, F))

cor_ph <- rean$corridor$d |> filter(Treatment != "Isotropic") |>
  mutate(phase = phase_of(Treatment, Replicate, Generation))
tru_ph <- rean$truth$d |> filter(Treatment != "Isotropic") |> mutate(phase = phase_of(Treatment, Replicate, Generation))
ss_ph  <- rean$sstest$tests |> filter(null == "N1w", scenario %in% SC[-1]) |>
  mutate(phase = phase_of(scenario, Replicate, generation))
ib_ph  <- ef$traj |> filter(status == "ok", Treatment %in% SC[-1]) |>
  mutate(Treatment = as.character(Treatment), phase = phase_of(Treatment, Replicate, generation),
         det = r2_pgd_dir > r2_cgd)
save(bounds, file = "data/derived/phase_boundaries.rda")

if (sys.nframe() == 0L) {
  options(width = 170)
  q <- function(x) if (all(is.na(x))) "—" else sprintf("%d/50: %s (%s–%s)", sum(!is.na(x)), median(x, na.rm = TRUE),
                                           quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
  cat("=== Boundary candidates: lineages with the marker / median (IQR) generation ===\n")
  print(bounds |> group_by(Scenario) |> summarise(across(c(R_gain, R_dR2, F_hzn, F_hzn_s, F_gain, F_diam, F_abs, F_nm), q)) |>
          as.data.frame())
  cat("\nReorganized (gain rise beyond isotropic 95th pct):\n"); print(table(bounds$Scenario, bounds$reorg))
  cat("\nAgreement between fall-apart markers (Redistributed; Spearman, pairwise complete; median difference vs F_hzn):\n")
  rb <- filter(bounds, Scenario == "Redistributed")
  mk <- c("F_hzn", "F_hzn_s", "F_gain", "F_diam", "F_abs", "F_nm")
  print(round(cor(rb[, mk], use = "pairwise.complete.obs", method = "spearman"), 2))
  print(sapply(mk, function(m) median(rb[[m]] - rb$F_hzn, na.rm = TRUE)))
  cat("\nPlateau length (F_hzn - R_gain), generations:\n")
  print(bounds |> filter(Scenario != "Isotropic") |> group_by(Scenario) |>
          summarise(n_both = sum(!is.na(R_gain) & !is.na(F_hzn)), med = median(F_hzn - R_gain, na.rm = TRUE),
                    neg = sum(F_hzn < R_gain, na.rm = TRUE)) |> as.data.frame())

  cat("\n=== Corridor direction by phase: fraction of links larger in the forward direction ===\n")
  print(cor_ph |> group_by(Treatment, phase) |>
          summarise(censuses = n(), w_all = mean(w_fwd_gt, na.rm = TRUE),
                    nm_censuses = sum(nm_ok), nm = mean(nm_fwd_gt, na.rm = TRUE),
                    w_same_censuses = mean(w_fwd_gt[nm_ok], na.rm = TRUE),
                    nm_minus_w_paired = mean((nm_fwd_gt - w_fwd_gt)[nm_ok], na.rm = TRUE), .groups = "drop") |>
          as.data.frame(), digits = 3)
  cat("Per-lineage (plateau only): lineages where Nm > w, mean over the lineage's plateau censuses\n")
  print(cor_ph |> filter(phase == "2 plateau", nm_ok) |> group_by(Treatment, Replicate) |>
          summarise(d = mean(nm_fwd_gt - w_fwd_gt), .groups = "drop") |> group_by(Treatment) |>
          summarise(lineages = n(), nm_better = sum(d > 0), sign_p = binom.test(sum(d > 0), n())$p.value) |> as.data.frame())
  cat("Isotropic (whole forward phase): w", round(mean(rean$corridor$d$w_fwd_gt[rean$corridor$d$Treatment == "Isotropic"]), 3),
      " Nm", round(mean(rean$corridor$d$nm_fwd_gt[rean$corridor$d$Treatment == "Isotropic"], na.rm = TRUE), 3), "\n")

  cat("\n=== Other inferences by phase ===\n")
  print(tru_ph |> group_by(Treatment, phase) |>
          summarise(censuses = n(), nm_retained = mean(Exclusion == ""), w_rho = median(w_rho_tr),
                    nm_rho = median(nm_rho_run, na.rm = TRUE), .groups = "drop") |> as.data.frame(), digits = 3)
  print(ss_ph |> group_by(scenario, phase) |>
          summarise(censuses = n(), ss_reject = mean(p < .05), ss_up = mean(r_S < 0), .groups = "drop") |> as.data.frame(), digits = 3)
  print(ib_ph |> group_by(Treatment, phase) |>
          summarise(censuses = n(), dR2_detect = mean(det), ratio_lt1 = mean(rel_fwd[det] < 1), .groups = "drop") |>
          as.data.frame(), digits = 3)

  cat("\n=== Robustness: corridor direction in the plateau under alternative fall-apart markers ===\n")
  for (F in c("F_hzn_s", "F_gain", "F_diam", "F_nm")) {
    z <- rean$corridor$d |> filter(Treatment != "Isotropic") |> mutate(phase = phase_of(Treatment, Replicate, Generation, F = F)) |>
      filter(phase == "2 plateau") |> group_by(Treatment) |>
      summarise(censuses = n(), w = mean(w_fwd_gt, na.rm = TRUE), nm = mean(nm_fwd_gt, na.rm = TRUE), .groups = "drop")
    cat(F, ":\n"); print(as.data.frame(z), digits = 3)
  }
}
