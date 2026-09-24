# Fig_IsotropicBaseline and Tab_IsotropicBaseline.
#
# Structural stability of the Population Graphs under symmetric migration
# (N=100) through the full simulated timeline -- the 2,000-generation burn-in
# (gens 0-1999) followed by the 1,000-generation isotropic forward phase
# (scenario 1, gens 2004-2999) -- for three statistics: average standardized
# degree centrality (C_D-bar; Dyer 2007), graph diameter, and the graph-mean
# absolute asymmetry index (mean |Delta_ij|). The signed graph mean
# (Delta-bar) is summarized for the text only; under symmetric migration the
# positive and negative edges cancel and it stays centred on zero.
#
# Each statistic is fit, per replicate, over the burn-in only with an offset
# exponential approach to equilibrium,
#     y(t) = y_inf + (y_0 - y_inf) * exp(-t / tau),
# which covers C_D-bar decaying from above and diameter / mean |Delta| rising
# from below. Replicates are the unit of replication (consecutive snapshots
# are strongly autocorrelated), so parameters are summarized as the median and
# 2.5-97.5% range across the per-replicate fits; t95 = 3 tau is the time to
# reach 95% of the approach. A stretched exponential (exp(-(t/tau)^beta)) was
# checked and not adopted: beta ~ 1.06 for diameter and ~ 0.93 for mean
# |Delta| (no AIC gain), and ~ 0.77 for C_D-bar (modest gain), so the single
# exponential is kept for all three to keep tau comparable across statistics.
# Stationarity of the forward phase is assessed as each replicate's forward
# mean relative to its own fitted plateau.
#
# Source data (data/isotropic_baseline.csv) is the N=100 burn-in and scenario-1
# subset of the private research repo's data/graph_summary.rda
# (R/graph_summarizer.R), with mean_delta and mean_abs_delta recomputed from the
# stored graphs with gstudio::graph_asymmetries(). Diameter is the weighted
# (conditional genetic distance) igraph diameter. C_D-bar is derived from the
# edge count, since the mean of d(v_i)/(K-1) over K nodes equals 2|E|/(K(K-1)).
#
# Run from the repository root:
#   Rscript R/isotropic_baseline.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

iso <- read.csv("data/isotropic_baseline.csv", stringsAsFactors = FALSE)

K <- 25
BURNIN_END <- 2000

stat_levels <- c("CD", "diameter", "mean_abs_delta")
stat_labels <- c("bar(C)[D]", "Diameter", "bar('|' * Delta * '|')")

long <- iso |>
  mutate(CD = 2 * n_edges / (K * (K - 1))) |>
  pivot_longer(all_of(stat_levels), names_to = "stat", values_to = "value")

## ---- Per-replicate exponential fits (burn-in only) ------------------------

fit_exp <- function(t, y) {
  tryCatch(
    nls(y ~ yinf + (y0 - yinf) * exp(-t / tau),
        start     = list(yinf = mean(y[t >= 1500]), y0 = y[1], tau = 300),
        algorithm = "port", lower = c(-Inf, -Inf, 1), upper = c(Inf, Inf, 1e5),
        control   = nls.control(maxiter = 500)),
    error = function(e) NULL)
}

fits <- long |>
  filter(phase == "burnin") |>
  group_by(stat, replicate) |>
  group_modify(function(df, key) {
    m <- fit_exp(df$generation, df$value)
    if (is.null(m)) return(tibble(ok = FALSE))
    cf <- coef(m)
    tibble(ok = TRUE, y0 = cf[["y0"]], yinf = cf[["yinf"]], tau = cf[["tau"]],
           R2 = 1 - sum(resid(m)^2) / sum((df$value - mean(df$value))^2))
  }) |>
  ungroup()

n_fail <- sum(!fits$ok)
fits   <- filter(fits, ok)

# Forward-phase stationarity: each replicate's forward mean vs its own plateau.
fwd <- long |>
  filter(phase == "forward") |>
  group_by(stat, replicate) |>
  summarise(fwd_mean = mean(value), .groups = "drop") |>
  inner_join(fits, by = c("stat", "replicate")) |>
  mutate(rel_dev = (fwd_mean - yinf) / yinf)

## ---- Tab_IsotropicBaseline ------------------------------------------------

q <- function(x, p) unname(quantile(x, p))
fmt <- function(x, d) sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"),
                             median(x), q(x, 0.025), q(x, 0.975))

tbl_iso <- fits |>
  left_join(select(fwd, stat, replicate, rel_dev), by = c("stat", "replicate")) |>
  group_by(stat) |>
  summarise(
    `Initial (y0)`             = fmt(y0,   if (first(stat) == "diameter") 2 else 3),
    `Plateau (y_inf)`          = fmt(yinf, if (first(stat) == "diameter") 2 else 3),
    `tau (gens)`               = fmt(tau, 0),
    `t95 = 3 tau (gens)`       = fmt(3 * tau, 0),
    `R2`                       = sprintf("%.2f", median(R2)),
    `Forward dev. from plateau (%)` = fmt(100 * rel_dev, 1),
    .groups = "drop") |>
  mutate(stat = factor(stat, stat_levels,
                       c("C_D-bar", "Diameter", "Mean |Delta|"))) |>
  arrange(stat) |>
  rename(Statistic = stat)

cat("=== Tab_IsotropicBaseline ===\n")
cat("Median (2.5-97.5% across replicates) of per-replicate burn-in fits\n")
print(as.data.frame(tbl_iso), row.names = FALSE)
cat(sprintf("Fits failing to converge: %d of %d\n",
            n_fail, n_fail + nrow(fits)))

## ---- Inline text values ---------------------------------------------------

# Signed Delta-bar: per-replicate window means, t-interval across replicates.
signed_window <- function(lo, hi) {
  r  <- iso |> filter(generation >= lo, generation <= hi) |>
    group_by(replicate) |> summarise(m = mean(mean_delta), .groups = "drop")
  tt <- t.test(r$m)
  c(mean = mean(r$m), lo = tt$conf.int[1], hi = tt$conf.int[2])
}
rv_iso_signed_all <- signed_window(0, 2999)
rv_iso_signed_fwd <- signed_window(BURNIN_END, 2999)
rv_iso_signed_max <- max(abs(tapply(iso$mean_delta, iso$generation, mean)))

cat(sprintf(paste0(
  "\nrv_iso_signed_all = %.5f (95%% CI %.5f to %.5f; gens 0-2999)",
  "\nrv_iso_signed_fwd = %.5f (95%% CI %.5f to %.5f; gens 2000-2999)",
  "\nrv_iso_signed_max = %.5f (largest |across-replicate mean| at any generation)\n"),
  rv_iso_signed_all[1], rv_iso_signed_all[2], rv_iso_signed_all[3],
  rv_iso_signed_fwd[1], rv_iso_signed_fwd[2], rv_iso_signed_fwd[3],
  rv_iso_signed_max))

## ---- Fig_IsotropicBaseline ------------------------------------------------

traj <- long |>
  group_by(stat, generation) |>
  summarise(mean = mean(value), sd = sd(value), .groups = "drop") |>
  mutate(Statistic = factor(stat, stat_levels, stat_labels))

# Fitted curve from the median per-replicate parameters: solid over the
# burn-in it was fit to, dashed where extrapolated across the forward phase.
med_par <- fits |> group_by(stat) |>
  summarise(y0 = median(y0), yinf = median(yinf), tau = median(tau), .groups = "drop")
# Reused as the isotropic reference curve by R/forward_scenarios.R.
write.csv(med_par, "data/derived/isotropic_fit_params.csv", row.names = FALSE)
curve <-expand_grid(stat = stat_levels, generation = 0:2999) |>
  left_join(med_par, by = "stat") |>
  mutate(fit = yinf + (y0 - yinf) * exp(-generation / tau),
         segment = ifelse(generation < BURNIN_END, "fit", "extrapolated"),
         Statistic = factor(stat, stat_levels, stat_labels))

fig_iso <- ggplot(traj, aes(generation, mean)) +
  geom_vline(xintercept = BURNIN_END, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "grey70", alpha = 0.5) +
  geom_line(data = curve, aes(generation, fit, linetype = segment),
            colour = "firebrick", linewidth = 0.7, show.legend = FALSE) +
  scale_linetype_manual(values = c(fit = "solid", extrapolated = "22")) +
  facet_grid(Statistic ~ ., scales = "free_y", labeller = label_parsed) +
  labs(x = "Generation", y = NULL) +
  theme_minimal(base_size = 11)

out_png <- "data/derived/fig-isotropic-baseline_reproduced.png"
ggsave(out_png, fig_iso, width = 7, height = 6.5, dpi = 150)
cat(sprintf("\nSaved Fig_IsotropicBaseline to %s\n", out_png))
