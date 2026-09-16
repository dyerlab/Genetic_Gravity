# Fig_StreamCScissors and Tab_DifferentiationVsAsymmetry.
#
# Reproduces the Stream C Arm A specificity check: across a gradient of
# symmetric (1:1) migration scenarios at progressively lower migration rate
# (N=100; conds 1, 4, 5, 6 = Isotropic, sym-mid, sym-low, sym-verylow), does
# the asymmetry index stay pinned at the isotropic null while neutral
# differentiation (mean cGD) climbs? Source data (data/graph_summary.csv)
# is the per-replicate/scenario/generation graph summary already computed
# from the individual-based simulation output by the private research repo's
# R/graph_summarizer.R; no raw simulation data is required to run this
# script. Logic below ports the private repo's R/specificity_analysis.R
# armA_forward() near-verbatim, only changing the data source (CSV instead of
# .rda) -- the scenario column is already a plain integer after the CSV
# round-trip so the as.integer(as.character(...)) coercion used against the
# factor column in the private repo is unnecessary here (checked with str()).
#
# Run from the repository root:
#   Rscript R/streamc_specificity.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

graph_summary <- read.csv("data/graph_summary.csv", stringsAsFactors = FALSE)

ARM_A_CONDS  <- c(1L, 4L, 5L, 6L)
ARM_A_LABELS <- c("1" = "Isotropic", "4" = "sym-mid",
                   "5" = "sym-low", "6" = "sym-verylow")

## ---- armA_forward() -------------------------------------------------------

armA_forward <- function(gs, conds = ARM_A_CONDS) {
  if (!"N" %in% names(gs)) gs$N <- 100L
  fw <- gs |> filter(phase == "forward", N == 100L, scenario %in% conds)
  present <- sort(unique(fw$scenario))
  fw <- fw |>
    mutate(
      scenario  = factor(ARM_A_LABELS[as.character(scenario)],
                          levels = unname(ARM_A_LABELS[as.character(present)])),
      step      = generation - 2000,
      replicate = factor(replicate))
  attr(fw, "present") <- present
  fw
}

fw <- armA_forward(graph_summary)

sa <- list(
  traj = fw |>
    group_by(scenario, generation, step) |>
    summarise(mean_cGD   = mean(mean_cGD,   na.rm = TRUE),
              mean_delta = mean(mean_delta, na.rm = TRUE), .groups = "drop"),
  summary = fw |>
    filter(generation == 2999) |>
    group_by(scenario) |>
    summarise(mean_cGD     = round(mean(mean_cGD, na.rm = TRUE), 3),
              signed_delta = round(mean(mean_delta, na.rm = TRUE), 4),
              abs_delta    = round(mean(abs(mean_delta), na.rm = TRUE), 4),
              .groups = "drop"))

## ---- Tab_DifferentiationVsAsymmetry ---------------------------------------

tbl_streamc <- sa$summary |>
  transmute(Scenario = scenario, `Mean cGD` = mean_cGD,
            `Signed Delta-bar` = signed_delta, `|Delta-bar|` = abs_delta)

cat("=== Tab_DifferentiationVsAsymmetry ===\n")
print(as.data.frame(tbl_streamc), row.names = FALSE)

## ---- Inline text values ---------------------------------------------------

s    <- sa$summary
iso  <- s[s$scenario == "Isotropic", ]
ext  <- s[nrow(s), ]  # lowest-m gradient scenario
fold <- round(ext$mean_cGD / iso$mean_cGD, 1)
max_signed <- signif(max(abs(s$signed_delta), na.rm = TRUE), 1)

cat(sprintf(
  "\nExtreme gradient scenario = %s, mean cGD fold-increase vs isotropic = %sx\n",
  as.character(ext$scenario), fold))
cat(sprintf("max |signed Delta-bar| across gradient = %s\n", max_signed))
cat(sprintf("|Delta-bar|: isotropic = %s, %s = %s\n",
            iso$abs_delta, as.character(ext$scenario), ext$abs_delta))

## ---- Fig_StreamCScissors --------------------------------------------------

fig_scissors <- sa$traj |>
  pivot_longer(c(mean_cGD, mean_delta), names_to = "measure", values_to = "value") |>
  mutate(measure = recode(measure,
                           mean_cGD   = "Differentiation (mean cGD)",
                           mean_delta = "Asymmetry index (Delta-bar)")) |>
  ggplot(aes(step, value, colour = scenario)) +
  geom_hline(data = data.frame(measure = "Asymmetry index (Delta-bar)", y = 0),
             aes(yintercept = y), linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ measure, scales = "free_y") +
  labs(x = "Generations since onset", y = NULL, colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")

out_png <- "data/derived/fig-streamc-scissors_reproduced.png"
ggsave(out_png, fig_scissors, width = 9, height = 3.6, dpi = 150)
cat(sprintf("\nSaved reproduced Fig_StreamCScissors to %s\n", out_png))
