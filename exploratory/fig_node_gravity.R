# exploratory/fig_node_gravity.R
# Phase strips use the heterozygosity-gradient reorganization boundary (2026-09-29); earlier file = gain-based boundary.
# 2026-10-01: strips use the phase policy (unflagged lineages held in ascent; lineage_phases.R phase_policy_R_gen).
#
# Fig_SourceSinkScore: the population-level source-sink score S through time.
#   Top row: median S across the 50 replicates, by deme (y) and generation (x),
#            for Isotropic, Redistributed and Obstructed (columns). The shared
#            burn-in (1904-1999) is drawn in every column, left of the onset line.
#            Red = source (high S), blue = sink (low S), grey = 0.
#   Bottom row: the source-sink gradient test at each census (50 replicates):
#            fraction rejected at alpha = 0.05 (locked spec: Spearman r(S, x),
#            spectral null on the census's own graph, B = 999) and fraction with
#            sources upstream (r(S, x) < 0), with binomial 95% intervals.
#
# S is computed from the saved graphs with gstudio::source_sink_scores() (every 5
# generations in 1904-1999, every 10 generations in 2004-2994). With GG_TAG=g05 (gamma = 1/2,
# the manuscript's bandwidth) both inputs are written by exploratory/half_bandwidth_results.R:
# node S (data/derived/node_gravity_heatmap_g05.rds) and source_sink_test() results
# (data/derived/sstest_g05.rds; censuses every 5 generations in burn-in, every 50 in the
# forward phase). Without a tag (gamma = 1, earlier figure versions) the test results are
# rean$sstest in data/derived/node_potential_reanalysis.rda.
#
# Run from the repository root:  Rscript exploratory/fig_node_gravity.R
# Outputs: data/derived/node_gravity_heatmap.rds (per-census node S),
#          media/fig-source-sink-score.png

suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(tidyr)
  library(parallel); library(ggplot2); library(patchwork) })
source("exploratory/gravity_convention.R")
source("exploratory/lineage_phases.R")
TAG <- Sys.getenv("GG_TAG"); SFX <- if (nzchar(TAG)) paste0("_", TAG) else ""   # GG_TAG=g05: gamma = 1/2 inputs (exploratory/half_bandwidth_results.R), next figure version
CACHE <- sprintf("data/derived/node_gravity_heatmap%s.rds", SFX)
FIG   <- if (nzchar(TAG)) "media/fig-source-sink-score-v6.png" else "media/fig-source-sink-score-v4.png"   # v6: gstudio 1.15 source_sink_test() (identical with rounding or tolerance ties; v5: earlier draws); v2 adds the phase strip
SCEN  <- c(Isotropic = 1L, Redistributed = 2L, Obstructed = 3L)

if (file.exists(CACHE)) nodes <- readRDS(CACHE) else {
  jobs <- rbind(
    expand.grid(Replicate = 1:50, scenario = "Burn-in", generation = seq(1904L, 1999L, 5L), stringsAsFactors = FALSE),
    expand.grid(Replicate = 1:50, scenario = names(SCEN), generation = seq(2004L, 2994L, 10L), stringsAsFactors = FALSE))
  res <- mclapply(seq_len(nrow(jobs)), function(i) {
    j <- jobs[i, ]
    g <- gc_load_graph(j$Replicate, j$generation, if (j$scenario == "Burn-in") NULL else SCEN[[j$scenario]])
    if (is.null(g)) return(NULL)
    sc <- tryCatch(source_sink_scores(g, gamma = if (nzchar(TAG)) 0.5 else 1), error = function(e) NULL)
    if (is.null(sc)) return(NULL)
    data.frame(j, deme = gc_x(sc$Stratum), S = sc$S, k = sc$degree)
  }, mc.cores = max(1L, detectCores() - 1L))
  nodes <- do.call(rbind, res)
  attr(nodes, "missing") <- nrow(jobs) - sum(!vapply(res, is.null, logical(1)))
  saveRDS(nodes, CACHE)
}
cat("Censuses without a graph:", attr(nodes, "missing"), "\n")

## ---- heatmap data: burn-in copied into each scenario column -------------------------
bi  <- filter(nodes, scenario == "Burn-in")
hm  <- bind_rows(lapply(names(SCEN), function(s) mutate(bi, scenario = s)), filter(nodes, scenario != "Burn-in")) |>
  group_by(scenario, generation, deme) |>
  summarise(S = median(S), n = n(), .groups = "drop") |>
  mutate(width = ifelse(generation < 2000, 5, 10), scenario = factor(scenario, names(SCEN)))
lim <- quantile(abs(hm$S), 0.98)
cat(sprintf("Colour limit +/- %.3f (98th percentile of |median S|); cells beyond it are squished\n", lim))
print(hm |> filter(generation %in% c(1999, 2454, 2824, 2994)) |> group_by(scenario, generation) |>
        summarise(S_deme1 = S[deme == 1], S_deme13 = S[deme == 13], S_deme25 = S[deme == 25],
                  r = cor(S, deme, method = "spearman"), .groups = "drop") |> as.data.frame(), digits = 3)

## ---- test data ------------------------------------------------------------------------
tt <- if (nzchar(TAG)) readRDS(sprintf("data/derived/sstest%s.rds", SFX)) else {
  e <- new.env(); load("data/derived/node_potential_reanalysis.rda", envir = e); e$rean$sstest$tests |> filter(null == "N1w") }
tt <- bind_rows(lapply(names(SCEN), function(s) filter(tt, scenario == "Burn-in") |> mutate(scenario = s)),
                filter(tt, scenario %in% names(SCEN)))
ci <- function(x) { b <- binom.test(sum(x), length(x))$conf.int; c(b[1], b[2]) }
pw <- tt |> group_by(scenario, generation) |>
  summarise(Rejected = mean(p < 0.05), rej_lo = ci(p < 0.05)[1], rej_hi = ci(p < 0.05)[2],
            `Sources upstream` = mean(r_S < 0), up_lo = ci(r_S < 0)[1], up_hi = ci(r_S < 0)[2], .groups = "drop")
# burn-in censuses are every 5 generations: pool them into 25-generation bins so
# every plotted point rests on a comparable number of censuses
pw_bi <- tt |> filter(generation < 2000) |> mutate(generation = 1900 + 25 * ((generation - 1900) %/% 25) + 12) |>
  group_by(scenario, generation) |>
  summarise(Rejected = mean(p < 0.05), rej_lo = ci(p < 0.05)[1], rej_hi = ci(p < 0.05)[2],
            `Sources upstream` = mean(r_S < 0), up_lo = ci(r_S < 0)[1], up_hi = ci(r_S < 0)[2], .groups = "drop")
pw <- bind_rows(pw_bi, filter(pw, generation > 2000))
pl <- bind_rows(
  transmute(pw, scenario, generation, series = "Rejected (α = 0.05)", value = Rejected, lo = rej_lo, hi = rej_hi),
  transmute(pw, scenario, generation, series = "Sources upstream, r(S, x) < 0", value = `Sources upstream`, lo = up_lo, hi = up_hi)) |>
  mutate(scenario = factor(scenario, names(SCEN)),
         series = factor(series, c("Rejected (α = 0.05)", "Sources upstream, r(S, x) < 0")))

## ---- figure -------------------------------------------------------------------------------
ink <- "#3d3d3a"; grid <- "#e6e5e1"
th <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid, linewidth = 0.3),
        axis.text = element_text(colour = ink), axis.title = element_text(colour = ink),
        strip.text = element_text(colour = ink, face = "bold", size = 10),
        legend.title = element_text(colour = ink, size = 9), legend.text = element_text(colour = ink, size = 9))
xs <- scale_x_continuous(limits = c(1900, 3000), breaks = seq(2000, 2750, 250), expand = c(0, 0))

top <- ggplot(hm, aes(generation, deme, fill = S, width = width)) +
  geom_tile(height = 1) +
  geom_vline(xintercept = 2000, colour = ink, linewidth = 0.4, linetype = "22") +
  scale_fill_gradient2(low = "#1c5cab", mid = "#f0efec", high = "#b3261e", midpoint = 0,
                       limits = c(-lim, lim), oob = scales::squish,
                       breaks = c(-lim, 0, lim), labels = c("sink", "0", "source"),
                       name = "Median S") +
  scale_y_continuous(breaks = c(1, 5, 10, 15, 20, 25), expand = c(0, 0)) + xs +
  facet_grid(. ~ scenario, labeller = as_labeller(c(Isotropic = "Symmetric", Redistributed = "Redistributed", Obstructed = "Obstructed"))) +
  labs(x = NULL, y = "Deme") + th +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), legend.position = "right",
        legend.key.height = unit(1.2, "cm"), legend.key.width = unit(0.35, "cm"), panel.spacing.x = unit(0.8, "lines"))

cols <- c("#1baf7a", "#eb6834")   # aqua, orange: blue is reserved for the phase strip
bot <- ggplot(pl, aes(generation, value, colour = series, fill = series)) +
  geom_hline(yintercept = 0.05, colour = "#9c9a92", linewidth = 0.3, linetype = "22") +
  geom_hline(yintercept = 0.5, colour = "#9c9a92", linewidth = 0.3, linetype = "13") +
  geom_vline(xintercept = 2000, colour = ink, linewidth = 0.4, linetype = "22") +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent, expand = c(0, 0.01)) + xs +
  facet_grid(. ~ scenario, labeller = as_labeller(c(Isotropic = "Symmetric", Redistributed = "Redistributed", Obstructed = "Obstructed"))) +
  labs(x = NULL, y = "Share of censuses") + th +
  theme(strip.text = element_blank(), legend.position = "right", panel.spacing.x = unit(0.8, "lines"),
        axis.text.x = element_blank())

strip <- lp_strip(names(SCEN), xs, facet = "col", base_size = 10) +
  geom_vline(xintercept = 2000, colour = ink, linewidth = 0.4, linetype = "22") +
  labs(x = "Generation") + theme(panel.spacing.x = unit(0.8, "lines"), legend.position = "right")

p <- top / bot / strip + plot_layout(heights = c(2.2, 1, 0.45))
ggsave(FIG, p, width = 11, height = 7.0, dpi = 200, bg = "white")
cat("Saved", FIG, "\n")
