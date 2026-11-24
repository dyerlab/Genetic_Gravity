# R/figures/fig_source_sink_examples.R
#
# Fig_SourceSinkExamples: four censuses drawn with plot(gravity_field(graph)) at the
# degree-neutral bandwidth (gamma = 1/2), as worked examples of what "no direction",
# "clear direction" and "direction collapsing downstream" look like on a map
# (local/wind_map_instructions.md, Part B; one bandwidth only, the bandwidth contrast is
# Fig_GravityFieldBandwidth).
#
# Selection rule: in each sample, the census whose r(S, x) at gamma = 1/2 is closest to
# that sample's median, ties broken by the lowest lineage number, then earliest generation.
#   Symmetric: all forward-phase test censuses (every 50 generations, 1,000 censuses).
#   Redistributed plateau, Obstructed plateau, Redistributed decay: the phase cores
#   (middle half of each lineage's phase, informative horizon 4.5; as
#   R/analysis/phase_core_rerun.R), i.e. the censuses behind the medians in the text.
# Exception (v3): the Redistributed plateau panel is the census distributed with gstudio as
#   data(gravity) and walked through in its vignette (lineage 4, generation 2504, plateau
#   core; chosen there as the strongest gradient detected by both directional tests), so the
#   figure shows the example readers can reproduce. v2 used the median census (lineage 31,
#   generation 2554, p = 0.175).
# r(S, x) and p come from the source-sink test results (data/derived/sstest_g05.rds, as
# written by R/analysis/half_bandwidth_results.R); r(S, x) is recomputed from the field
# as a check.
#
# Node positions: demes in order along a horseshoe (a 270-degree arc, deme 1 at lower left,
# deme 25 at lower right), the same in every panel, so long edges separate as chords; a
# gently curving line (the "river" of the instructions) bundled them. The Kamada-Kawai
# layout of Fig_GravityFieldBandwidth is not used here because it differs between censuses.
#
# Run from the repository root:  Rscript R/figures/fig_source_sink_examples.R
# Output: media/fig-source-sink-examples-v3.png, local/source_sink_examples.md (selection)

dir.create("output", showWarnings = FALSE)   # run reports (not tracked)
suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(ggplot2) })
source("R/common/gravity_convention.R")
FIG <- "media/fig-source-sink-examples-v3.png"
ld <- function(f) { e <- new.env(); load(f, envir = e); e }
SCEN <- c(Symmetric = 1L, Redistributed = 2L, Obstructed = 3L); AS <- c("Redistributed", "Obstructed")

## ---- phase cores (as R/analysis/phase_core_rerun.R, horizon 4.5) ---------------------------
IV <- ld("data/derived/phase_core_rerun.rda")$IV$F45
core_phase <- function(d) {
  x <- d |> left_join(IV |> rename(sample = Scenario), by = c("sample", "Replicate")); g <- x$generation
  inn <- function(a, b) { lo <- a + (b - a) / 4; hi <- b - (b - a) / 4; g >= lo & g <= hi }
  d$phase <- ifelse(!x$sample %in% AS, "forward phase",
             ifelse(inn(x$asc_a, x$asc_b) %in% TRUE, "ascent",
             ifelse(x$pl_ok & inn(x$pl_a, x$pl_b) %in% TRUE, "plateau",
             ifelse(x$de_ok & inn(x$de_a, x$de_b) %in% TRUE, "decay", NA))))
  d }

## ---- selection --------------------------------------------------------------------------
tt <- readRDS("data/derived/sstest_g05.rds") |> mutate(sample = recode(scenario, Isotropic = "Symmetric")) |>
  filter(sample %in% names(SCEN)) |> core_phase()
WANT <- c("Symmetric forward phase", "Redistributed plateau", "Obstructed plateau", "Redistributed decay")
sel <- tt |> mutate(key = paste(sample, phase)) |> filter(key %in% WANT) |> group_by(key) |>
  mutate(n = n(), med = median(r_S), dist = round(abs(r_S - med), 9), pct = 100 * (rank(r_S) - 0.5) / n()) |>
  arrange(dist, Replicate, generation, .by_group = TRUE) |>
  filter(if (key[1] == "Redistributed plateau") Replicate == 4 & generation == 2504 else row_number() == 1) |>   # data(gravity)
  ungroup() |> mutate(key = factor(key, WANT)) |> arrange(key)
stopifnot(nrow(sel) == 4)
print(sel |> select(key, Replicate, generation, r_S, med, pct, p, n), width = 200)

## ---- fields -------------------------------------------------------------------------------
gs <- lapply(seq_len(nrow(sel)), function(i) {
  g <- gc_load_graph(sel$Replicate[i], sel$generation[i], SCEN[[sel$sample[i]]])
  V(g)$name <- as.character(gc_x(V(g)$name)); g })
f0 <- gravity_field(setNames(gs, WANT), gamma = 0.5)
rx <- vapply(WANT, function(p) { n <- f0$nodes[f0$nodes$panel == p, ]
  cor(signif(n$S, 9), as.numeric(n$Stratum), method = "spearman") }, numeric(1))   # signif: near-ties are ties, as in source_sink_test()
stopifnot(max(abs(rx - sel$r_S)) < 1e-6)                        # field agrees with the test's r(S, x)
lab <- sprintf("%s: lineage %d, generation %d\nr(S, x) %.2f (median %.2f), p = %.3f, degree share %.2f",
               c("Symmetric", "Redistributed, plateau", "Obstructed, plateau", "Redistributed, decay"),
               sel$Replicate, sel$generation, sel$r_S, sel$med, sel$p, f0$diagnostics$degree_share)
f <- gravity_field(setNames(gs, lab), gamma = 0.5)

th <- seq(225, -45, length.out = 25) * pi / 180
xy <- cbind(x = 10 * cos(th), y = 10 * sin(th)); rownames(xy) <- as.character(1:25)
p <- plot(f, layout = xy, node_size = 4, node_labels = "name", min_delta = 0.5, subtitle_stats = FALSE) +
  facet_wrap(~panel, ncol = 2) +
  geom_text(data = data.frame(x = xy[c(1, 25), 1], y = xy[c(1, 25), 2] - 2.2, label = c("upstream", "downstream")),
            aes(x, y, label = label), inherit.aes = FALSE, size = 3, colour = "grey35")
if (file.exists(FIG)) cat(FIG, "exists; not overwritten\n") else {
  ggsave(FIG, p, width = 10.5, height = 10, dpi = 200, bg = "white"); cat("Wrote", FIG, "\n") }

md <- c("# Fig_SourceSinkExamples: selected censuses", "",
        "Written by `R/figures/fig_source_sink_examples.R`. Rule: r(S, x) at γ = ½ closest to the sample median (ties: lowest lineage, then earliest generation), except the *Redistributed* plateau panel, which is gstudio's `data(gravity)` census (lineage 4, generation 2504). Percentile = rank of the chosen census's r(S, x) within its sample.", "",
        "| Sample | n | Median r(S, x) | Lineage | Generation | r(S, x) | Percentile | p | Degree share |",
        "| :-- | --: | --: | --: | --: | --: | --: | --: | --: |",
        sprintf("| %s | %d | %.3f | %d | %d | %.3f | %.0f | %.3f | %.3f |", sel$key, sel$n, sel$med, sel$Replicate,
                sel$generation, sel$r_S, sel$pct, sel$p, f0$diagnostics$degree_share))
writeLines(md, "output/source_sink_examples.md")
