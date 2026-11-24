# Figure S6 (Sampling intensity) -> media/fig-sampling-intensity-v3.png
# (v3: regular x-axis breaks 5-25 by 5, no vertical grid; v2: one merged legend for n; v1 had duplicated n legends)
#
# Power and false-positive rates against the number of populations sampled (K), one line
# per number of individuals per population (n), 95% lineage-bootstrap intervals.
#   Top left:     Redistributed plateau power, source-sink gradient test (alpha = 0.05)
#   Top right:    Redistributed plateau power, directional isolation at the held-out (K, n) cutoff
#   Bottom left:  source-sink false-positive rate, Symmetric (solid, band) and sym-low (dashed)
#   Bottom right: false-positive rate of the uncalibrated rule Delta R^2 > 0, same encoding
# Input: data/derived/density_power_tables.rda (R/analysis/density_power_results.R; tables
# `ss` and `dr`). Run from the repository root.
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })
load("data/derived/density_power_tables.rda")
ink <- "#3d3d3a"
ncol_ <- c(`10` = "#86b6ef", `20` = "#3987e5", `30` = "#256abf", `50` = "#184f95", `100` = "#0d366b")  # ordinal blue ramp (validated)
th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), axis.text = element_text(colour = ink), axis.title = element_text(colour = ink),
        plot.title = element_text(colour = ink, size = 11, face = "bold"), legend.text = element_text(colour = ink),
        legend.title = element_text(colour = ink))
prep <- function(d) d |> mutate(n = factor(n, levels = c(10, 20, 30, 50, 100)))
xs <- scale_x_continuous(breaks = seq(5, 25, 5), limits = c(5, 25))
cs <- list(scale_colour_manual(values = ncol_, name = "Individuals per population"),
           scale_fill_manual(values = ncol_, guide = "none"),
           guides(colour = guide_legend(order = 1, override.aes = list(linetype = "solid", fill = NA))))
power_panel <- function(d, title) ggplot(prep(d), aes(K, est, colour = n, fill = n, group = n)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.12) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  xs + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) + cs +
  labs(x = NULL, y = "Plateau power", title = title) + th
fpr_panel <- function(d, title, ymax) {
  d <- prep(d) |> mutate(sample = factor(as.character(sample), c("Symmetric", "sym-low")))
  ggplot(d, aes(K, est, colour = n, group = interaction(n, sample))) +
    geom_hline(yintercept = 0.05, colour = "grey55", linewidth = 0.4, linetype = "dotted") +
    geom_ribbon(data = filter(d, sample == "Symmetric"), aes(ymin = lo, ymax = hi, fill = n), colour = NA, alpha = 0.12) +
    geom_line(aes(linetype = sample), linewidth = 0.7) + geom_point(data = filter(d, sample == "Symmetric"), size = 1.4) +
    scale_linetype_manual(values = c(Symmetric = "solid", `sym-low` = "22"), name = "Migration", guide = guide_legend(order = 2, override.aes = list(colour = ink)),
                          labels = c(Symmetric = "Symmetric", `sym-low` = "Strong drift (sym-low)")) +
    xs + scale_y_continuous(limits = c(0, ymax)) + cs +
    labs(x = "Populations sampled (K)", y = "False-positive rate", title = title) + th
}
a <- power_panel(filter(ss, sample == "Redistributed"), "Source–sink gradient test")
b <- power_panel(filter(dr, rule == "det_cell", sample == "Redistributed"), "Directional isolation, calibrated cutoff")
c_ <- fpr_panel(filter(ss, sample %in% c("Symmetric", "sym-low")), "Source–sink gradient test", 0.2)
d_ <- fpr_panel(filter(dr, rule == "det0", sample %in% c("Symmetric", "sym-low")), "Directional isolation, rule ΔR² > 0", 0.55)
p <- (a + b) / (c_ + d_) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.box = "horizontal", plot.tag = element_text(colour = ink, face = "bold"))
FIG <- "media/fig-sampling-intensity-v3.png"
stopifnot(!file.exists(FIG))                                   # new filename per variant; never overwrite
ggsave(FIG, p, width = 9, height = 8, dpi = 150, bg = "white")
cat("Wrote", FIG, "\n")
