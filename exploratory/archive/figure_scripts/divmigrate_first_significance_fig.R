# divmigrate_first_significance_fig.R
# [dropped thread] time to first significance, 50 reps (OUR BH z-test) -> media/fig-divmigrate-first-significance.png
# Run from the manuscript repo root.
# Recovered from the 2026-09-24..27 working session scratchpad.

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
x <- read.csv("data/time_to_first_significance.csv") |>
  mutate(Treatment = factor(Treatment, c("Isotropic", "Redistributed", "Obstructed")),
         dist = abs(as.integer(sub("Pop", "", from)) - as.integer(sub("Pop", "", to))))
cols <- c(Isotropic = "grey45", Redistributed = "firebrick", Obstructed = "steelblue")
th <- theme_minimal(base_size = 11) + theme(legend.position = "none")
cnt <- x |> count(Treatment, first_sig_gen)
pA <- ggplot(cnt, aes(factor(first_sig_gen), n, fill = Treatment)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.3, size = 3.3) +
  scale_fill_manual(values = cols) + scale_y_continuous(limits = c(0, 55), breaks = c(0, 10, 25, 50)) +
  scale_x_discrete(labels = c("1999" = "1999\n(before treatment)", "2004" = "2004\n(first forward census)")) +
  labs(x = "Generation of first significant asymmetry", y = "Replicates (of 50)", fill = NULL, title = "A  Time to first significance") +
  th + theme(legend.position = "top")
pB <- ggplot(x, aes(Treatment, n_sig_pairs, colour = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.12, height = 0, size = 1.2, alpha = 0.7) +
  scale_colour_manual(values = cols) + scale_y_continuous(limits = c(0, 300)) +
  labs(x = NULL, y = "Significant pairs (of 300)", title = "B  Pairs flagged at that census") + th
dd <- x |> distinct(Replicate, first_sig_gen, from, to, dist) |> count(dist)
pC <- ggplot(dd, aes(dist, n)) + geom_col(fill = "grey40", width = 0.8) +
  geom_col(data = filter(dd, dist == 1), fill = "darkorange", width = 0.8) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 24)) +
  labs(x = "Chain distance of the top pair (1 = true stepping-stone link)", y = "Distinct results", title = "C  Where the top pair sits (each distinct result counted once; 51 in total)") + th
p <- (pA | pB) / pC + plot_layout(heights = c(1.1, 1)) +
  plot_annotation(caption = paste0("divMigrate (Nm) single-snapshot test: bootstrap of individuals (1,000 resamples), Benjamini–Hochberg across the 300 pairs at q = 0.05.\n",
                                   "Scan starts at generation 1999 (shared symmetric burn-in); 49 of 50 replicates are already significant there, so all three treatments share that result."))
ggsave("media/fig-divmigrate-first-significance.png", p, width = 10, height = 7.5, dpi = 150)
u <- distinct(x, Replicate, first_sig_gen, from, to, dist); cat("distinct results:", nrow(u), " top-pair distance median", median(u$dist), "; adjacent (d=1):", sum(u$dist == 1), "\n")
print(x |> group_by(Treatment) |> summarise(median_pairs = median(n_sig_pairs), range = paste(range(n_sig_pairs), collapse = "-")))
