# R/figures/fig_gravity_field_bandwidth.R
#
# Fig_GravityFieldBandwidth: one census drawn as a gravity field at three bandwidths, to
# show what gamma does. The census is the one shipped with gstudio as `data(gravity)`
# (Redistributed scenario, lineage 4, generation 2504: the plateau-core census with the
# strongest source-sink gradient detected by both tests); its saved Population Graph is
# identical to popgraph(to_mv(gravity), gravity$Population).
#
# Panels: gamma = Inf (flat kernel, w = 1/k: S is degree alone), 1 (local mean) and
# 1/2 (degree-neutral, the manuscript's bandwidth); v2 also had 1/4, left to the sweep
# figure (Fig_BandwidthGammaSweep). Every panel uses
# the same node positions: gstudio's Kamada-Kawai layout (layout = "kk", computed on the
# pGD graph of the first panel; v1 placed nodes by deme index instead).
# Surface = interpolated S (red = source, blue = sink); arrows point
# source -> sink on each edge, length proportional to |Delta|, the weakest half hidden.
# Facet labels give each bandwidth's degree share (R^2 of Delta on the degree-only
# Delta0), cor(S, k) and r(S, x) (Spearman, S against deme number; negative = sources
# upstream, the test statistic of the source-sink gradient test).
#
# Run from the repository root:  Rscript R/figures/fig_gravity_field_bandwidth.R
# v4: node labels placed with ggrepel (gstudio's fixed labels overlapped in the dense
# upstream cluster).
# Output: media/fig-gravity-field-bandwidth-v4.png

suppressPackageStartupMessages({ library(igraph); library(gstudio); library(ggplot2) })
source("R/common/gravity_convention.R")
FIG <- "media/fig-gravity-field-bandwidth-v4.png"

g <- gc_load_graph(4, 2504, 2L)
V(g)$name <- as.character(gc_x(V(g)$name))                     # label nodes by deme number
PAN <- c(`γ = ∞ (degree only)` = Inf, `γ = 1 (local mean)` = 1, `γ = ½ (degree-neutral)` = 0.5)
f0 <- gravity_field(setNames(rep(list(g), length(PAN)), names(PAN)), gamma = PAN)
rx <- vapply(names(PAN), function(p) { n <- f0$nodes[f0$nodes$panel == p, ]
  cor(n$S, as.numeric(n$Stratum), method = "spearman") }, numeric(1))
lab <- sprintf("%s\ndegree share %.2f, cor(S, k) %.2f, r(S, x) %.2f", names(PAN),
               f0$diagnostics$degree_share, f0$diagnostics$cor_S_degree, rx)
f <- gravity_field(setNames(rep(list(g), length(PAN)), lab), gamma = PAN)
print(f)

p <- plot(f, layout = "kk", node_size = 4, node_labels = "none", min_delta = 0.5, subtitle_stats = FALSE) +
  facet_wrap(~panel, nrow = 1)
# labels: gstudio draws them as fixed text above each node, which collides in the dense
# upstream cluster; repel them instead, from the plot's own node layer (same positions)
nd <- Filter(function(l) inherits(l$geom, "GeomPoint") && is.data.frame(l$data) && "node" %in% names(l$data), p$layers)[[1]]$data
p <- p + ggrepel::geom_text_repel(data = nd, aes(x, y, label = node), size = 2.7, fontface = "bold", colour = "grey10",
                                  point.size = 4, box.padding = 0.4, point.padding = 0.15, min.segment.length = 0.3,
                                  segment.colour = "grey30", segment.size = 0.3, max.overlaps = Inf, seed = 1,
                                  bg.colour = "white", bg.r = 0.12)
if (file.exists(FIG)) cat(FIG, "exists; not overwritten\n") else {
  ggsave(FIG, p, width = 12, height = 7.5, dpi = 200, bg = "white"); cat("Wrote", FIG, "\n") }
