# exploratory/lineage_phases.R
#
# The one place a census is assigned to a phase. Each replicate lineage x
# treatment has its own boundaries (exploratory/phase_boundaries.R):
#   R_gen  end of reorganization: t90 of the lineage's smoothed within-deme
#          heterozygosity gradient along the chain (OLS slope of He on deme index,
#          interior demes), from its burn-in level to its maximum
#          (exploratory/he_gradient_boundary.R). Gravity-free, like F_gen.
#          (Until 2026-09-29 this was R_gain, t90 of the within-pGD directional
#          gain, a gravity statistic; kept as column R_gain for reference.)
#   F_gen  start of falling apart: first census at which the lineage's demes
#          carry fewer than 4.5 polymorphic loci on average (the informative
#          horizon, R/information_loss.R). NA = never crossed by generation 2999.
#          R_gen NA (the gain never rose; 1 Obstructed lineage) = still reorganizing.
# Phases of a forward census at generation t:
#   "Reorganize"  t <  R_gen
#   "Plateau"     R_gen <= t < F_gen
#   "Fall apart"  t >= F_gen
# Redistributed and Obstructed lineages get phases. Isotropic and the symmetric
# gradient get "Symmetric" (no imposed direction, hence no reorganization; F_hzn
# never fires in an isotropic lineage), and burn-in censuses get "Burn-in".
# `reorganized` flags lineages whose gain rise exceeds the 95th percentile of
# the isotropic rises (44/50 Redistributed, 39/50 Obstructed); lineages that did
# not reorganize still get phases from R_gen, and can be dropped with it.
#
# Usage:  source("exploratory/lineage_phases.R")
#         d$phase <- lp_phase(d$scenario, d$replicate, d$generation)
# Table:  data/derived/lineage_phases.csv (written by running this file)

LP_LEVELS <- c("Burn-in", "Symmetric", "Reorganize", "Plateau", "Fall apart")

lp_bounds <- local({
  f <- "data/derived/lineage_phases.csv"
  if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
})

# Phase policy (Rodney, 2026-09-30): lineages whose He-gradient rise does not
# clear the isotropic 95th percentile (reorganized == FALSE) are held in
# ascent until decay; column phase_policy_R_gen = R_gen if reorganized, else
# F_gen (NA if none). lp_phase() uses it by default; rcol = "R_gen" gives the
# earlier all-t90 policy.
lp_phase <- function(scenario, replicate, generation, bounds = lp_bounds, rcol = "phase_policy_R_gen") {
  scenario <- as.character(scenario)
  b <- bounds[match(paste(scenario, replicate), paste(bounds$Scenario, bounds$Replicate)), ]
  R <- ifelse(is.na(b[[rcol]]), Inf, b[[rcol]])    # no R: ascent never completed
  F <- ifelse(is.na(b$F_gen), Inf, b$F_gen)
  out <- ifelse(generation < 2000, "Burn-in",
         ifelse(!scenario %in% c("Redistributed", "Obstructed"), "Symmetric",
         ifelse(generation >= F, "Fall apart", ifelse(generation >= R, "Plateau", "Reorganize"))))
  factor(out, LP_LEVELS)
}

if (sys.nframe() == 0L) {
  e <- new.env(); load("data/derived/phase_boundaries.rda", envir = e)
  h <- new.env(); load("data/derived/he_gradient_boundary.rda", envir = h)
  b <- e$bounds[e$bounds$Scenario %in% c("Redistributed", "Obstructed"),
                c("Scenario", "Replicate", "R_gain", "F_hzn", "reorg")]
  names(b) <- c("Scenario", "Replicate", "R_gain", "F_gen", "reorganized_gain")
  b <- merge(b, h$hb[, c("Scenario", "Replicate", "R_he", "reorg_he")], by = c("Scenario", "Replicate"))
  b$R_gen <- b$R_he; b$reorganized <- b$reorg_he
  b$phase_policy_R_gen <- ifelse(b$reorganized, b$R_gen, b$F_gen)
  b <- b[, c("Scenario", "Replicate", "R_gen", "F_gen", "reorganized", "phase_policy_R_gen", "R_gain", "reorganized_gain")]
  write.csv(b, "data/derived/lineage_phases.csv", row.names = FALSE)
  cat("Wrote data/derived/lineage_phases.csv:", nrow(b), "lineages\n")
  lp_bounds <- b
  g <- expand.grid(Replicate = 1:50, Scenario = c("Redistributed", "Obstructed"), generation = seq(2004, 2999, 5),
                   stringsAsFactors = FALSE)
  g$phase <- lp_phase(g$Scenario, g$Replicate, g$generation, b)
  cat("\nCensuses per phase (every 5 generations):\n"); print(table(g$Scenario, droplevels(g$phase)))
  cat("\nLineages in each phase, by generation (Redistributed):\n")
  print(with(subset(g, Scenario == "Redistributed" & generation %in% seq(2004, 2999, 100)),
             table(generation, droplevels(phase))))
  cat("\nObstructed:\n")
  print(with(subset(g, Scenario == "Obstructed" & generation %in% seq(2004, 2999, 100)),
             table(generation, droplevels(phase))))
}

## ---- phase-share strip for trajectory figures ------------------------------------------
# One ordinal blue ramp for the three phases (validated with the dataviz palette
# checker: monotone lightness, single hue, light end >= 2:1 on white).
LP_COL <- c(Reorganize = "#86b6ef", Plateau = "#2a78d6", `Fall apart` = "#104281")
LP_LAB <- c(Reorganize = "Ascent", Plateau = "Plateau", `Fall apart` = "Decay")   # manuscript phase names

#' Share of the 50 lineages in each phase at each census
lp_share <- function(scenarios, gens = seq(2004, 2999, 5), bounds = lp_bounds, rcol = "phase_policy_R_gen") {
  g <- expand.grid(Replicate = 1:50, Scenario = scenarios, generation = gens, stringsAsFactors = FALSE)
  g$phase <- as.character(lp_phase(g$Scenario, g$Replicate, g$generation, bounds, rcol))
  g <- g[g$phase %in% names(LP_COL), ]
  if (!nrow(g)) return(data.frame(Scenario = character(), generation = numeric(), phase = factor(), share = numeric()))
  out <- as.data.frame(table(Scenario = g$Scenario, generation = g$generation, phase = factor(g$phase, names(LP_COL))),
                       stringsAsFactors = FALSE)
  out$generation <- as.numeric(out$generation); out$share <- out$Freq / 50
  out$phase <- factor(out$phase, names(LP_COL))
  out
}

#' Strip layer: stacked area of lineages per phase, faceted like the figure above it.
#' `levels` are the panel levels (may include scenarios without phases, e.g.
#' "Isotropic", which get a label instead of a strip). `facet` = "col" or "row".
lp_strip <- function(levels, xscale, facet = "col", ink = "#3d3d3a", base_size = 11,
                     none_label = "no imposed direction", strip_labels = FALSE) {
  sh <- lp_share(intersect(levels, c("Redistributed", "Obstructed")))
  sh$Scenario <- factor(sh$Scenario, levels)
  none <- data.frame(Scenario = factor(setdiff(levels, c("Redistributed", "Obstructed")), levels))
  p <- ggplot2::ggplot(sh, ggplot2::aes(generation, share, fill = phase)) +
    ggplot2::geom_area(stat = "identity", position = ggplot2::position_stack(reverse = TRUE), colour = NA) +
    ggplot2::scale_fill_manual(values = LP_COL, labels = LP_LAB, name = "Lineages in phase", drop = FALSE) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "25", "50"), limits = c(0, 1), expand = c(0, 0)) +
    xscale + ggplot2::labs(y = "Lineages") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = ink),
                   legend.position = "bottom", legend.title = ggplot2::element_text(colour = ink, size = base_size - 1),
                   legend.text = ggplot2::element_text(colour = ink, size = base_size - 1),
                   legend.key.size = ggplot2::unit(0.4, "cm"))
  if (nrow(none)) p <- p + ggplot2::geom_text(data = none, ggplot2::aes(x = 2500, y = 0.5, label = none_label),
                                              inherit.aes = FALSE, colour = "#8a8984", size = (base_size - 2) / ggplot2::.pt)
  if (facet == "col") p <- p + ggplot2::facet_grid(. ~ Scenario, drop = FALSE)
  else p <- p + ggplot2::facet_grid(Scenario ~ ., drop = FALSE)
  if (!strip_labels) p <- p + ggplot2::theme(strip.text = ggplot2::element_blank())
  p
}
