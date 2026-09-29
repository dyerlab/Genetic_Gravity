# exploratory/ss_bandwidth_by_phase.R
#
# Bandwidth sensitivity of the source-sink gradient test, by each lineage's own
# phase (replaces the calendar-peak-window version in
# exploratory/node_potential_reanalysis.R, section `bw`, and proposal 9).
#
# Censuses: the round-1 grid, generations 2004-2954 every 50, 50 replicates, for
# Isotropic, Redistributed and Obstructed (3,000 censuses). Phase from
# exploratory/lineage_phases.R (Isotropic = "Symmetric", the FPR reference).
# For each census, S is recomputed under five bandwidths (canonical local,
# fixed global = mean edge weight, local x1/2, local x2, perplexity search with
# perplexity 4), each through exploratory/gravity_convention.R, and tested with
# the locked spec: Spearman r(S, x), Moran spectral randomization
# (adespatial::msr, "singleton") on the census's own graph with weights
# exp(-e^2 / 2 bbar^2), B = 999, two-sided, alpha = 0.05. The null weights do not
# depend on the bandwidth.
#
# Run from the repository root:  Rscript exploratory/ss_bandwidth_by_phase.R
# Output: data/derived/ss_bandwidth_by_phase.rda (`bwp`: per census x bandwidth; `meta`)

suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(tidyr); library(parallel) })
source("exploratory/gravity_convention.R")
source("exploratory/lineage_phases.R")
OUT <- "data/derived/ss_bandwidth_by_phase.rda"
SEED <- 20260930L; B <- 999L
SC <- c(Isotropic = 1L, Redistributed = 2L, Obstructed = 3L)
ALTS <- list(canonical = list(), global = "global", `local x1/2` = list(scale = 0.5),
             `local x2` = list(scale = 2), perplexity = list(perplexity = 4))

if (file.exists(OUT)) load(OUT) else {
  jobs <- expand.grid(Replicate = 1:50, scenario = names(SC), generation = seq(2004L, 2954L, 50L), stringsAsFactors = FALSE)
  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(jobs)), function(i) {
    set.seed(SEED + i)
    g <- gc_load_graph(jobs$Replicate[i], jobs$generation[i], SC[[jobs$scenario[i]]]); if (is.null(g)) return(NULL)
    nodes <- sort(V(g)$name); x <- gc_x(nodes); bbar <- mean(E(g)$weight)
    Aw <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)[nodes, nodes]
    Aw[Aw > 0] <- exp(-Aw[Aw > 0]^2 / (2 * bbar^2))
    lw <- spdep::mat2listw(Aw, style = "B", zero.policy = TRUE)
    out <- list()
    for (nm in names(ALTS)) {
      a <- ALTS[[nm]]
      cs <- tryCatch(if (identical(a, "global")) gc_census(g, bandwidth = bbar) else do.call(gc_census, c(list(g), a)),
                     error = function(e) NULL)
      if (is.null(cs)) { out[[nm]] <- data.frame(jobs[i, ], bandwidth = nm, r_S = NA, p = NA, failed = TRUE); next }
      S <- cs$nodes$S[match(nodes, cs$nodes$node)]
      r_obs <- suppressWarnings(cor(S, x, method = "spearman"))
      dr <- as.matrix(adespatial::msr(S, lw, nrepet = B, method = "singleton", simplify = TRUE))
      r_null <- as.numeric(cor(dr, x, method = "spearman"))
      out[[nm]] <- data.frame(jobs[i, ], bandwidth = nm, r_S = r_obs,
                              p = (1 + sum(abs(r_null) >= abs(r_obs) - 1e-12)) / (B + 1), failed = FALSE)
    }
    do.call(rbind, out)
  }, mc.cores = max(1L, detectCores() - 1L), mc.preschedule = FALSE)
  bwp <- do.call(rbind, res)
  meta <- list(seed = SEED, B = B, wall_min = as.numeric(difftime(Sys.time(), t0, units = "mins")),
               gstudio = as.character(packageVersion("gstudio")), adespatial = as.character(packageVersion("adespatial")),
               censuses = nrow(jobs), date = as.character(Sys.Date()))
  save(bwp, meta, file = OUT)
}
bwp$phase <- lp_phase(bwp$scenario, bwp$Replicate, bwp$generation)
bwp$bandwidth <- factor(bwp$bandwidth, names(ALTS))

if (sys.nframe() == 0L) {
  options(width = 170)
  cat(sprintf("Censuses: %d; wall time %.1f min; failures by bandwidth:\n", meta$censuses, meta$wall_min))
  print(with(bwp, tapply(failed, list(bandwidth, phase), sum))[, c("Symmetric", "Reorganize", "Plateau", "Fall apart")])
  ok <- filter(bwp, !failed)
  tab <- ok |> group_by(scenario, phase, bandwidth) |>
    summarise(censuses = n(), lineages = n_distinct(Replicate), reject = mean(p < .05), med_r = median(r_S),
              frac_up = mean(r_S < 0), .groups = "drop") |>
    left_join(ok |> group_by(scenario, phase, bandwidth, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |>
                group_by(scenario, phase, bandwidth) |> summarise(lin_up = mean(m < 0), .groups = "drop"),
              by = c("scenario", "phase", "bandwidth"))
  print(as.data.frame(tab |> filter(!(scenario == "Obstructed" & phase == "Fall apart"))), digits = 3)
  cat("\nCanonical row vs round-1 locked test (same censuses, different null seeds): rejection by phase\n")
  e <- new.env(); load("data/derived/node_potential_reanalysis.rda", envir = e)
  r1 <- e$rean$sstest$tests |> filter(null == "N1w", scenario %in% names(SC)) |>
    mutate(phase = lp_phase(scenario, Replicate, generation)) |> group_by(scenario, phase) |> summarise(r1 = mean(p < .05), .groups = "drop")
  print(left_join(filter(tab, bandwidth == "canonical") |> select(scenario, phase, reject), r1, by = c("scenario", "phase")) |> as.data.frame(), digits = 3)
  cat("\nPaired agreement with canonical on rejection (plateau, both asymmetric scenarios):\n")
  w <- ok |> filter(phase == "Plateau") |> select(scenario, Replicate, generation, bandwidth, p) |>
    pivot_wider(names_from = bandwidth, values_from = p)
  for (nm in names(ALTS)[-1]) { z <- w[!is.na(w[[nm]]) & !is.na(w$canonical), ]
    cat(sprintf("  %-11s both reject %.3f, canonical only %.3f, alternative only %.3f (n = %d)\n", nm,
                mean(z$canonical < .05 & z[[nm]] < .05), mean(z$canonical < .05 & z[[nm]] >= .05),
                mean(z$canonical >= .05 & z[[nm]] < .05), nrow(z))) }
}
