# R/simulation/fresh_replicates.R
#
# Part 3 of local/gstudio_1145_verification_and_gamma_rule_instructions.md: 50 new
# lineages, simulated with the manuscript's driver and parameters but new seeds, to
# test gamma = 1, gamma = 0.5 and gamma* out of sample.
#
# Simulation: the functions and constants of the private research repo's
# R/simulate.R (read, not modified): K = 25 demes, N = 100, 20 biallelic loci, symmetric
# burn-in at m = 0.025 for 2,000 generations, then forward conditions 1 (Symmetric),
# 2 (Redistributed), 3 (Obstructed) and 5 (sym-low) for 1,000 generations each, branching
# from the generation-1999 state. Lineages are numbered 1001-1050, so the driver's seed
# rules (burn-in rep * 10000; forward rep * 10000 + cond) give seeds disjoint from the
# manuscript's lineages 1-70.
# Saved per lineage x condition: every 5 generations the interior He slope and the mean
# number of loci polymorphic within a deme (the two phase markers; burn-in 1904-1999 and
# forward 2004-2999), and Population Graphs (gstudio::popgraph, as build_graph())
# on the 50-generation grid 2004-2954 only.
# Phases: the current rule (R/simulation/he_gradient_boundary.R, lineage_phases.R): R =
# first generation the penalized-spline smooth (k = 10) of the He slope completes 90%
# of its rise from the lineage's burn-in mean; reorganized if the rise exceeds the 95th
# percentile of the fresh Symmetric lineages' rises; F = first census with < 4.5
# polymorphic loci per deme; unflagged lineages held in ascent until F.
# Tests on the 50-generation grid at gamma = 1, 0.5 and gamma*, with stock gstudio:
# source_sink_test() (B = 999, seed 20261002 + census index) and ibgd(mode = "dcGD") (Delta R^2).
#
# Run from the repository root:  Rscript R/simulation/fresh_replicates.R
#   (resumable: each lineage x condition is cached under data/derived/fresh_replicates/)
# Output: data/derived/fresh_replicates.rda (`fr_stats`, `fr_bounds`, `fr_tests`, `meta`)

suppressPackageStartupMessages({ library(gstudio); library(igraph); library(dplyr); library(parallel); library(mgcv) })
source("R/common/gravity_convention.R")      # gc_x(), gc_xg()
SIMR  <- file.path(Sys.getenv("GG_SIM_REPO", "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs"), "R", "simulate.R")
FDIR  <- "data/derived/fresh_replicates"; OUTF <- "data/derived/fresh_replicates.rda"
REPS  <- 1001:1050; CONDS <- c(Symmetric = 1L, Redistributed = 2L, Obstructed = 3L, `sym-low` = 5L)
CORES <- as.integer(Sys.getenv("MC_CORES", unset = max(1L, detectCores() - 1L)))
dir.create(FDIR, showWarnings = FALSE)

## ---- simulator functions and constants, as in R/simulate.R ----------------------------------
local({
  src <- readLines(SIMR)
  keep <- c(1:(grep("Run controls", src)[1] - 1), grep("^stepping_stone_1d_asym <- ", src):(grep("^# Save genotype snapshot", src) - 1))
  suppressMessages(eval(parse(text = src[keep]), envir = globalenv()))
})
burnin_seed  <- function(rep) rep * 10000L
forward_seed <- function(rep, cond) rep * 10000L + cond

markers <- function(pop) {                                    # he_census() and n_poly_pop on a population in memory
  L <- grep("^L[0-9]", names(pop))
  M <- vapply(L, function(j) { a <- as.character(pop[[j]]); (substr(a, 1, 2) == "01") + (substr(a, 4, 5) == "01") }, numeric(nrow(pop)))
  pops <- sort(unique(pop$Population)); n <- as.vector(table(pop$Population)[pops])
  p <- rowsum(M, pop$Population)[pops, , drop = FALSE] / (2 * n)
  he <- rowMeans(2 * p * (1 - p)); x <- as.integer(sub("\\D+", "", pops)); int <- x >= 3 & x <= 23
  c(slope_int = unname(coef(lm(he[int] ~ x[int]))[2]), n_poly_pop = mean(rowSums(p > 0 & p < 1)))
}

run_burnin <- function(rep) {
  f <- file.path(FDIR, sprintf("burnin_%d.rds", rep)); if (file.exists(f)) return(invisible())
  set.seed(burnin_seed(rep))
  pop <- make_populations(make_initial_freqs(K, L, POP_NAMES, LOCUS_NAMES), N = 100L)
  m <- stepping_stone_1d_asym(K, M_BURN, M_BURN, POP_NAMES); st <- list()
  for (gen in seq_len(N_GEN_BURNIN - 1L)) {
    pop <- one_generation(pop, m, 100L)
    if ((gen + 1L) %% 5L == 0L && gen >= 1904L) st[[length(st) + 1]] <- data.frame(generation = gen, t(markers(pop)))
  }
  saveRDS(list(stats = do.call(rbind, st), pop1999 = pop), f)
}

run_forward <- function(rep, cond) {
  f <- file.path(FDIR, sprintf("fwd_%d_%d.rds", rep, cond)); if (file.exists(f)) return(invisible())
  b <- readRDS(file.path(FDIR, sprintf("burnin_%d.rds", rep))); pop <- b$pop1999
  row <- FORWARD[FORWARD$cond == cond, ]
  m <- stepping_stone_1d_asym(K, row$m_fwd, row$m_rev, POP_NAMES)
  set.seed(forward_seed(rep, cond)); st <- list(); graphs <- list()
  for (step in seq_len(N_GEN_FWD)) {
    pop <- one_generation(pop, m, 100L); gen <- 1999L + step
    if (step %% 5L == 0L) st[[length(st) + 1]] <- data.frame(generation = gen, t(markers(pop)))
    if ((gen - 2004L) %% 50L == 0L && gen <= 2954L)
      graphs[[as.character(gen)]] <- tryCatch(build_graph(pop), error = function(e) NULL)
  }
  saveRDS(list(stats = do.call(rbind, st), graphs = graphs), f)
}

if (!file.exists(OUTF)) {
  t0 <- Sys.time()
  invisible(mclapply(REPS, run_burnin, mc.cores = CORES, mc.preschedule = FALSE))
  fj <- expand.grid(rep = REPS, cond = CONDS)
  invisible(mclapply(seq_len(nrow(fj)), function(i) run_forward(fj$rep[i], fj$cond[i]), mc.cores = CORES, mc.preschedule = FALSE))
  t_sim <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  ## ---- phase markers and boundaries ----------------------------------------------------------
  fr_stats <- bind_rows(
    bind_rows(lapply(REPS, function(r) cbind(Replicate = r, sample = "Burn-in", readRDS(file.path(FDIR, sprintf("burnin_%d.rds", r)))$stats))),
    bind_rows(lapply(seq_len(nrow(fj)), function(i) cbind(Replicate = fj$rep[i], sample = names(CONDS)[CONDS == fj$cond[i]],
                                                          readRDS(file.path(FDIR, sprintf("fwd_%d_%d.rds", fj$rep[i], fj$cond[i])))$stats))))
  grid <- seq(2004, 2999, 5)
  sm <- function(g, y) as.numeric(predict(gam(y ~ s(g, k = 10), data = data.frame(g = g, y = y), method = "REML"), data.frame(g = grid)))
  base <- fr_stats |> filter(sample == "Burn-in") |> group_by(Replicate) |> summarise(b = mean(slope_int), .groups = "drop")
  fr_bounds <- bind_rows(lapply(c("Symmetric", "Redistributed", "Obstructed"), function(s) bind_rows(lapply(REPS, function(r) {
    x <- fr_stats |> filter(sample == s, Replicate == r) |> arrange(generation)
    f <- sm(x$generation, x$slope_int); b0 <- base$b[base$Replicate == r]; amp <- max(f) - b0; i90 <- which(f >= b0 + 0.9 * amp)[1]
    data.frame(Scenario = s, Replicate = r, he_amp = amp, R_gen = if (is.na(i90)) NA else grid[i90],
               F_gen = x$generation[x$n_poly_pop < 4.5][1])
  }))))
  thr <- quantile(fr_bounds$he_amp[fr_bounds$Scenario == "Symmetric"], 0.95)
  thr_ms <- local({ h <- new.env(); load("data/derived/he_gradient_boundary.rda", envir = h); h$thr })
  fr_bounds <- fr_bounds |> mutate(reorganized = he_amp > thr, reorganized_ms_thr = he_amp > thr_ms,
                                   phase_policy_R_gen = ifelse(reorganized, R_gen, F_gen))

  ## ---- tests on the 50-generation grid -------------------------------------------------------
  # gamma* (the zero-degree-correlation bandwidth): root in log2(gamma) in [-2, 1] of cor(S(gamma), S0);
  # the 31-point grid is scanned first and the largest sign change refined by uniroot (as
  # gamma_star_rule.R, not in the repository)
  GRID <- seq(-2, 1, 0.1)
  find_root <- function(fn) {
    v <- vapply(GRID, fn, 0)
    if (all(is.na(v))) return(list(lg = NA_real_, flag = "na"))
    s <- sign(v); cr <- which(s[-1] * s[-length(s)] <= 0 & !is.na(s[-1]) & !is.na(s[-length(s)]))
    if (!length(cr)) return(list(lg = GRID[which.min(abs(v))], flag = "no_root"))
    i <- max(cr)
    lg <- if (v[i + 1] == 0) GRID[i + 1] else if (v[i] == 0) GRID[i] else
      uniroot(fn, c(GRID[i], GRID[i + 1]), f.lower = v[i], f.upper = v[i + 1], tol = 1e-4)$root
    list(lg = lg, flag = if (length(cr) > 1) "multi_root" else "")
  }
  tj <- expand.grid(Replicate = REPS, sample = names(CONDS), generation = seq(2004L, 2954L, 50L), stringsAsFactors = FALSE)
  tj$seed <- 20261002L + seq_len(nrow(tj))
  fr_tests <- bind_rows(mclapply(seq_len(nrow(tj)), function(i) tryCatch({
    j <- tj[i, ]; g <- readRDS(file.path(FDIR, sprintf("fwd_%d_%d.rds", j$Replicate, CONDS[[j$sample]])))$graphs[[as.character(j$generation)]]
    if (is.null(g)) return(NULL)
    x <- gc_xg(g); sc0 <- source_sink_scores(g); ok <- sc0$degree > 0; S0c <- sc0$S0
    S_at <- function(gam) source_sink_scores(g, gam)$S
    rp <- if (max(abs(S0c[ok])) < 1e-12) list(lg = log2(0.5), flag = "flat_degree") else
      find_root(function(lg) suppressWarnings(cor(S_at(2^lg)[ok], S0c[ok])))
    gset <- c(`1` = 1, `0.5` = 0.5, star = 2^rp$lg)
    bind_rows(lapply(names(gset), function(nm) { gam <- gset[[nm]]
      set.seed(j$seed); t1 <- source_sink_test(g, x, gam)
      ib <- tryCatch(ibgd(g, x, mode = "dcGD", gamma = gam)$estimate, error = function(e) c(`R2 cGD` = NA, `adj R2 dcGD` = NA, `slope ratio` = NA))
      data.frame(j[, 1:3], setting = nm, gamma = gam, flag = rp$flag, edges = ecount(g), r_S = unname(t1$statistic), p = t1$p.value,
                 dR2 = unname(ib["adj R2 dcGD"] - ib["R2 cGD"]), rel_fwd = unname(ib["slope ratio"])) }))
  }, error = function(e) data.frame(tj[i, 1:3], err = conditionMessage(e))), mc.cores = CORES, mc.preschedule = FALSE))
  meta <- list(date = as.character(Sys.Date()), reps = REPS, conds = CONDS, thr_fresh = thr, thr_manuscript = thr_ms,
               sim_min = t_sim, total_min = as.numeric(difftime(Sys.time(), t0, units = "mins")),
               gstudio = as.character(packageVersion("gstudio")))
  save(fr_stats, fr_bounds, fr_tests, meta, file = OUTF)
  cat(sprintf("Wrote %s (%.0f min)\n", OUTF, meta$total_min))
}
