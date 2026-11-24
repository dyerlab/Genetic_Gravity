# R/simulation/density_power_full.R
#
# FULL PANEL: power and false-positive rates of the two asymmetry tests against sampling
# design (number of populations K, individuals per population n, loci L), by subsampling
# saved simulation genotypes and rebuilding the Population Graph. Follows the prototype
# (density_power_prototype.R, not in the repository).
#
# Censuses (seed 20261123):
#   Symmetric (scenario 1)  50 lineages x N_SYM censuses, forward phase 2004-2954 (50-gen grid)
#   sym-low   (scenario 5)  50 lineages x N_LOW censuses, same grid (drift check of the cutoff)
#   Redistributed, Obstructed: every lineage with a plateau, N_ASY censuses from the plateau
#                           core (middle 50% of [R', F); data/derived/lineage_phases.csv)
# Schemes per census (x = original deme index throughout):
#   spread      K in K_MAIN x n in N_MAIN, all loci. Stratified draw: demes 1..25 cut into
#               K contiguous blocks, one deme per block (spans the chain).
#   contiguous  K in K_CONT x n in N_CONT, all loci: a block of K adjacent demes at a random
#               start (same count, narrower span: separates population count from spatial
#               resolution). Not run for sym-low.
#   loci        L in L_SUB x K in K_LOCI, n = 30, spread draw; L loci drawn from those
#               polymorphic in the census. Not run for sym-low.
#   N_DRAW draws per cell (one for the unsubsampled K = 25, n = 100, all loci).
# Per graph: edges, mean degree, components, largest-component size, leaves, terminal demes
#   kept; degree share and cor(S, S0) (gamma = 1/2); source_sink_test() (MSR on the graph,
#   B = 999) and ibgd(mode = "dcGD") (Delta R^2, b_fwd/b_rev) on the whole graph, and, when
#   the graph is disconnected, the same two on its largest component (>= 4 nodes).
# Per census: multilocus Fst (genetic_structure(mode = "Fst"), all 25 demes x 100), loci
#   polymorphic, for the Falconer (1981) sampling optimum (gstudio::optimal_sampling()).
#
# Run from the repository root (resumable; each census cached):
#   Rscript R/simulation/density_power_full.R          (env: CORES, N_SYM, N_LOW, N_ASY, N_DRAW)
# Output: data/derived/density_power_full.rda (`res` = one row per graph, `cens`, `meta`)

suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(parallel) })
source("R/common/gravity_convention.R")

GAMMA  <- 0.5
CORES  <- as.integer(Sys.getenv("CORES", max(1L, detectCores() - 1L)))
N_SYM  <- as.integer(Sys.getenv("N_SYM", "6")); N_LOW <- as.integer(Sys.getenv("N_LOW", "4"))
N_ASY  <- as.integer(Sys.getenv("N_ASY", "4")); N_DRAW <- as.integer(Sys.getenv("N_DRAW", "3"))
K_MAIN <- c(25L, 21L, 18L, 15L, 12L, 10L, 8L, 5L); N_MAIN <- c(100L, 50L, 30L, 20L, 10L)
K_CONT <- c(15L, 10L, 8L);  N_CONT <- c(100L, 30L)
K_LOCI <- c(25L, 18L, 12L, 8L); L_SUB <- c(10L, 8L)
SCN    <- c(Symmetric = 1L, Redistributed = 2L, Obstructed = 3L, `sym-low` = 5L)
CACHE  <- Sys.getenv("CACHE", "data/derived/density_power_full_cache"); dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
OUT    <- Sys.getenv("OUT", "data/derived/density_power_full.rda")

## ---- censuses -------------------------------------------------------------------------------
set.seed(20261123L)
lp <- read.csv("data/derived/lineage_phases.csv", stringsAsFactors = FALSE)
core <- lp |> mutate(Rp = ifelse(reorganized, R_gen, F_gen), Fi = ifelse(is.na(F_gen), 2999, F_gen)) |>
  filter(!is.na(Rp), Rp < Fi) |> mutate(a = Rp + (Fi - Rp) / 4, b = Fi - (Fi - Rp) / 4)
grid5 <- seq(2004L, 2999L, 5L); grid50 <- seq(2004L, 2954L, 50L)
pick_core <- function(s) bind_rows(lapply(core$Replicate[core$Scenario == s], function(r) {
  w <- core[core$Scenario == s & core$Replicate == r, ]; g <- grid5[grid5 >= w$a & grid5 <= w$b]
  if (!length(g)) return(NULL)
  data.frame(scenario = s, Replicate = r, generation = if (length(g) == 1) g else sort(sample(g, min(N_ASY, length(g))))) }))
pick_sym <- function(s, m) bind_rows(lapply(1:50, function(r)
  data.frame(scenario = s, Replicate = r, generation = sort(sample(grid50, m)))))
cens <- bind_rows(pick_sym("Symmetric", N_SYM), pick_sym("sym-low", N_LOW),
                  pick_core("Redistributed"), pick_core("Obstructed"))

## ---- graph construction and tests -------------------------------------------------------------
load_geno <- function(rep, gen, scn) {
  e <- new.env(); load(sprintf("%s/replicate%d/rep%d-genotypes-scenario%d-%d.rda", GC_SIM_DIR, rep, rep, scn, gen), envir = e)
  e$genotypes
}
draw_spread <- function(K) if (K >= 25L) 1:25 else
  vapply(split(1:25, cut(1:25, K, labels = FALSE)), function(b) b[sample.int(length(b), 1L)], 1L)
draw_contig <- function(K) { s <- sample.int(25L - K + 1L, 1L); s:(s + K - 1L) }

tests_on <- function(g) {
  x <- gc_xg(g)
  v <- c(deg_share = NA, cor_S_S0 = NA, r_S = NA, p_S = NA, r2_cgd = NA, r2_dcgd = NA, dR2 = NA, ratio = NA)
  if (vcount(g) < 4 || ecount(g) < 2) return(v)
  t_s <- tryCatch(source_sink_test(g, x, GAMMA), error = function(e) NULL)
  if (!is.null(t_s)) v[1:4] <- c(attr(t_s$scores, "degree_share"), attr(t_s$scores, "cor_S_S0"), unname(t_s$statistic), t_s$p.value)
  ib <- tryCatch(ibgd(g, x, mode = "dcGD", gamma = GAMMA), error = function(e) NULL)
  if (!is.null(ib)) v[5:8] <- c(ib$estimate[c("R2 cGD", "adj R2 dcGD")], unname(ib$statistic), ib$estimate[["slope ratio"]])
  v
}

one_graph <- function(gt, loci, demes, n) {
  pops <- sprintf("Pop%02d", demes)
  keep <- unlist(lapply(pops, function(p) { i <- which(gt$Population == p); if (n < length(i)) sample(i, n) else i }))
  sub <- gt[keep, ]
  g <- popgraph(to_mv(sub[, loci, drop = FALSE]), groups = as.factor(as.character(sub$Population)))
  k <- degree(g); cm <- components(g)
  v <- c(n_edges = ecount(g), mean_deg = mean(k), n_comp = cm$no, lc_size = max(cm$csize), n_leaf = sum(k == 1),
         ends = all(c(1, 25) %in% demes), span = diff(range(demes)), tests_on(g))
  lc <- setNames(rep(NA_real_, 8), paste0(names(tests_on(make_empty_graph())), "_lc"))
  if (cm$no > 1 && max(cm$csize) >= 4) {
    h <- induced_subgraph(g, which(cm$membership == which.max(cm$csize)))
    lc[] <- tests_on(h)
  }
  c(v, lc)
}

run_census <- function(i) {
  cz <- cens[i, ]; f <- file.path(CACHE, sprintf("%s_%02d_%d.rds", gsub("[^A-Za-z]", "", cz$scenario), cz$Replicate, cz$generation))
  if (file.exists(f)) return(readRDS(f))
  set.seed(20261123L + i)
  gt <- load_geno(cz$Replicate, cz$generation, SCN[[cz$scenario]])
  all_loci <- setdiff(names(gt), c("ID", "Population"))
  poly <- all_loci[vapply(all_loci, function(l) length(unique(as.character(gt[[l]]))) > 1, TRUE)]
  fst <- tryCatch({ s <- genetic_structure(gt, mode = "Fst"); s$Fst[s$Locus == "Multilocus"] }, error = function(e) NA_real_)
  cells <- expand.grid(scheme = "spread", K = K_MAIN, n = N_MAIN, L = 20L, stringsAsFactors = FALSE)
  if (cz$scenario != "sym-low") cells <- rbind(cells,
    expand.grid(scheme = "contiguous", K = K_CONT, n = N_CONT, L = 20L, stringsAsFactors = FALSE),
    expand.grid(scheme = "loci", K = K_LOCI, n = 30L, L = L_SUB, stringsAsFactors = FALSE))
  rows <- list()
  for (c_ in seq_len(nrow(cells))) {
    ce <- cells[c_, ]; nd <- if (ce$scheme == "spread" && ce$K == 25L && ce$n == 100L) 1L else N_DRAW
    for (d in seq_len(nd)) {
      demes <- if (ce$scheme == "contiguous") draw_contig(ce$K) else draw_spread(ce$K)
      loci  <- if (ce$scheme == "loci") sample(poly, min(ce$L, length(poly))) else all_loci
      v <- tryCatch(one_graph(gt, loci, demes, ce$n), error = function(e) NULL)
      if (!is.null(v)) rows[[length(rows) + 1]] <- data.frame(cz, ce, draw = d, L_eff = length(intersect(loci, poly)),
                                                             demes = paste(demes, collapse = "-"), t(v))
    }
  }
  out <- bind_rows(rows) |> mutate(fst = fst, n_poly = length(poly))
  saveRDS(out, f); out
}

t0 <- Sys.time()
res <- bind_rows(mclapply(seq_len(nrow(cens)), function(i) tryCatch(run_census(i), error = function(e)
  data.frame(cens[i, ], err = conditionMessage(e))), mc.cores = CORES, mc.preschedule = FALSE))
meta <- list(gamma = GAMMA, K_MAIN = K_MAIN, N_MAIN = N_MAIN, K_CONT = K_CONT, N_CONT = N_CONT, K_LOCI = K_LOCI, L_SUB = L_SUB,
             N_SYM = N_SYM, N_LOW = N_LOW, N_ASY = N_ASY, N_DRAW = N_DRAW, date = as.character(Sys.Date()),
             gstudio = as.character(packageVersion("gstudio")), cores = CORES,
             wall_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
save(res, cens, meta, file = OUT)
cat(sprintf("Wrote %s: %d graphs from %d censuses in %.1f min\n", OUT, nrow(res), nrow(cens), meta$wall_min))
