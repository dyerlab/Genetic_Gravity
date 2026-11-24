# R/simulation/permutation_nulls_1145.R
#
# Supplement "Permutation nulls for Delta", rerun under gstudio 1.14.5
# (local/gstudio_1145_verification_and_gamma_rule_instructions.md, Part 1d).
# The numbers in the supplement came from the private research repo
# (AsymmetricPopGraphs: R/gravity_existence_profile.R, R/gravity_location_profile.R;
# write-up local/archive/gravity_significance_tests.md). Both tests are the stock gstudio
# functions, called with return_null = TRUE so the null draws (needed for z-scores and
# null SDs) come back with the result.
#
# Rewiring null (asymmetry_network(rewire = "degree")): statistic mean |Delta|;
#   null = randomize_graph(mode = "degree") (1.14.5: igraph::rewire(keeping_degseq()))
#   with the observed edge weights shuffled on; one-sided upper tail; nperm = 99 as in
#   the original run. Samples: burn-in 1904-1999 every 5 generations (1,000 censuses)
#   and Symmetric / Redistributed / Obstructed 2004-2954 every 50 (1,000 each).
# Label-permutation null (asymmetry_permutation()): individuals' population labels
#   permuted, popgraph refit (alpha = 1), weights carried onto the fixed edges;
#   nperm = 499 and seeds 1e5 * rep + gen as in the original window run, so Delta and
#   the null draws match it. 1.14.5 p is centred on each edge's null mean; the
#   uncentred p (1.14.4) is computed from the same draws for comparison.
#   Sample: burn-in 1904-1999 every 5 generations, 50 replicates (1,000 censuses).
#
# Run from the repository root:  Rscript R/simulation/permutation_nulls_1145.R
#   CHECK=1 ... only recomputes a few cached censuses and compares (same seeds).
#   PART=rewiring / PART=label ... runs one part. Each census is cached under
#   data/derived/permutation_nulls_1145_cache/, so an interrupted run resumes.
# Output: data/derived/permutation_nulls_1145.rda (`rew`, `loc_edges`, `loc_cens`,
#   `checks`, `meta`), assembled from data/derived/permutation_nulls_1145_{rewiring,label}.rda.

suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(parallel) })
source("R/common/gravity_convention.R")      # GC_SIM_DIR, gc_load_graph()

OUT   <- "data/derived/permutation_nulls_1145.rda"
SEED  <- 20261001L
CORES <- max(1L, detectCores() - 1L)
SCN   <- c(Symmetric = 1L, Redistributed = 2L, Obstructed = 3L)
geno_file <- function(rep, gen, sc = NULL) if (is.null(sc))
  sprintf("%s/replicate%d/rep%d-genotypes-%d.rda", GC_SIM_DIR, rep, rep, gen) else
  sprintf("%s/replicate%d/rep%d-genotypes-scenario%d-%d.rda", GC_SIM_DIR, rep, rep, sc, gen)

## ---- the two tests, with their null draws ----------------------------------------------------
#' Rewiring null: asymmetry_network(); obs = mean |Delta|, null = mean |Delta| of each valid rewired graph
rewire_null <- function(graph, nperm = 99L, rewire = "degree") {
  a <- asymmetry_network(graph, nperm = nperm, rewire = rewire, return_null = TRUE)
  list(obs = unname(a$statistic), null = a$null_distribution, p = a$p.value)
}
#' Label-permutation null: asymmetry_permutation(); per-edge null mean, SD and p (centred, as
#' gstudio reports it, and uncentred, as gstudio 1.14.4 did) from the returned null matrix
label_null <- function(graph, data, groups, nperm = 499L) {
  a <- suppressWarnings(asymmetry_permutation(graph, data, groups, nperm = nperm, return_null = TRUE))
  null <- attr(a, "null_distribution"); B <- rowSums(!is.na(null)); centre <- rowMeans(null, na.rm = TRUE)
  data.frame(from = a$from, to = a$to, delta = a$delta, null_mean = centre,
             null_sd = apply(null, 1, sd, na.rm = TRUE), B = B, p_centred = a$p.value,
             p_uncentred = (1 + rowSums(abs(null) >= abs(a$delta), na.rm = TRUE)) / (1 + B),
             stringsAsFactors = FALSE)
}

load_geno <- function(f) { e <- new.env(); load(f, envir = e)
  list(mv = e$genotypes |> select(-any_of(c("Population", "ID"))) |> to_mv(), grp = as.factor(e$genotypes$Population)) }

## ---- jobs ------------------------------------------------------------------------------------
rew_jobs <- rbind(
  expand.grid(Replicate = 1:50, sample = "Burn-in", generation = seq(1904L, 1999L, 5L), stringsAsFactors = FALSE),
  expand.grid(Replicate = 1:50, sample = names(SCN), generation = seq(2004L, 2954L, 50L), stringsAsFactors = FALSE))
loc_jobs <- expand.grid(Replicate = 1:50, generation = seq(1904L, 1999L, 5L))
rew_graph <- function(j) gc_load_graph(j$Replicate, j$generation, if (j$sample == "Burn-in") NULL else SCN[[j$sample]])

## ---- reproduction check against cached censuses (same seeds) ---------------------------------
CACHE <- "data/derived/permutation_nulls_1145_cache"   # one file per census; resumable
check_cache <- function() {
  set.seed(1); ri <- sample(nrow(rew_jobs), 4); li <- sample(nrow(loc_jobs), 2)
  r <- lapply(ri, function(i) { f <- file.path(CACHE, "rewiring", sprintf("%05d.rds", i)); if (!file.exists(f)) return(NULL)
    old <- readRDS(f); g <- rew_graph(rew_jobs[i, ]); set.seed(SEED + i); b <- rewire_null(g, 99L)
    data.frame(test = "rewiring", job = i, max_diff = max(abs(c(b$obs - old$obs, b$p - old$p, mean(b$null) - old$null_mean)))) })
  l <- lapply(li, function(i) { f <- file.path(CACHE, "label", sprintf("%05d.rds", i)); if (!file.exists(f)) return(NULL)
    old <- readRDS(f); j <- loc_jobs[i, ]; g <- gc_load_graph(j$Replicate, j$generation)
    G <- load_geno(geno_file(j$Replicate, j$generation)); set.seed(1e5 * j$Replicate + j$generation); b <- label_null(g, G$mv, G$grp, 499L)
    data.frame(test = "label", job = i, max_diff = max(abs(as.matrix(b[, 3:8]) - as.matrix(old[, names(b)[3:8]])))) })
  do.call(rbind, c(r, l))
}
checks <- check_cache()
print(checks, digits = 3)
if (nzchar(Sys.getenv("CHECK"))) quit(save = "no")
if (!is.null(checks)) stopifnot(all(checks$max_diff < 1e-12))

## ---- run (each part cached separately; PART=rewiring or PART=label runs one) --------------
OUT_REW <- "data/derived/permutation_nulls_1145_rewiring.rda"
OUT_LOC <- "data/derived/permutation_nulls_1145_label.rda"
part <- Sys.getenv("PART")
dir.create(file.path(CACHE, "rewiring"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CACHE, "label"), recursive = TRUE, showWarnings = FALSE)
cached <- function(sub, i, fn) { f <- file.path(CACHE, sub, sprintf("%05d.rds", i))
  if (file.exists(f)) return(readRDS(f)); v <- fn(); saveRDS(v, f); v }
if (!file.exists(OUT_REW) && part %in% c("", "rewiring")) {
  t0 <- Sys.time()
  rew <- bind_rows(mclapply(seq_len(nrow(rew_jobs)), function(i) cached("rewiring", i, function() {
    j <- rew_jobs[i, ]; g <- rew_graph(j); if (is.null(g)) return(NULL)
    set.seed(SEED + i)
    r <- tryCatch(rewire_null(g, 99L), error = function(e) NULL); if (is.null(r)) return(cbind(j, failed = TRUE))
    ok <- r$null[!is.na(r$null)]
    cbind(j, edges = ecount(g), obs = r$obs, n_null = length(ok), null_mean = mean(ok), null_sd = sd(ok),
          z = (r$obs - mean(ok)) / sd(ok), p = r$p, failed = FALSE)
  }), mc.cores = CORES, mc.preschedule = FALSE))
  wall_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  save(rew, checks, wall_min, file = OUT_REW)
}
if (!file.exists(OUT_LOC) && part %in% c("", "label")) {
  t0 <- Sys.time()
  loc_edges <- bind_rows(mclapply(seq_len(nrow(loc_jobs)), function(i) cached("label", i, function() {
    j <- loc_jobs[i, ]; g <- gc_load_graph(j$Replicate, j$generation); if (is.null(g)) return(NULL)
    G <- load_geno(geno_file(j$Replicate, j$generation))
    set.seed(1e5 * j$Replicate + j$generation)
    d <- tryCatch(label_null(g, G$mv, G$grp, 499L), error = function(e) NULL); if (is.null(d)) return(NULL)
    k <- degree(g)
    cbind(Replicate = j$Replicate, generation = j$generation, d, k_from = k[d$from], k_to = k[d$to])
  }), mc.cores = CORES, mc.preschedule = FALSE))
  wall_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  save(loc_edges, checks, wall_min, file = OUT_LOC)
}
if (file.exists(OUT_REW) && file.exists(OUT_LOC) && !file.exists(OUT)) {
  e1 <- new.env(); load(OUT_REW, envir = e1); e2 <- new.env(); load(OUT_LOC, envir = e2)
  rew <- e1$rew; loc_edges <- e2$loc_edges
  loc_cens <- loc_edges |> group_by(Replicate, generation) |> summarise(edges = n(), .groups = "drop")
  meta <- list(date = as.character(Sys.Date()), gstudio = as.character(packageVersion("gstudio")),
               igraph = as.character(packageVersion("igraph")), R = R.version.string,
               seed_rewiring = SEED, nperm_rewiring = 99L, nperm_label = 499L,
               wall_min = c(rewiring = e1$wall_min, label = e2$wall_min))
  save(rew, loc_edges, loc_cens, checks, meta, file = OUT)
  cat(sprintf("rewiring: %d censuses; label: %d censuses, %d edges\n", nrow(rew), nrow(loc_cens), nrow(loc_edges)))
}
