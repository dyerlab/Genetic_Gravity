# exploratory/node_potential_reanalysis.R
#
# Reanalysis in the manuscript's sign convention (local/node_potential_reanalysis.md).
# Every signed quantity comes from exploratory/gravity_convention.R, or from the
# round-1 potential phi via S = -phi. Existing simulated data only.
#
# Sections (each cached in the output list; delete OUT to recompute):
#   truth     Task 3a  w_{j|i} and divMigrate Nm against the true migration
#                      matrix, as run (m_{i->j}) and transposed (m_{j->i}).
#   simil     Task 3c  sigma_ij = 1/(1 + P_ij) against Nm, as run and with Nm
#                      transposed.
#   ibgd_sym           Delta R^2 on the symmetric gradient (sym-mid/low/verylow),
#                      every 5 generations 2004-2999, for Tasks 4-5 and the
#                      supplement's Differentiation without Direction.
#   overlap   Task 4   source-sink gradient test vs Delta R^2 > 0, per census.
#   sstest    Task 5   calibration, power, sign consistency (round-1 output, new sign).
#   ident     Task 7b  identities in the new notation.
#   bw        Task 8   bandwidth sensitivity of the locked test.
#
# Run from the repository root:  Rscript exploratory/node_potential_reanalysis.R
# Output: data/derived/node_potential_reanalysis.rda (list `rean`).

suppressPackageStartupMessages({ library(igraph); library(gstudio); library(dplyr); library(tidyr)
  library(parallel); library(ggplot2) })
source("exploratory/gravity_convention.R")
source("exploratory/ibd_cgd_vs_pgd.R")          # ibd_one(); main block guarded

SEED  <- 20260930L
B     <- 999L
OUT   <- "data/derived/node_potential_reanalysis.rda"
DM    <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data/divMigrate"
POPS  <- sprintf("Pop%02d", 1:25)
CORES <- max(1L, detectCores() - 1L)
LAB   <- c(iso = "Isotropic", flux = "Redistributed", rate = "Obstructed")
RATES <- list(iso = c(0.025, 0.025), flux = c(0.040, 0.010), rate = c(0.025, 0.010))
BLOCK <- function(g) cut(g, c(1900, 2000, 2250, 2500, 2750, 3000),
                         labels = c("1904-1999", "2004-2249", "2254-2499", "2504-2749", "2754-2999"))

set.seed(SEED)
rean <- if (file.exists(OUT)) { e <- new.env(); load(OUT, envir = e); e$rean } else list()
save_rean <- function() save(rean, file = OUT)
timed <- function(name, expr) {
  if (!is.null(rean[[name]])) { cat(sprintf("[%s] cached\n", name)); return(invisible()) }
  t0 <- Sys.time(); val <- expr
  val$wall_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  rean[[name]] <<- val; save_rean()
  cat(sprintf("[%s] %.1f min\n", name, val$wall_min))
}

true_migration <- function(forward, reverse) {           # [i, j] = m_{i->j}
  K <- length(POPS); M <- matrix(0, K, K, dimnames = list(POPS, POPS))
  i <- seq_len(K - 1L); M[cbind(i, i + 1L)] <- forward; M[cbind(i + 1L, i)] <- reverse
  diag(M) <- NA; M
}
rho  <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))
rmse <- function(a, b) sqrt(mean((a - b)^2))

## ---- Task 3a: truth, as run and transposed ----------------------------------------
timed("truth", {
  off <- row(diag(25)) != col(diag(25))
  rows <- mclapply(1:50, function(r) {
    ew <- new.env(); load(file.path(DM, sprintf("divMig.%d.w.rda", r)), envir = ew)
    em <- new.env(); load(file.path(DM, sprintf("divMig.%d.matrices.rda", r)), envir = em)
    out <- list()
    for (lab in names(LAB)) {
      T0 <- true_migration(RATES[[lab]][1], RATES[[lab]][2]); Tr <- T0 / max(T0[off]); TrT <- t(Tr)
      Warr <- ew$gravity_W[[lab]]; Marr <- em$divmig_M[[lab]]
      gens <- dimnames(Warr)[[3]]
      for (j in seq_along(gens)) {
        W <- Warr[POPS, POPS, j]; M <- Marr[POPS, POPS, gens[j]]
        reason <- if (any(!is.finite(M[off]))) "nonfinite" else if (any(M[off] < 0)) "negative" else ""
        Ws <- W / max(W[off])
        v <- c(w_rho_run = rho(Ws[off], Tr[off]), w_rmse_run = rmse(Ws[off], Tr[off]),
               w_rho_tr = rho(Ws[off], TrT[off]), w_rmse_tr = rmse(Ws[off], TrT[off]),
               nm_rho_run = NA, nm_rmse_run = NA, nm_rho_tr = NA, nm_rmse_tr = NA)
        if (!nzchar(reason)) { Ms <- M / max(M[off])
          v[5:8] <- c(rho(Ms[off], Tr[off]), rmse(Ms[off], Tr[off]), rho(Ms[off], TrT[off]), rmse(Ms[off], TrT[off])) }
        out[[length(out) + 1]] <- data.frame(Replicate = r, Treatment = LAB[[lab]], Generation = as.integer(gens[j]),
                                             Exclusion = reason, t(v))
      }
    }
    do.call(rbind, out)
  }, mc.cores = CORES)
  d <- do.call(rbind, rows)
  # reproduction check against the stored w_raw metrics (all replicates)
  st <- bind_rows(lapply(1:50, function(r) { e <- new.env()
    load(file.path(DM, sprintf("divMig.%d.truth_w.rda", r)), envir = e); e$gravity_truth_w })) |>
    filter(Score == "w_raw") |> mutate(Treatment = LAB[as.character(Treatment)])
  chk <- inner_join(d, st, by = c("Replicate", "Treatment", "Generation"))
  list(d = d, repro_max_err = max(abs(chk$w_rho_run - chk$rho_all), abs(chk$w_rmse_run - chk$rmse_all)),
       repro_n = nrow(chk))
})

## ---- Task 3c: similarity vs Nm, as run and Nm transposed -----------------------
pij_similarity <- function(pij) {                       # as R/divmigrate_pij_correlation.R
  A <- pij; diag(A) <- 0
  g <- graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE)
  D <- distances(g, mode = "out", weights = E(g)$weight)[rownames(pij), colnames(pij)]
  D[!is.finite(D)] <- NA; S <- 1 / (1 + D); diag(S) <- NA; S
}
timed("simil", {
  rows <- mclapply(1:50, function(r) {
    ep <- new.env(); load(file.path(DM, sprintf("divMig.%d.pij.rda", r)), envir = ep)
    em <- new.env(); load(file.path(DM, sprintf("divMig.%d.matrices.rda", r)), envir = em)
    ec <- new.env(); load(file.path(DM, sprintf("divMig.%d.correlation.rda", r)), envir = ec)
    out <- list()
    for (lab in names(LAB)) {
      Parr <- ep$divmig_P[[lab]]; Marr <- em$divmig_M[[lab]]; gens <- dimnames(Parr)[[3]]
      for (j in seq_along(gens)) {
        S <- pij_similarity(Parr[POPS, POPS, j]); M <- Marr[POPS, POPS, gens[j]]
        off <- row(M) != col(M); M[!is.finite(M)] <- NA; Mt <- t(M)
        k1 <- off & !is.na(S) & !is.na(M); k2 <- off & !is.na(S) & !is.na(Mt)
        sa <- (S - t(S))[upper.tri(S)]; ma <- (M - t(M))[upper.tri(M)]; ka <- !is.na(sa) & !is.na(ma)
        out[[length(out) + 1]] <- data.frame(Replicate = r, Treatment = LAB[[lab]], Generation = as.integer(gens[j]),
          rho_run = rho(S[k1], M[k1]), rho_tr = rho(S[k2], Mt[k2]),
          rho_asym_run = rho(sa[ka], ma[ka]), rho_asym_tr = rho(sa[ka], -ma[ka]))
      }
    }
    d <- do.call(rbind, out)
    stored <- ec$divmig_cor |> mutate(Treatment = LAB[as.character(Treatment)]) |>
      select(Replicate, Treatment, Generation, Correlation, Correlation_asym)
    left_join(d, stored, by = c("Replicate", "Treatment", "Generation"))
  }, mc.cores = CORES)
  d <- do.call(rbind, rows)
  list(d = d, repro_max_err = max(abs(d$rho_run - d$Correlation), na.rm = TRUE),
       repro_max_err_asym = max(abs(d$rho_asym_run - d$Correlation_asym), na.rm = TRUE))
})

## ---- Delta R^2 on the symmetric gradient --------------------------------------------
timed("ibgd_sym", {
  sym <- c(`sym-mid` = 4L, `sym-low` = 5L, `sym-verylow` = 6L)
  jobs <- expand.grid(Replicate = 1:50, Treatment = names(sym), generation = seq(2004L, 2999L, 5L),
                      stringsAsFactors = FALSE)
  stats <- c("edges", "rel_fwd", "r2_cgd", "r2_pgd_blind", "r2_pgd_dir")
  res <- mclapply(seq_len(nrow(jobs)), function(i) {
    g <- gc_load_graph(jobs$Replicate[i], jobs$generation[i], sym[[jobs$Treatment[i]]])
    st <- setNames(rep(NA_real_, length(stats)), stats); status <- "missing"
    if (!is.null(g)) { o <- tryCatch(ibd_one(g), error = function(e) NULL)
      if (is.null(o)) status <- "failed" else { st <- o[stats]; status <- "ok" } }
    data.frame(jobs[i, ], t(st), status = status)
  }, mc.cores = CORES)
  list(traj = do.call(rbind, res))
})

## ---- Tasks 4-5: source-sink test (round 1, new sign) and overlap ----------------
timed("sstest", {
  e1 <- new.env(); load("data/derived/node_potential_test.rda", envir = e1)
  tt <- e1$tests_tbl |> filter(field == "phi", nodes == "all", cor == "spearman") |>
    transmute(Replicate, scenario = as.character(scenario), generation, null, r_S = -r_obs, p) # S = -phi
  ct <- e1$census_tbl |> transmute(Replicate, scenario = as.character(scenario), generation, n_comp,
                                   resid_frac = cycle_frac)
  list(tests = tt, census = ct)
})

timed("overlap", {
  ef <- new.env(); load("data/derived/ibgd_trajectory_fit.rda", envir = ef)
  ib <- bind_rows(ef$traj |> mutate(Treatment = as.character(Treatment)),
                  rean$ibgd_sym$traj) |>
    filter(status == "ok") |>
    transmute(Replicate, scenario = Treatment, generation, dR2 = r2_pgd_dir - r2_cgd, rel_fwd)
  s <- rean$sstest$tests |> filter(null == "N1w")
  d <- inner_join(s, ib, by = c("Replicate", "scenario", "generation")) |>
    mutate(block = BLOCK(generation), ss_rej = p < 0.05, ib_pos = dR2 > 0)
  list(d = d, n_ss = nrow(s), n_joined = nrow(d))
})

## ---- Task 7b: identities in the new notation ------------------------------------
timed("ident", {
  pick <- rbind(data.frame(Replicate = sample(1:50, 10), gen = 1999L, sc = NA_integer_),
                data.frame(Replicate = sample(1:50, 10), gen = 2454L, sc = 2L))
  rows <- lapply(seq_len(nrow(pick)), function(k) {
    g <- gc_load_graph(pick$Replicate[k], pick$gen[k], if (is.na(pick$sc[k])) NULL else pick$sc[k])
    cs <- gc_census(g)
    # global bandwidth: Delta_{i->j} = K_ij (1/d_j - 1/d_i), d_i = sum_k K_ik
    bg <- mean(E(g)$weight); csg <- gc_census(g, bandwidth = bg)
    Kmat <- matrix(0, 25, 25, dimnames = list(POPS, POPS)); el <- as_edgelist(g)
    Kmat[el] <- Kmat[el[, 2:1]] <- exp(-E(g)$weight^2 / (2 * bg^2))
    dsum <- rowSums(Kmat); ed <- csg$edges
    Kij <- Kmat[cbind(ed$i, ed$j)]
    idg <- Kij * (1 / dsum[ed$j] - 1 / dsum[ed$i])
    # spanning tree: gradient fraction 1
    tr <- mst(g); ctr <- gc_census(tr)
    W <- cs$W
    data.frame(rowsum_err = max(abs(rowSums(W, na.rm = TRUE)[rowSums(W > 0, na.rm = TRUE) > 0] - 1)),
               balance_err = max(abs(cs$nodes$balance - (cs$nodes$g - 1))),
               sum_balance = abs(sum(cs$nodes$g - 1)),
               global_id_err = max(abs(ed$delta - idg)),
               global_id_wrong_sign_err = max(abs(ed$delta + idg)),
               tree_grad_frac_err = abs(ctr$census[["grad_frac"]] - 1),
               telescoping_err = abs(sum(cs$B %*% cs$nodes$S) - sum(cs$nodes$S * cs$nodes$n)),
               resid_frac = cs$census[["resid_frac"]])
  })
  list(d = cbind(pick, do.call(rbind, rows)))
})

## ---- Task 8: bandwidth sensitivity of the locked test ----------------------------
BW_ALTS <- list(canonical = list(), global = "global", `local x1/2` = list(scale = 0.5),
                `local x2` = list(scale = 2), perplexity = list(perplexity = 4))
timed("bw", {
  jobs <- rbind(
    expand.grid(Replicate = 1:50, scenario = c("Isotropic", "Redistributed"), generation = seq(2254L, 2454L, 50L),
                stringsAsFactors = FALSE),
    expand.grid(Replicate = 1:50, scenario = c("Isotropic", "Obstructed"), generation = seq(2754L, 2954L, 50L),
                stringsAsFactors = FALSE))
  sc <- c(Isotropic = 1L, Redistributed = 2L, Obstructed = 3L)
  res <- mclapply(seq_len(nrow(jobs)), function(i) {
    set.seed(SEED + i)
    g <- gc_load_graph(jobs$Replicate[i], jobs$generation[i], sc[[jobs$scenario[i]]])
    if (is.null(g)) return(NULL)
    nodes <- sort(V(g)$name); x <- gc_x(nodes)
    bbar <- mean(E(g)$weight)
    Aw <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)[nodes, nodes]
    Aw[Aw > 0] <- exp(-Aw[Aw > 0]^2 / (2 * bbar^2))
    lw <- spdep::mat2listw(Aw, style = "B", zero.policy = TRUE)
    out <- list()
    for (nm in names(BW_ALTS)) {
      a <- BW_ALTS[[nm]]
      cs <- tryCatch(if (identical(a, "global")) gc_census(g, bandwidth = bbar) else do.call(gc_census, c(list(g), a)),
                     error = function(e) NULL)
      if (is.null(cs)) next
      S <- cs$nodes$S[match(nodes, cs$nodes$node)]
      r_obs <- rho(S, x)
      dr <- as.matrix(adespatial::msr(S, lw, nrepet = B, method = "singleton", simplify = TRUE))
      r_null <- as.numeric(cor(dr, x, method = "spearman"))
      out[[nm]] <- data.frame(jobs[i, ], bandwidth = nm, r_S = r_obs,
                              p = (1 + sum(abs(r_null) >= abs(r_obs) - 1e-12)) / (B + 1),
                              dbar = cs$census[["dbar"]], resid_frac = cs$census[["resid_frac"]])
    }
    do.call(rbind, out)
  }, mc.cores = CORES, mc.preschedule = FALSE)
  list(d = do.call(rbind, res))
})

rean$meta <- list(seed = SEED, B = B, gstudio = as.character(packageVersion("gstudio")),
                  adespatial = as.character(packageVersion("adespatial")),
                  igraph = as.character(packageVersion("igraph")), R = R.version.string,
                  date = as.character(Sys.Date()))
save_rean()

## ---- Summaries --------------------------------------------------------------------
if (sys.nframe() == 0L) {
  options(width = 160)
  cat("\n=== Task 3a: truth, reproduction check: max err", signif(rean$truth$repro_max_err, 3),
      "over", rean$truth$repro_n, "censuses\n")
  td <- rean$truth$d |> mutate(win = cut(Generation, c(2000, 2250, 2500, 2750, 2900, 3000)))
  print(td |> group_by(Treatment) |>
          summarise(censuses = n(), excluded = sum(Exclusion != ""), pct_excl = round(100 * mean(Exclusion != ""), 1),
                    w_rho_run = median(w_rho_run), w_rho_tr = median(w_rho_tr),
                    w_rho_run_paired = median(w_rho_run[Exclusion == ""]), w_rho_tr_paired = median(w_rho_tr[Exclusion == ""]),
                    nm_rho_run = median(nm_rho_run, na.rm = TRUE), nm_rho_tr = median(nm_rho_tr, na.rm = TRUE),
                    w_rmse_run = median(w_rmse_run), w_rmse_tr = median(w_rmse_tr),
                    nm_rmse_run = median(nm_rmse_run, na.rm = TRUE), nm_rmse_tr = median(nm_rmse_tr, na.rm = TRUE)) |>
          as.data.frame(), digits = 3)
  cat("\nBy window (median; paired = censuses where Nm is defined):\n")
  print(td |> group_by(Treatment, win) |>
          summarise(pct_retained = round(100 * mean(Exclusion == ""), 1),
                    w_tr_all = median(w_rho_tr), w_tr_paired = median(w_rho_tr[Exclusion == ""]),
                    w_run_all = median(w_rho_run), nm_run = median(nm_rho_run, na.rm = TRUE),
                    nm_tr = median(nm_rho_tr, na.rm = TRUE), .groups = "drop") |> as.data.frame(), digits = 3)
  cat("\nPaired per-census: fraction with w (transposed) rho > Nm rho, w RMSE < Nm RMSE\n")
  print(td |> filter(Exclusion == "") |> group_by(Treatment) |>
          summarise(n = n(), w_higher_rho = mean(w_rho_tr > nm_rho_run), w_lower_rmse = mean(w_rmse_tr < nm_rmse_run),
                    w_run_higher_rho = mean(w_rho_run > nm_rho_run)) |> as.data.frame(), digits = 3)
  cat("\nRetained at the final census (2999):\n")
  print(td |> filter(Generation == max(Generation)) |> group_by(Treatment) |>
          summarise(retained = mean(Exclusion == "")) |> as.data.frame())

  cat("\n=== Task 3c: similarity; reproduction max err", signif(rean$simil$repro_max_err, 3),
      "(asym", signif(rean$simil$repro_max_err_asym, 3), ")\n")
  print(rean$simil$d |> mutate(win = cut(Generation, c(2000, 2250, 2500, 2750, 3000))) |>
          group_by(Treatment, win) |>
          summarise(rho_run = median(rho_run, na.rm = TRUE), rho_tr = median(rho_tr, na.rm = TRUE),
                    asym_run = median(rho_asym_run, na.rm = TRUE), asym_tr = median(rho_asym_tr, na.rm = TRUE),
                    .groups = "drop") |> as.data.frame(), digits = 3)

  cat("\n=== Delta R^2 > 0 on the symmetric gradient (FPR) ===\n")
  print(rean$ibgd_sym$traj |> group_by(Treatment) |>
          summarise(n = n(), ok = sum(status == "ok"),
                    fpr = mean(r2_pgd_dir > r2_cgd, na.rm = TRUE),
                    lo = binom.test(sum(r2_pgd_dir > r2_cgd, na.rm = TRUE), sum(status == "ok"))$conf.int[1],
                    hi = binom.test(sum(r2_pgd_dir > r2_cgd, na.rm = TRUE), sum(status == "ok"))$conf.int[2],
                    med_dR2 = median(r2_pgd_dir - r2_cgd, na.rm = TRUE)) |> as.data.frame(), digits = 3)

  cat("\n=== Task 5: source-sink test (S, Spearman), FPR at alpha 0.01/0.05/0.10 ===\n")
  s5 <- rean$sstest$tests |> mutate(block = BLOCK(generation))
  print(s5 |> filter(scenario %in% c("Burn-in", "Isotropic", "sym-mid", "sym-low", "sym-verylow")) |>
          group_by(null, scenario) |>
          summarise(n = n(), a01 = mean(p < .01), a05 = mean(p < .05), a10 = mean(p < .10), .groups = "drop") |>
          pivot_wider(names_from = null, values_from = c(a01, a05, a10)) |> as.data.frame(), digits = 3)
  cat("\nN1w p-value histogram (10 bins), burn-in and isotropic:\n")
  print(s5 |> filter(null == "N1w", scenario %in% c("Burn-in", "Isotropic")) |>
          group_by(scenario, bin = cut(p, seq(0, 1, .1), include.lowest = TRUE)) |> summarise(f = n(), .groups = "drop") |>
          group_by(scenario) |> mutate(f = round(f / sum(f), 3)) |> pivot_wider(names_from = bin, values_from = f) |> as.data.frame())
  cat("\nPower (N1w, alpha 0.05), median r(S,x), fraction of censuses with r(S,x)<0:\n")
  print(s5 |> filter(null == "N1w", scenario %in% c("Isotropic", "Redistributed", "Obstructed")) |>
          group_by(scenario, block) |>
          summarise(n = n(), reject = mean(p < .05), med_r = median(r_S), frac_neg = mean(r_S < 0), .groups = "drop") |>
          as.data.frame(), digits = 3)
  cat("\nLineage sign consistency: fraction of lineages whose mean r(S,x) over the block is < 0:\n")
  print(s5 |> filter(null == "N1w", scenario != "Burn-in") |> group_by(scenario, block, Replicate) |>
          summarise(m = mean(r_S), .groups = "drop") |> group_by(scenario, block) |>
          summarise(lineages = n(), frac_neg = mean(m < 0), .groups = "drop") |>
          pivot_wider(names_from = block, values_from = frac_neg) |> as.data.frame(), digits = 3)
  lean <- s5 |> filter(null == "N1w") |> group_by(scenario, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |>
    group_by(scenario) |> summarise(mean = mean(m), pos = sum(m > 0), n = n(),
                                    sign_p = binom.test(sum(m > 0), n())$p.value, .groups = "drop")
  cat("\nLineage lean, r(S,x) by lineage (all censuses):\n"); print(as.data.frame(lean), digits = 3)
  cat("\nCensuses with >1 connected component:", sum(rean$sstest$census$n_comp > 1), "of", nrow(rean$sstest$census), "\n")
  cat("Non-directional residual share, median by scenario:\n")
  print(rean$sstest$census |> group_by(scenario) |>
          summarise(med = median(resid_frac), q25 = quantile(resid_frac, .25), q75 = quantile(resid_frac, .75)) |>
          as.data.frame(), digits = 3)

  cat("\n=== Task 4: overlap (joined", rean$overlap$n_joined, "of", rean$overlap$n_ss, "censuses) ===\n")
  ov <- rean$overlap$d
  print(ov |> group_by(scenario, block) |>
          summarise(n = n(), both = sum(ss_rej & ib_pos), ss_only = sum(ss_rej & !ib_pos),
                    ib_only = sum(!ss_rej & ib_pos), neither = sum(!ss_rej & !ib_pos),
                    both_of_either = both / max(1, both + ss_only + ib_only),
                    same_dir_of_both = mean((r_S < 0 & rel_fwd < 1)[ss_rej & ib_pos]),
                    .groups = "drop") |> as.data.frame(), digits = 3)
  cat("\nPooled by scenario:\n")
  print(ov |> group_by(scenario) |>
          summarise(n = n(), ss = mean(ss_rej), ib = mean(ib_pos), both = sum(ss_rej & ib_pos),
                    either = sum(ss_rej | ib_pos), both_of_either = both / either,
                    same_dir_of_both = mean((r_S < 0 & rel_fwd < 1)[ss_rej & ib_pos]),
                    phi = suppressWarnings(cor(ss_rej, ib_pos)), .groups = "drop") |> as.data.frame(), digits = 3)

  cat("\n=== Task 7b: identities (max over 20 censuses) ===\n")
  print(sapply(rean$ident$d[, -(1:3)], max), digits = 3)
  cat("resid_frac range:", range(rean$ident$d$resid_frac), "\n")

  cat("\n=== Task 8: bandwidth sensitivity (N1w, B = 999) ===\n")
  bw <- rean$bw$d |> mutate(peak = ifelse(generation < 2700, "2254-2499", "2754-2999"))
  bw_tab <- bw |> filter((scenario == "Redistributed") | (scenario == "Obstructed") | scenario == "Isotropic") |>
    group_by(scenario, peak, bandwidth) |>
    summarise(n = n(), reject = mean(p < .05), med_r = median(r_S), frac_neg = mean(r_S < 0), .groups = "drop")
  lin <- bw |> group_by(scenario, peak, bandwidth, Replicate) |> summarise(m = mean(r_S), .groups = "drop") |>
    group_by(scenario, peak, bandwidth) |> summarise(lin_neg = mean(m < 0), .groups = "drop")
  print(left_join(bw_tab, lin, by = c("scenario", "peak", "bandwidth")) |> as.data.frame(), digits = 3)
  cat("\nWall times (min):", paste(names(rean), sapply(rean, function(z) if (is.list(z) && !is.null(z$wall_min)) round(z$wall_min, 1) else NA)), "\n")
}
