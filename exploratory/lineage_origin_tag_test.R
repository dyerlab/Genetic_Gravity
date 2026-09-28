# exploratory/lineage_origin_tag_test.R
#
# Task 9 of local/node_potential_reanalysis.md: does one generation of the
# simulation driver's migrate -> mate -> resample-to-N cycle, under the
# symmetric burn-in matrix, favour parents from the upstream (i - 1) or the
# downstream (i + 1) neighbour?
#
# From saved burn-in states at generation 1999, every individual is tagged with
# its deme of origin. One generation follows AsymmetricPopGraphs/R/simulate.R
# one_generation(): gstudio::migrate() with stepping_stone_1d_asym(25, 0.025,
# 0.025), then per deme the driver's mixed_mating(N = 1, s = 0) parent draws
# (mothers and fathers each sample(seq_len(n), n, replace = TRUE)), then
# resampling of offspring to N = 100 when n != 100. Alleles are not needed to
# count parental origins, so the parent indices are drawn exactly as
# mixed_mating() draws them and the origins of both parents are recorded.
#
# Spec: 10 states, one generation each. Also run: 20 independent generations
# from each of all 50 states, for a tighter CI (same code, more draws).
#
# Run from the repository root:  Rscript exploratory/lineage_origin_tag_test.R
# Output: data/derived/lineage_origin_tag_test.rds

suppressPackageStartupMessages({ library(gstudio); library(dplyr) })
set.seed(20260930)
SIM <- "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"
K <- 25L; N <- 100L; POPS <- sprintf("Pop%02d", 1:K)

stepping_stone_1d_asym <- function(K, m_fwd, m_rev, pop_names) {   # copied from the driver
  mat <- matrix(0.0, nrow = K, ncol = K)
  for (i in seq_len(K)) {
    has_left <- i > 1L; has_right <- i < K
    if (has_left)  mat[i, i - 1L] <- m_rev
    if (has_right) mat[i, i + 1L] <- m_fwd
    mat[i, i] <- 1.0 - (has_left * m_rev) - (has_right * m_fwd)
  }
  rownames(mat) <- colnames(mat) <- pop_names
  mat
}
burn <- stepping_stone_1d_asym(K, 0.025, 0.025, POPS)

one_gen_origins <- function(pop) {
  pop <- migrate(pop, stratum = "Population", m = burn)
  out <- list()
  for (pn in POPS) {
    d <- pop[pop$Population == pn, , drop = FALSE]; n <- nrow(d)
    moms <- sample(seq_len(n), size = n, replace = TRUE)       # mixed_mating(), s = 0
    dads <- sample(seq_len(n), size = n, replace = TRUE)
    keep <- seq_len(n)
    if (n != N) keep <- sample(seq_len(n), size = N, replace = TRUE)   # resample offspring to N
    x <- as.integer(sub("\\D+", "", pn))
    o <- c(d$Origin[moms[keep]], d$Origin[dads[keep]])   # both parents of every offspring
    out[[pn]] <- data.frame(deme = x, n_post = n, imm_up = sum(d$Origin == x - 1L), imm_down = sum(d$Origin == x + 1L),
                            up = sum(o == x - 1L), down = sum(o == x + 1L),
                            resident = sum(o == x), other = sum(abs(o - x) > 1L))
  }
  do.call(rbind, out)
}

reps <- sample(1:50, 10)                       # the spec's 10 states
res <- list()
for (r in c(reps, setdiff(1:50, reps))) {       # then the other 40, for power
  e <- new.env(); load(sprintf("%s/replicate%d/rep%d-genotypes-1999.rda", SIM, r, r), envir = e)
  pop <- e$genotypes[, c("ID", "Population")]
  pop$Origin <- as.integer(sub("\\D+", "", pop$Population))
  for (k in 1:20) res[[length(res) + 1]] <- cbind(rep = r, draw = k, one_gen_origins(pop))
}
d <- do.call(rbind, res)
saveRDS(d, "data/derived/lineage_origin_tag_test.rds")

summ <- function(x, lab) {
  u <- sum(x$up); v <- sum(x$down); bt <- binom.test(u, u + v)
  data.frame(sample = lab, upstream = u, downstream = v, frac_up = u / (u + v),
             lo = bt$conf.int[1], hi = bt$conf.int[2], p = bt$p.value)
}
int <- filter(d, deme >= 2, deme <= 24)
cat("Replicates:", paste(sort(reps), collapse = ", "), "\n\n")
print(rbind(summ(filter(int, draw == 1, rep %in% reps), "spec: 10 states x 1 generation"),
            summ(filter(int, rep %in% reps), "10 states x 20 generations"),
            summ(int, "50 states x 20 generations")), digits = 4, row.names = FALSE)
cat("\nImmigrants (individuals, before mating), interior demes:\n")
print(rbind(summ(transmute(filter(int, draw == 1, rep %in% reps), up = imm_up, down = imm_down), "spec: immigrants, 10 x 1"),
            summ(transmute(int, up = imm_up, down = imm_down), "immigrants, 50 x 20")), digits = 4, row.names = FALSE)
# Gametes cluster within immigrants, so the gamete binomial is too narrow. The
# independent unit is one generation from one state (200 units): test the mean
# per-unit upstream fraction, and the per-unit gametes-per-immigrant by side.
unit <- int |> group_by(rep, draw) |>
  summarise(up = sum(up), down = sum(down), imm_up = sum(imm_up), imm_down = sum(imm_down), .groups = "drop") |>
  mutate(frac_up = up / (up + down), gpi_up = up / imm_up, gpi_down = down / imm_down)
tt <- t.test(unit$frac_up, mu = 0.5)
cat(sprintf("\nPer-generation upstream gamete fraction: mean %.4f (95%% CI %.4f-%.4f), t-test p = %.3g, n = %d\n",
            mean(unit$frac_up), tt$conf.int[1], tt$conf.int[2], tt$p.value, nrow(unit)))
st <- unit |> group_by(rep) |> summarise(f = mean(frac_up))
cat(sprintf("Per-state means (50 states): mean %.4f, %d/50 above 0.5; t-test p = %.3g\n", mean(st$f),
            sum(st$f > 0.5), t.test(st$f, mu = 0.5)$p.value))
u10 <- filter(unit, rep %in% reps)
cat(sprintf("Spec's 10 states only (200 units): mean %.4f, t-test p = %.3g\n", mean(u10$frac_up),
            t.test(u10$frac_up, mu = 0.5)$p.value))
cat(sprintf("Gametes per immigrant: upstream %.3f, downstream %.3f (paired t p = %.3g)\n",
            mean(unit$gpi_up), mean(unit$gpi_down), t.test(unit$gpi_up, unit$gpi_down, paired = TRUE)$p.value))
cat("\nPer deme (20 generations pooled, gamete binomial; too narrow, see above):\n")
pd <- int |> group_by(deme) |> group_modify(~ summ(.x, "")) |> ungroup() |> select(-sample)
print(as.data.frame(pd), digits = 3)
cat(sprintf("\nDemes with p < 0.05: %d of %d (expected ~%.1f by chance)\n", sum(pd$p < .05), nrow(pd), 0.05 * nrow(pd)))
cat("Post-migration deme size: mean", round(mean(d$n_post), 2), "; fraction != 100:", round(mean(d$n_post != 100), 3), "\n")
