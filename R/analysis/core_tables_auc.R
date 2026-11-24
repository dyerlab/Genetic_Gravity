# R/analysis/core_tables_auc.R
#
# (1) Tab_TestOverlap and Tab_BandwidthSensitivity on phase cores (middle 50% of each
#     lineage's phase, F_4.5 horizon; definitions in R/analysis/phase_core_rerun.R), gamma = 1/2.
# (2) Slope-ratio AUC at gamma = 1/2 and gamma = 1: the rank-based separation of the forward-to-reverse
#     slope ratio b_fwd/b_rev (lower = forward) between asymmetric censuses and every Symmetric
#     forward-phase census, AUC = P(ratio_asym < ratio_sym) + 0.5 P(tie), 5-generation grid,
#     by scenario and phase (full phase and core). The original 0.92-0.99 (gamma = 1) was computed
#     by a script no longer in the repository; the gamma = 1 column checks the reconstruction.
#
# Run from the repository root:  Rscript R/analysis/core_tables_auc.R
# Output: data/derived/core_tables_auc.rda (printed tables)

suppressPackageStartupMessages(library(dplyr))
src <- readLines("R/analysis/phase_core_rerun.R")
eval(parse(text = src[1:(grep("^## ---- summaries", src) - 1)]))          # data, boundaries, intervals(), assign_phase()
IV <- intervals("R_gen", "F_gen")

## ---- Tab_TestOverlap, cores ------------------------------------------------------------------
ov <- inner_join(T5 |> select(sample, Replicate, generation, r_S, p), G5 |> select(sample, Replicate, generation, r2_cgd, r2_pgd_dir, rel_fwd),
                 by = c("sample", "Replicate", "generation")) |> filter(!is.na(r2_pgd_dir), sample %in% AS)
ovc <- assign_phase(ov, IV, core = TRUE) |> filter(!is.na(phase)) |> mutate(ss = p < .05, ib = r2_pgd_dir > r2_cgd)
tov <- ovc |> group_by(sample, phase) |> summarise(cens = n(), lin = n_distinct(Replicate), both = sum(ss & ib), sso = sum(ss & !ib), ibo = sum(!ss & ib),
                                                   neither = sum(!ss & !ib), boe = both / max(1, both + sso + ibo), same = sum((r_S < 0 & rel_fwd < 1)[ss & ib]), .groups = "drop")
print(as.data.frame(tov), digits = 3)
jt <- tov |> filter(!(sample == "Obstructed" & phase == "decay")); cat(sprintf("joint agreeing (excl. Obstructed decay): %d of %d\n", sum(jt$same), sum(jt$both)))

## ---- Tab_BandwidthSensitivity, cores -----------------------------------------------------------
eb <- assign_phase(BW$edge_bw |> rename(sample = Scenario), IV, core = TRUE) |> filter(!is.na(phase))
tbw <- eb |> group_by(bandwidth, phase) |> summarise(cens = n(), concord = mean(concord), rs = mean(rho_signed, na.rm = TRUE),
                                                     ra = mean(rho_abs, na.rm = TRUE), kept = mean(dbar_kept), .groups = "drop")
print(as.data.frame(tbw), digits = 3)
allc <- eb |> distinct(sample, Replicate, generation, phase); cat("core censuses by phase:\n"); print(table(allc$phase)); print(table(allc$sample, allc$phase))
pf <- anti_join(allc, eb |> filter(bandwidth == "perplexity") |> distinct(sample, Replicate, generation), by = c("sample", "Replicate", "generation"))
cat("perplexity failures in cores:\n"); print(table(pf$phase))

## ---- slope-ratio AUC ------------------------------------------------------------------------------
auc <- function(pos, null) { pos <- pos[is.finite(pos)]; null <- null[is.finite(null)]
  r <- rank(c(pos, null)); n1 <- length(pos); n0 <- length(null)
  1 - (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0) }          # P(pos < null) + 0.5 P(tie)
AUC <- bind_rows(lapply(list(`γ = ½` = G5, `γ = 1` = G1), function(G) {
  g <- G |> filter(!is.na(rel_fwd)); sym <- g$rel_fwd[g$sample == "Symmetric"]
  a <- g |> filter(sample %in% AS)
  bind_rows(lapply(c(full = FALSE, core = TRUE), function(cr) assign_phase(a, IV, core = cr) |> filter(!is.na(phase)) |>
    group_by(sample, phase) |> summarise(n = n(), auc = auc(rel_fwd, sym), .groups = "drop") |> mutate(window = if (cr) "core" else "full")))
}), .id = "bandwidth")
print(as.data.frame(AUC |> arrange(bandwidth, window, sample, phase)), digits = 3)
save(tov, tbw, AUC, file = "data/derived/core_tables_auc.rda")
