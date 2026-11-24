# R/analysis/permutation_nulls_1145_report.R
#
# Tables for Part 1d of local/gstudio_1145_verification_and_gamma_rule_instructions.md
# from data/derived/permutation_nulls_1145.rda, with an edge-by-edge comparison to the
# earlier window runs in the private research repo (uncentred, 1.14.4; centred,
# 1.14.4 + the 1.14.5 patch). Writes local/permutation_nulls_1145_tables.md.
#
# Run from the repository root:  Rscript R/analysis/permutation_nulls_1145_report.R

dir.create("output", showWarnings = FALSE)   # run reports (not tracked)
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
source("R/common/lineage_phases.R")
load("data/derived/permutation_nulls_1145.rda")                 # rew, loc_edges, loc_cens, checks, meta
FS <- file.path(Sys.getenv("GG_SIM_DIR", "/Users/rodney/Documents/Research/PopulationGenetics/AsymmetricPopGraphs/data"), "divMigrate", "firstsig")
OUTMD <- "output/permutation_nulls_1145_tables.md"
md <- function(df) {
  df <- as.data.frame(df); hdr <- paste("|", paste(names(df), collapse = " | "), "|")
  sep <- paste("|", paste(c(":----", rep(":---:", ncol(df) - 1)), collapse = " | "), "|")
  c(hdr, sep, apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")), "")
}
pct <- function(v, dig = 1) ifelse(is.na(v), "—", sprintf(paste0("%.", dig, "f%%"), 100 * v))
L <- character(0); add <- function(...) L <<- c(L, ...)

## ---- rewiring null ---------------------------------------------------------------------------
PH_LAB <- c(Reorganize = "ascent", Plateau = "plateau", `Fall apart` = "decay")
rw <- rew |> filter(!failed) |>
  mutate(phase = ifelse(sample %in% c("Redistributed", "Obstructed"), PH_LAB[as.character(lp_phase(sample, Replicate, generation))], "all"),
         col = ifelse(phase == "all", sample, paste(sample, phase)),
         band = cut(edges, c(0, 60, 80, 100, Inf), labels = c("≤ 60", "61–80", "81–100", "> 100")))
add("## Part 1d. Supplement \"Permutation nulls for Δ\" under gstudio 1.14.5", "",
    sprintf("Mirror checks (same seed, mirror vs the gstudio function): largest |Δp| = %g, largest |Δ statistic| = %g over %d censuses.",
            max(checks$p_diff), max(checks$stat_diff, na.rm = TRUE), nrow(checks)), "",
    "### Graph-level rewiring null (asymmetry_network, rewire = \"degree\", 99 draws)", "",
    sprintf("%d censuses (%d failed outright).", nrow(rew), sum(rew$failed)), "",
    "**Valid null draws by edge count (25 nodes)**", "",
    md(rw |> group_by(band) |> summarise(censuses = n(), `valid draws` = pct(sum(n_null) / (99 * n())),
                                         `censuses with 0 valid` = sum(n_null == 0), .groups = "drop") |> rename(Edges = band)))
COLS <- c("Burn-in", "Symmetric", "Redistributed ascent", "Redistributed plateau", "Redistributed decay", "Obstructed ascent", "Obstructed plateau", "Obstructed decay")
rz <- rw |> mutate(col = factor(col, COLS)) |> group_by(col) |>
  summarise(n = n(), edges = median(edges), rej = pct(mean(p < .05)), below = pct(mean(obs < null_mean)),
            zq = sprintf("%.2f to %.2f", quantile(z, .05), quantile(z, .95)), z = sprintf("%.2f", median(z)),
            obs = sprintf("%.3f", median(obs)), nm = sprintf("%.3f", median(null_mean)), minp = sprintf("%.2f", min(p)), .groups = "drop")
add("**Position of the observed mean |Δ| in the null, and rejection (one-sided, α = 0.05)**", "",
    md(rz |> rename(`Sample / phase` = col, Censuses = n, `median edges` = edges, `rejected` = rej, `obs below null mean` = below,
                    `median z` = z, `z 5–95%` = zq, `median obs` = obs, `median null mean` = nm, `min p` = minp)))
zr <- rw |> group_by(col) |> summarise(z = median(z), .groups = "drop")

## ---- label-permutation null ------------------------------------------------------------------
le <- loc_edges |> mutate(d0 = 1 / k_from - 1 / k_to, ad0 = abs(1 / k_from - 1 / k_to),
                          bin = cut(ad0, c(-Inf, 0, 0.02, 0.05, 0.1, Inf), labels = c("0", "(0, 0.02]", "(0.02, 0.05]", "(0.05, 0.1]", "> 0.1")))
old_u <- local({ e <- new.env(); load(file.path(FS, "gravity.location_window.burnin_1900_2000.rda"), envir = e); e$window_edges })
old_c <- local({ e <- new.env(); load(file.path(FS, "gravity.location_window.burnin_1900_2000_centred_gstudio1.14.5.rda"), envir = e); e$window_edges })
k3 <- c("Replicate", "generation", "from", "to")
cmp <- le |> select(all_of(k3), delta, p_centred, p_uncentred) |>
  inner_join(old_u |> select(all_of(k3), delta_u = delta, p_u = p_value), by = k3) |>
  inner_join(old_c |> select(all_of(k3), delta_c = delta, p_c = p_value), by = k3)
add("### Edge-level label-permutation null (asymmetry_permutation, 499 permutations, original seeds)", "",
    sprintf("%d censuses, %d edges (%d failed censuses). Edge-by-edge against the earlier runs (%d edges matched): max |ΔΔ| = %g; centred p vs the 1.14.4 + patch run: %d edges differ (max %g); uncentred p (same draws) vs the 1.14.4 run: %d differ (max %g).",
            nrow(loc_cens), nrow(le), 1000 - nrow(loc_cens), nrow(cmp), max(abs(cmp$delta - cmp$delta_u), abs(cmp$delta - cmp$delta_c)),
            sum(abs(cmp$p_centred - cmp$p_c) > 1e-12), max(abs(cmp$p_centred - cmp$p_c)),
            sum(abs(cmp$p_uncentred - cmp$p_u) > 1e-12), max(abs(cmp$p_uncentred - cmp$p_u))), "")
fpr <-data.frame(α = c(0.002, 0.01, 0.05, 0.10), `Uncentred (1.14.4)` = pct(sapply(c(0.002, 0.01, 0.05, 0.10), function(a) mean(le$p_uncentred <= a + 1e-9))),
                  `Centred (1.14.5)` = pct(sapply(c(0.002, 0.01, 0.05, 0.10), function(a) mean(le$p_centred <= a + 1e-9))), check.names = FALSE)
add("**Edge false-positive rate (symmetric migration; every rejection is false; rejection = p ≤ α, the convention of the supplement's 8.4% / 24.6%)**", "", md(fpr))
h <- function(p) pct(as.numeric(table(cut(p, seq(0, 1, .1), include.lowest = TRUE))) / length(p))
hist_tab <- rbind(c("Uncentred", h(le$p_uncentred)), c("Centred (1.14.5)", h(le$p_centred)))
colnames(hist_tab) <- c("p bin", sprintf("%.1f–%.1f", seq(0, .9, .1), seq(.1, 1, .1)))
add("**p-value histogram (uniform = 10% per bin)**", "", md(hist_tab))
bh <- le |> group_by(Replicate, generation) |> summarise(bh_c = any(p.adjust(p_centred, "BH") <= .05), bh_u = any(p.adjust(p_uncentred, "BH") <= .05), .groups = "drop")
byrep <- le |> group_by(Replicate) |> summarise(f = mean(p_centred <= .05 + 1e-9))
add(sprintf("Censuses with ≥ 1 edge significant after BH (q = 0.05): uncentred %s, centred %s. Per-replicate centred FPR at 0.05: %.2f–%.2f (median %.2f).",
            pct(mean(bh$bh_u)), pct(mean(bh$bh_c)), min(byrep$f), max(byrep$f), median(byrep$f)), "")
add("**By degree difference |1/kᵢ − 1/kⱼ|**", "",
    md(le |> group_by(bin) |> summarise(Edges = n(), `FPR@0.05 uncentred` = pct(mean(p_uncentred <= .05 + 1e-9), 1), `FPR@0.05 centred` = pct(mean(p_centred <= .05 + 1e-9), 1),
                                        `p > 0.9 uncentred` = pct(mean(p_uncentred > .9), 1), `p > 0.9 centred` = pct(mean(p_centred > .9), 1), .groups = "drop") |>
         rename(`|1/kᵢ − 1/kⱼ|` = bin)))
sl_obs <- coef(lm(delta ~ d0, le))[2]; sl_null <- coef(lm(null_mean ~ d0, le))[2]
ce <- le |> mutate(ends = pmin(as.integer(sub("\\D+", "", from)), as.integer(sub("\\D+", "", to))) <= 3 |
                     pmax(as.integer(sub("\\D+", "", from)), as.integer(sub("\\D+", "", to))) >= 23)
add(sprintf("cor(Δ, 1/kᵢ − 1/kⱼ) = %.2f; slope of observed Δ on the degree term = %.2f; slope of the null mean on it = %.2f (ratio %.2f). Chain-end edges (an endpoint within 2 of either end) centred FPR@0.05 %s vs interior %s.",
            cor(le$delta, le$d0), sl_obs, sl_null, sl_obs / sl_null, pct(mean(ce$p_centred[ce$ends] <= .05 + 1e-9), 0), pct(mean(ce$p_centred[!ce$ends] <= .05 + 1e-9), 0)), "",
    "(gstudio's delta for an edge stored as (i, j) is −Δ_{i→j}, so the degree term in its sign is 1/kᵢ − 1/kⱼ; the slopes are sign-invariant.)", "")

## ---- current text vs 1.14.5 ------------------------------------------------------------------
pl <- rw |> filter(col %in% c("Symmetric", "Redistributed plateau", "Obstructed plateau")) |> group_by(col) |> summarise(r = mean(p < .05), .groups = "drop")
zz <- zr |> filter(col %in% c("Burn-in", "Symmetric", "Redistributed ascent", "Redistributed plateau", "Obstructed ascent", "Obstructed plateau"))
zdec <- zr$z[zr$col == "Redistributed decay"]; rdec <- mean(rw$p[rw$col == "Redistributed decay"] < .05)
eq0 <- le |> filter(ad0 == 0)
cmp_tab <- data.frame(
  Item = c("Symmetric reference: censuses / edges",
           "Rewiring: valid draws above 100 edges",
           "Rewiring: rejections (α = 0.05), Symmetric / Redistributed plateau / Obstructed plateau (Redistributed decay)",
           "Rewiring: rejections, burn-in",
           "Rewiring: observed position relative to the null mean (median z by sample/phase)",
           "Label: edge FPR at α = 0.05, uncentred",
           "Label: excess of p near 0 and near 1 (uncentred; share in 0–0.1 / 0.9–1)",
           "Label: centred FPR at α = 0.05",
           "Label: centring removes the excess near 1 (centred share in 0.9–1)",
           "Label: observed Δ as a fraction of the degree term (slope ratio obs / null centre)",
           "Label: centred FPR on equal-degree edges",
           "Sign-exchangeability null (81% / 25% / 14–100% / 24%)",
           "Moran spectral randomization FPR on S (4–6%; true-landscape 11%)"),
  Current = c("1,000 / 95,903", "0.3%", "never rejected", "never rejected", "2.4–3.4 SD below", "8.4%", "yes", "24.6%", "yes",
              "about 0.65", "15.7%", "as stated", "as stated"),
  `gstudio 1.14.5` = c(sprintf("%s / %s", format(nrow(loc_cens), big.mark = ","), format(nrow(le), big.mark = ",")),
                       pct(with(rw[rw$edges > 100, ], sum(n_null) / (99 * length(n_null)))),
                       sprintf("%s (%s)", paste(pct(pl$r[pl$col == "Symmetric"]), pct(pl$r[pl$col == "Redistributed plateau"]), pct(pl$r[pl$col == "Obstructed plateau"]), sep = " / "), pct(rdec)),
                       pct(mean(rw$p[rw$sample == "Burn-in"] < .05)),
                       sprintf("%.1f–%.1f SD below (burn-in, Symmetric, ascent, plateau); Redistributed decay %.1f SD below", -max(zz$z), -min(zz$z), -zdec),
                       pct(mean(le$p_uncentred <= .05 + 1e-9)),
                       sprintf("yes (%s / %s)", pct(mean(le$p_uncentred <= .1)), pct(mean(le$p_uncentred > .9))),
                       pct(mean(le$p_centred <= .05 + 1e-9)),
                       sprintf("yes (%s)", pct(mean(le$p_centred > .9))),
                       sprintf("%.2f", sl_obs / sl_null),
                       sprintf("%s with p < 0.05 (the text's convention for this number); %s with p ≤ 0.05", pct(mean(eq0$p_centred < .05)), pct(mean(eq0$p_centred <= .05 + 1e-9))),
                       "not affected: uses graph_asymmetries() only (unchanged)",
                       "not affected: no 1.14.5 code path (graph_asymmetries(), adespatial)"),
  check.names = FALSE)
add("### Current supplement text vs gstudio 1.14.5", "", md(cmp_tab))
writeLines(L, OUTMD); cat("Wrote", OUTMD, "\n")
