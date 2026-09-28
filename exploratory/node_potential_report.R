# exploratory/node_potential_report.R
# Tables for Steps 3-5 of local/node_potential_test_instructions.md, from
# data/derived/node_potential_test.rda (exploratory/node_potential_test.R).
# Run from the repository root:  Rscript exploratory/node_potential_report.R

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
load("data/derived/node_potential_test.rda")
nulls <- c("N0", "N1", "N1w", "N1true")
tt <- tests_tbl |> mutate(null = factor(null, nulls), field = factor(field, c("phi", "div", "phi_res")))
fmt <- function(d) print(as.data.frame(d), digits = 3, row.names = FALSE)

cat("meta:", unlist(meta), "\n")

## ---- Step 3: calibration --------------------------------------------------------
cat("\n== Step 3a: FPR, burn-in 1904-1999 (1,000 censuses), all nodes ==\n")
cal <- tt |> filter(scenario == "Burn-in")
fmt(cal |> filter(nodes == "all") |> group_by(field, cor, null) |>
      summarise(a01 = mean(p <= .01), a05 = mean(p <= .05), a10 = mean(p <= .10), .groups = "drop") |>
      pivot_wider(names_from = null, values_from = c(a01, a05, a10), names_glue = "{null}_{.value}") |>
      select(field, cor, starts_with("N0"), starts_with("N1_"), starts_with("N1w"), starts_with("N1true")))

cat("\n== Step 3b: p-value histogram, burn-in, Spearman, all nodes (uniform = 0.10) ==\n")
fmt(cal |> filter(nodes == "all", cor == "spearman") |> group_by(field, null) |>
      summarise(h = list(as.numeric(table(cut(p, seq(0, 1, .1), include.lowest = TRUE)) / n())), .groups = "drop") |>
      unnest_wider(h, names_sep = "_"))

cat("\n== Step 3c: FPR at 0.05, burn-in, interior nodes only ==\n")
fmt(cal |> filter(nodes == "interior") |> group_by(field, cor, null) |> summarise(a05 = mean(p <= .05), .groups = "drop") |>
      pivot_wider(names_from = null, values_from = a05))

cat("\n== Step 3d: FPR at 0.05 by graph-size tercile, burn-in, Spearman, all nodes ==\n")
cal_c <- cal |> left_join(census_tbl |> select(Replicate, scenario, generation, edges),
                          by = c("Replicate", "scenario", "generation"))
cal_c$size <- cut(cal_c$edges, quantile(census_tbl$edges[census_tbl$scenario == "Burn-in"], c(0, 1/3, 2/3, 1)),
                  include.lowest = TRUE)
fmt(cal_c |> filter(nodes == "all", cor == "spearman") |> group_by(field, size, null) |>
      summarise(a05 = mean(p <= .05), .groups = "drop") |> pivot_wider(names_from = null, values_from = a05))

cat("\n== Step 3e: per-replicate FPR at 0.05 (burn-in, Spearman, all nodes): quantiles ==\n")
fmt(cal |> filter(nodes == "all", cor == "spearman") |> group_by(field, null, Replicate) |>
      summarise(f = mean(p <= .05), .groups = "drop") |> group_by(field, null) |>
      summarise(min = min(f), q25 = quantile(f, .25), med = median(f), q75 = quantile(f, .75),
                max = max(f), reps_zero = sum(f == 0), .groups = "drop"))

cat("\n== Step 3f: FPR at 0.05 under symmetric scenarios (forward phase; Spearman, all nodes) ==\n")
sym <- c("Burn-in", "Isotropic", "sym-mid", "sym-low", "sym-verylow")
fmt(tt |> filter(scenario %in% sym, nodes == "all", cor == "spearman") |>
      mutate(scenario = factor(scenario, sym)) |> group_by(field, scenario, null) |>
      summarise(a05 = mean(p <= .05), .groups = "drop") |> pivot_wider(names_from = null, values_from = a05))

## ---- Step 4: power and consistency ----------------------------------------------
fw <- tt |> filter(scenario %in% c("Isotropic", "Redistributed", "Obstructed")) |>
  mutate(block = cut(generation, c(2000, 2250, 2500, 2750, 3000),
                     labels = c("2004-2249", "2254-2499", "2504-2749", "2754-2999")))
cat("\n== Step 4a: rejection rate at 0.05 (Spearman, all nodes) ==\n")
fmt(fw |> filter(nodes == "all", cor == "spearman") |> group_by(field, scenario, block, null) |>
      summarise(rej = mean(p <= .05), .groups = "drop") |> pivot_wider(names_from = null, values_from = rej))
cat("\n== Step 4b: median r and sign consistency (Spearman, all nodes; one r per census) ==\n")
fmt(fw |> filter(nodes == "all", cor == "spearman", null == "N0") |> group_by(field, scenario, block) |>
      summarise(med_r = median(r_obs), frac_pos = mean(r_obs > 0), .groups = "drop"))

cat("\n== Step 4c: across the Redistributed reversal: node-level r vs chain-oriented sum, by generation ==\n")
load_nodes <- readRDS("data/derived/node_potential_nodes.rds")
# chain-oriented A reconstructed from node fields: A_grad = sum phi_i n_i is not stored; use
# div-weighted: A = sum_edges o*Delta. Here we report the gradient telescoping proxy
# A_ends = contribution of chain ends, via phi at x = 1 and x = 25.
rev <- fw |> filter(scenario == "Redistributed", nodes == "all", cor == "spearman", null == "N1", field == "phi") |>
  group_by(generation) |> summarise(med_r = median(r_obs), frac_pos = mean(r_obs > 0), rej = mean(p <= .05))
ends <- load_nodes |> filter(scenario == "Redistributed") |> group_by(Replicate, generation) |>
  summarise(phi_end1 = phi[x == 1][1], phi_end25 = phi[x == 25][1],
            r_int = cor(phi[x >= 3 & x <= 23], x[x >= 3 & x <= 23], method = "spearman"), .groups = "drop") |>
  group_by(generation) |> summarise(med_phi_end1 = median(phi_end1, na.rm = TRUE),
                                    med_phi_end25 = median(phi_end25, na.rm = TRUE),
                                    med_r_interior = median(r_int, na.rm = TRUE))
fmt(left_join(rev, ends, by = "generation"))

## ---- Step 5: diagnostics ------------------------------------------------------------
cat("\n== Step 5: census diagnostics (medians) ==\n")
fmt(census_tbl |> group_by(scenario) |>
      summarise(n = n(), edges = median(edges), disconnected = mean(n_comp > 1),
                cor_phi_phi0 = median(cor_phi_phi0, na.rm = TRUE), grad_local = median(grad_frac),
                grad_global = median(grad_frac_global, na.rm = TRUE), cycle_local = median(cycle_frac),
                has_leaf = mean(has_leaf), argmax_leaf = mean(argmax_leaf), argmin_leaf = mean(argmin_leaf),
                .groups = "drop"))
