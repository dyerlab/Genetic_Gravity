# exploratory/node_potential_step0.R
#
# Step 0 of local/node_potential_test_instructions.md: identity checks on 20
# burn-in censuses and 20 forward (Redistributed) censuses, plus the
# telescoping decomposition of the chain-oriented sum.
#
# Run from the repository root:  Rscript exploratory/node_potential_step0.R

source("exploratory/node_potential_helpers.R")
set.seed(20260928)

cases <- rbind(data.frame(rep = sample(1:50, 20), gen = sample(seq(1904L, 1999L, 5L), 20, TRUE),
                          scenario = NA_integer_, set = "burn-in"),
               data.frame(rep = sample(1:50, 20), gen = sample(seq(2254L, 2754L, 50L), 20, TRUE),
                          scenario = 2L, set = "Redistributed"))

one <- function(i) {
  cs <- cases[i, ]
  g <- load_graph(cs$rep, cs$gen, if (is.na(cs$scenario)) NULL else cs$scenario)
  nf <- node_fields(g); nd <- nf$nodes

  ## 1-2: kernel row sums and div = 1 - c
  id1 <- max(abs(nd$rowsum - 1)); id2 <- max(abs(nd$div - (1 - nd$c)))

  ## 3: global bandwidth: Delta = K (1/d_u - 1/d_v); conductance-weighted fit exact
  b <- mean(E(g)$weight)
  ga <- graph_asymmetries(g, bandwidth = b)
  el <- as_edgelist(ga, names = TRUE); nodes <- V(ga)$name
  K <- exp(-E(ga)$weight^2 / (2 * b^2))
  d <- as.numeric(tapply(c(K, K), factor(match(c(el[, 1], el[, 2]), nodes), seq_along(nodes)), sum))
  u <- match(el[, 1], nodes); v <- match(el[, 2], nodes)
  id3a <- max(abs(E(ga)$delta - K * (1 / d[u] - 1 / d[v])))
  Bm <- incidence(el, nodes); D <- K * Bm
  psi <- as.numeric(MASS::ginv(D) %*% E(ga)$delta)
  id3b <- max(abs(E(ga)$delta - D %*% psi))
  comp <- components(ga)$membership
  id3c <- max(tapply(psi - 1 / d, comp, function(z) diff(range(z))))   # psi = 1/d + const

  ## 4: spanning tree gradient fraction = 1 (local and global bandwidth)
  tr <- mst(g, weights = E(g)$weight)
  id4a <- abs(node_fields(tr)$census[["grad_frac"]] - 1)
  id4b <- abs(node_fields(tr, bandwidth = mean(E(tr)$weight))$census[["grad_frac"]] - 1)

  ## 5: telescoping. A = sum over edges of Delta oriented toward the higher deme
  x <- nd$x; xe_u <- x[match(nf$el[, 1], nd$node)]; xe_v <- x[match(nf$el[, 2], nd$node)]
  o <- ifelse(xe_u < xe_v, 1, -1)
  A <- sum(o * nf$delta)
  A_grad_edges <- sum(o * as.numeric(nf$B %*% nd$phi))
  n_i <- as.numeric(crossprod(nf$B, o))            # (# edges as lower) - (# edges as higher)
  A_grad_nodes <- sum(nd$phi * n_i)
  ends <- nd$x %in% c(1, 25)
  contrib <- nd$phi * n_i
  c(id1 = id1, id2 = id2, id3a = id3a, id3b = id3b, id3c = id3c, id4a = id4a, id4b = id4b,
    id5 = abs(A_grad_edges - A_grad_nodes), A = A, A_grad = A_grad_nodes,
    ends_share = sum(contrib[ends]) / A_grad_nodes,
    ends_abs_share = sum(abs(contrib[ends])) / sum(abs(contrib)),
    ends_nweight = sum(abs(n_i[ends])) / sum(abs(n_i)),
    edges = gsize(g), grad_frac = nf$census[["grad_frac"]])
}
res <- cbind(cases, do.call(rbind, parallel::mclapply(seq_len(nrow(cases)), one, mc.cores = 13)))

tol <- c(id1 = 1e-10, id2 = 1e-10, id3a = 1e-10, id3b = 1e-8, id3c = 1e-8, id4a = 1e-8, id4b = 1e-8, id5 = 1e-10)
cat("Identity checks (max over 40 censuses; tolerance):\n")
for (nm in names(tol)) cat(sprintf("  %-5s max %.2e  tol %.0e  %s\n", nm, max(res[[nm]]), tol[[nm]],
                                   if (max(res[[nm]]) <= tol[[nm]]) "PASS" else "FAIL"))
cat("\nTelescoping: A, gradient part, and chain-end share of the gradient part (median [range]):\n")
print(res |> dplyr::group_by(set) |> dplyr::summarise(
  A = median(A), A_grad = median(A_grad), grad_share_of_A = median(A_grad / A),
  ends_share = median(ends_share), ends_abs_share = median(ends_abs_share),
  ends_nweight = median(ends_nweight), grad_frac = median(grad_frac)) |> as.data.frame(), digits = 3)
saveRDS(res, "data/derived/node_potential_step0.rds")
