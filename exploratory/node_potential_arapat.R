# exploratory/node_potential_arapat.R
# Step 6 of local/node_potential_test_instructions.md: illustrative node-potential
# test on gstudio::arapat (Araptus attenuatus, Baja California). Population Graph
# built as in the gstudio vignettes; covariates latitude and longitude (population
# means). Illustrative only.
# Run from the repository root:  Rscript exploratory/node_potential_arapat.R

source("exploratory/node_potential_helpers.R")
suppressPackageStartupMessages({ library(adespatial); library(dplyr) })
set.seed(20260928)
data(arapat, package = "gstudio")
g <- popgraph(to_mv(arapat), groups = arapat$Population)
cat(sprintf("arapat Population Graph: %d nodes, %d edges, %d component(s)\n",
            vcount(g), ecount(g), components(g)$no))
nf <- node_fields(g); nd <- nf$nodes
xy <- arapat |> group_by(Population) |> summarise(lat = mean(Latitude), lon = mean(Longitude))
nd <- left_join(nd, xy, by = c(node = "Population"))
A <- as_adjacency_matrix(nf$ga, sparse = FALSE)[nd$node, nd$node]
Aw <- as_adjacency_matrix(nf$ga, attr = "weight", sparse = FALSE)[nd$node, nd$node]
bbar <- mean(E(nf$ga)$weight); Aw[Aw > 0] <- exp(-Aw[Aw > 0]^2 / (2 * bbar^2))
listws <- list(N1 = as_listw(A), N1w = as_listw(Aw))
cat(sprintf("gradient fraction %.3f; cor(phi, phi0) %.3f; leaves %d\n",
            nf$census[["grad_frac"]], nf$census[["cor_phi_phi0"]], sum(nd$leaf)))
out <- list()
for (fld in c("phi", "div", "phi_res")) {
  dr <- null_draws(nd[[fld]], listws, 999L)
  for (cv in c("lat", "lon")) out[[length(out) + 1]] <- cbind(field = fld, covariate = cv,
                                                              field_tests(nd[[fld]], nd[[cv]], dr))
}
print(do.call(rbind, out), digits = 3, row.names = FALSE)
saveRDS(list(nodes = nd, tests = do.call(rbind, out)), "data/derived/node_potential_arapat.rds")
