#' Computer Note (Heredity, gstudio 1.16): worked example with the Araptus attenuatus data
#' distributed with gstudio as `arapat` (Garrick et al. 2013).
#'
#' Every number, table and figure in the note's worked example comes from this
#' script.  It uses all 363 beetles, 39 sites and three lineages (Cape,
#' Mainland, Peninsula).  The directional analyses are shown descriptively on
#' the full graph; the hypothesis tests for the peninsular expansions are those
#' of Dyer (in review), R/case_study_arapat.R, and are not repeated here.
#'
#' Output: computer_note/data/*.csv, computer_note/media/*.png and
#' computer_note/data/note_arapat_log.txt (all values quoted in the text).
#' Run from the computer_note/ folder:  Rscript R/note_arapat.R

suppressPackageStartupMessages({
  library(gstudio)
  library(dplyr)
  library(igraph)
  library(ggplot2)
})

set.seed(2009)
dir.create("data", showWarnings = FALSE)
dir.create("media", showWarnings = FALSE)
log_file <- "data/note_arapat_log.txt"
sink(log_file, split = TRUE)
cat("gstudio", as.character(packageVersion("gstudio")),
    "| igraph", as.character(packageVersion("igraph")),
    "| ggplot2", as.character(packageVersion("ggplot2")),
    "|", R.version.string, "\n\n")

# ---- 1. Data and the locus class -------------------------------------------
data(arapat)
loci <- column_class(arapat, "locus")
cat(sprintf("arapat: %d individuals, %d sites, %d loci (%s)\n",
            nrow(arapat), n_distinct(arapat$Population), length(loci),
            paste(loci, collapse = ", ")))
print(table(arapat$Species))
print(arapat |>
        group_by(Species) |>
        summarise(sites = n_distinct(Population), n = n(),
                  lat_min = min(Latitude), lat_max = max(Latitude)))
# A locus is a column of the data frame, so ordinary data-frame tools apply.
cat("\nfirst genotype at LTRS:", format(arapat$LTRS[1]),
    "| heterozygote?", is_heterozygote(arapat$LTRS[1]), "\n")
cat("missing genotypes per locus:\n")
print(sapply(arapat[loci], function(l) sum(is.na(l))))

# ---- 2. Genetic diversity by lineage (Table 2) -------------------------------
div <- arapat |>
  group_by(Species) |>
  genetic_diversity(mode = c("A", "Ae", "He", "Ho", "Fis"))
div_tab <- div |>
  group_by(Stratum) |>
  summarise(Sites = NA_integer_,
            N = NA_integer_,
            Polymorphic = sum(A > 1),
            A = mean(A), Ae = mean(Ae), He = mean(He), Ho = mean(Ho),
            Fis = 1 - mean(Ho) / mean(He)) |>
  rename(Lineage = Stratum)
sz <- arapat |> group_by(Species) |>
  summarise(Sites = n_distinct(Population), N = n())
div_tab$Sites <- sz$Sites[match(div_tab$Lineage, sz$Species)]
div_tab$N     <- sz$N[match(div_tab$Lineage, sz$Species)]
cat("\nTable 2: multilocus diversity by lineage (means over 8 loci;",
    "Fis from the means)\n")
print(as.data.frame(div_tab), digits = 3, row.names = FALSE)
write.csv(div, "data/table2_diversity_by_locus.csv", row.names = FALSE)
write.csv(div_tab, "data/table2_diversity_by_lineage.csv", row.names = FALSE)

# ---- 3. Differentiation among lineages (Table 3) ------------------------------
modes <- c("Gst", "Gst_prime", "Dest")
str_tab <- lapply(modes, function(m) {
  s <- genetic_structure(arapat, stratum = "Species", mode = m)
  data.frame(Locus = s$Locus, value = s[[2]], mode = m)
}) |>
  bind_rows() |>
  tidyr::pivot_wider(names_from = mode, values_from = value)
cat("\nTable 3: differentiation among the three lineages\n")
print(as.data.frame(str_tab), digits = 3, row.names = FALSE)
write.csv(str_tab, "data/table3_structure.csv", row.names = FALSE)

pw <- genetic_structure(arapat, stratum = "Species", mode = "Gst",
                        pairwise = TRUE)
cat("\npairwise Gst between lineages\n")
print(round(pw, 3))

perm <- genetic_structure(arapat, stratum = "Species", mode = "Gst",
                          nperm = 999)
cat("\nmultilocus Gst among lineages with 999 permutations: Gst =",
    round(perm$Gst[perm$Locus == "Multilocus"], 3), ", p =",
    perm$p.value[perm$Locus == "Multilocus"], "\n")

# Within the peninsular lineage, among sites
pen <- droplevels(filter(arapat, Species == "Peninsula"))
gst_pen <- genetic_structure(pen, stratum = "Population", mode = "Gst")
cat("multilocus Gst among peninsular sites:",
    round(gst_pen$Gst[gst_pen$Locus == "Multilocus"], 3), "\n")

# ---- 4. Population Graph of all 39 sites (Figure 1) ------------------------------
t_graph <- system.time(
  graph <- population_graph(arapat, stratum = "Population", decorate = TRUE)
)
# Majority lineage at each site (some sites hold beetles of two lineages).
lineage <- arapat |>
  count(Population, Species) |>
  group_by(Population) |>
  summarise(Lineage = as.character(Species[which.max(n)]),
            Mixed = n_distinct(Species) > 1)
V(graph)$Lineage <- lineage$Lineage[match(V(graph)$name,
                                           as.character(lineage$Population))]
comp <- components(graph)
cat(sprintf("\nPopulation Graph: %d nodes, %d edges, %d components (sizes %s), built in %.2f s\n",
            vcount(graph), ecount(graph), comp$no,
            paste(sort(comp$csize, decreasing = TRUE), collapse = ", "),
            t_graph[["elapsed"]]))
cat("mixed-lineage sites:",
    paste(lineage$Population[lineage$Mixed], collapse = ", "), "\n")
cat("component membership by majority lineage:\n")
print(table(Component = comp$membership, Lineage = V(graph)$Lineage))
cat("mean degree", round(mean(degree(graph)), 2), "; leaves",
    sum(degree(graph) == 1), "\n")

# Edges between lineages
el <- as_edgelist(graph)
lin_e <- cbind(V(graph)$Lineage[match(el[, 1], V(graph)$name)],
               V(graph)$Lineage[match(el[, 2], V(graph)$name)])
cat("edges joining sites of different lineages:",
    sum(lin_e[, 1] != lin_e[, 2]), "of", nrow(el), "\n")
print(table(apply(lin_e, 1, function(z) paste(sort(z), collapse = "-"))))

coords <- strata_coordinates(arapat)

# Base map: a Natural Earth land outline when rnaturalearth and
# rnaturalearthdata are installed (install.packages() both once), as in the
# Visualization chapter of the documentation; otherwise the graphs are drawn
# on plain longitude-latitude axes.  Needs gstudio >= 1.15.8 and ggplot2 >= 3.5.
base <- NULL
if (requireNamespace("rnaturalearth", quietly = TRUE) &&
    requireNamespace("rnaturalearthdata", quietly = TRUE)) {
  land <- rnaturalearth::ne_countries(
    country = c("Mexico", "United States of America"),
    scale = "medium", returnclass = "sf")
  base <- ggplot() +
    geom_sf(data = land, fill = "grey92", colour = "grey55", linewidth = 0.3) +
    coord_sf(xlim = c(-117.5, -108.5), ylim = c(22.5, 32.5), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.background = element_rect(fill = "aliceblue", colour = NA))
}
cat("\nbase map:", if (is.null(base)) "none (rnaturalearth not installed)"
    else "Natural Earth land outline", "\n")
plain_axes <- function(p) {
  if (!is.null(base)) return(p)
  p + labs(x = "Longitude", y = "Latitude") + theme_minimal(base_size = 11)
}

p1 <- plot(graph, layout = coords, node_fill = "Lineage",
           node_labels = "none", base = base) |>
  plain_axes() +
  labs(fill = "Lineage")
ggsave("media/fig1-popgraph-arapat.png", p1, width = 6, height = 6.5,
       dpi = 300, bg = "white")

# Conditional genetic distance (cGD) and its relation to geography
cgd <- genetic_distance(arapat, stratum = "Population", mode = "cGD")
geo <- strata_distance(coords, mode = "Circle")
nm  <- intersect(rownames(as.matrix(cgd)), rownames(as.matrix(geo)))
cg  <- as.matrix(cgd)[nm, nm]; gg <- as.matrix(geo)[nm, nm]
fin <- is.finite(cg[upper.tri(cg)])
cat(sprintf("\ncGD: %d site pairs, %d connected (finite), Spearman rho(cGD, km) = %.2f\n",
            sum(upper.tri(cg)), sum(fin),
            cor(cg[upper.tri(cg)][fin], gg[upper.tri(gg)][fin],
                method = "spearman")))

# ---- 5. Sampling diagnostics -----------------------------------------------------
diag_all <- popgraph_diagnostics(arapat, graph = graph)
cat("\n"); print(diag_all)
write.csv(diag_all$strata, "data/diagnostics_strata.csv", row.names = FALSE)

# ---- 6. Gravity field at both bandwidths (Figure 2) --------------------------
field <- c(`gamma = 1 (local mean)` = gravity_field(graph, gamma = 1),
           `gamma = 1/2 (degree-neutral)` = gravity_field(graph, gamma = 0.5))
cat("\n"); print(field$diagnostics, digits = 3, row.names = FALSE)
print(gravity_field(graph, gamma = 0.5))
print(gravity_field(graph, gamma = 1))
nodes <- as.data.frame(field, what = "nodes")
edges <- as.data.frame(field, what = "edges")
write.csv(nodes, "data/gravity_nodes.csv", row.names = FALSE)
write.csv(edges, "data/gravity_edges.csv", row.names = FALSE)

# Edges whose direction changes between bandwidths
ew <- edges |>
  mutate(key = paste(pmin(from, to), pmax(from, to)),
         sgn = sign(delta) * ifelse(from < to, 1, -1)) |>
  select(panel, key, sgn) |>
  tidyr::pivot_wider(names_from = panel, values_from = sgn)
flip <- sum(ew[[2]] != ew[[3]], na.rm = TRUE)
cat(sprintf("edges whose sign differs between gamma = 1 and 1/2: %d of %d\n",
            flip, nrow(ew)))
# Population-level ordering at the two bandwidths
sw <- nodes |> select(panel, Stratum, S) |>
  tidyr::pivot_wider(names_from = panel, values_from = S)
cat(sprintf("Spearman rho between S at gamma = 1 and 1/2: %.2f\n",
            cor(sw[[2]], sw[[3]], method = "spearman")))
cat(sprintf("cor(S, degree): gamma = 1 %.2f; gamma = 1/2 %.2f\n",
            with(filter(nodes, panel == unique(nodes$panel)[1]), cor(S, degree)),
            with(filter(nodes, panel == unique(nodes$panel)[2]), cor(S, degree))))

p2 <- plot(field, layout = coords, node_labels = "none", base = base) |>
  plain_axes()
ggsave("media/fig2-gravity-field-arapat.png", p2, width = 9, height = 5.5,
       dpi = 300, bg = "white")

# ---- 7. Directed output: dcGD graph ---------------------------------------------
dgraph <- dcgd(graph, gamma = 0.5)
cat(sprintf("\ndcgd(): directed = %s, %d arcs over %d edges\n",
            is_directed(dgraph), ecount(dgraph), ecount(graph)))
print(head(as_data_frame(dgraph, what = "edges"), 4))

# ---- 8. Export -------------------------------------------------------------------
gsf <- to_sf(graph)
cat("\nto_sf():", paste(class(gsf), collapse = "/"), "\n")
print(lapply(gsf, function(z) class(z)[1]))

sessionInfo()
sink()
