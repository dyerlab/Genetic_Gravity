#' Case study: directional gene flow in the peninsular lineage of the
#' Sonoran Desert bark beetle, Araptus attenuatus (Garrick et al. 2013).
#'
#' Hypotheses (a priori, from Garrick et al. 2013, Fig. 6): two range
#' expansions out of a single mid-peninsular refuge near 27.5 N (between
#' their local populations 9 and 10), inferred from separate regressions of
#' nuclear allelic richness on latitude: the central regional population
#' (CBP) expanded southward and the northern one (NBP) northward.
#'
#' The Population Graph is estimated for the whole peninsular lineage, as the
#' regional populations were identified from all of the data together.  Its
#' connected components recapitulate those regional populations: component 1
#' is central (CBP), component 3 northern (NBP), and the two dyads are Cape
#' Region sites.  Each expansion is tested within its own component, as
#' Garrick et al. fit a separate regression to each regional population;
#' gravity is computed within a component, so the scores are those of the
#' full graph.  Covariate for the source-sink gradient test is distance from
#' the refuge along the long axis of the peninsula; directional IBGD treats a
#' pair (i, j) as forward when j lies farther from the refuge than i.
#'
#' Output: data/case_study_arapat.rda and the numbers printed below.

rm(list=ls())
library( tidyverse )
library( gstudio )
library( igraph )

set.seed( 2013 )
refuge_lat <- 27.5

# Peninsular lineage, sites with at least four individuals.
data( arapat )
arapat |>
  filter( Species == "Peninsula" ) |>
  droplevels() |>
  group_by( Population ) |>
  filter( n() >= 4 ) |>
  ungroup() |>
  droplevels() -> peninsula

graph <- population_graph( peninsula, stratum = "Population", decorate = TRUE )
nodes <- V( graph )$name
comp  <- components( graph )
V( graph )$Component <- factor( comp$membership )


# Long axis of the peninsula: first principal axis of the site coordinates
# in km (equirectangular projection about the mean latitude), oriented so
# that position increases to the north.
lat0 <- mean( V( graph )$Latitude )
xy   <- cbind( x = V( graph )$Longitude * 111.32 * cos( lat0 * pi / 180 ),
               y = V( graph )$Latitude * 110.57 )
pc   <- prcomp( xy, center = TRUE, scale. = FALSE )
v    <- pc$rotation[, 1]
if( v["y"] < 0 ) v <- -v
a    <- setNames( as.vector( sweep( xy, 2, pc$center ) %*% v ), nodes )

# The refuge is the point on the axis at 27.5 N; d is distance from it.
a_ref <- ( refuge_lat * 110.57 - pc$center["y"] ) / v["y"]
d     <- abs( a - a_ref )
V( graph )$axis   <- a
V( graph )$refuge <- d


# One expansion per component: central (1) to the south, northern (3) to the
# north.  Every site in each component lies on its expansion's side of the
# refuge, so distance from the refuge is position along the expansion.
expansions <- c( South = "1", North = "3" )

test_expansion <- function( k ) {
  sub <- induced_subgraph( graph, which( V( graph )$Component == k ) )
  dk  <- d[ V( sub )$name ]
  sep <- abs( outer( dk, dk, "-" ) )
  fwd <- outer( dk, dk, function( di, dj ) dj > di )
  list( graph   = sub,
        field   = c( `gamma = 1/2` = gravity_field( sub, gamma = 0.5 ),
                     `gamma = 1` = gravity_field( sub, gamma = 1 ) ),
        ss_half = source_sink_test( sub, x = dk, gamma = 0.5, nperm = 999 ),
        ss_one  = source_sink_test( sub, x = dk, gamma = 1, nperm = 999 ),
        fit_cgd = ibgd( sub, mode = "cgd", distance = sep, nperm = 999 ),
        fit_pgd = ibgd( sub, mode = "pgd", distance = sep, forward = fwd, gamma = 0.5 ) )
}
results <- lapply( expansions, test_expansion )
stopifnot( all( V( results$South$graph )$Latitude < refuge_lat ),
           all( V( results$North$graph )$Latitude > refuge_lat ) )


# Summary
cat( sprintf( "%d individuals, %d sites, %.1f-%.1f N\n", nrow( peninsula ), length( nodes ),
              min( V( graph )$Latitude ), max( V( graph )$Latitude ) ) )
cat( sprintf( "%d edges; %d components of size %s\n", ecount( graph ), comp$no,
              paste( comp$csize, collapse = ", " ) ) )
cat( sprintf( "axis: PC1 explains %.1f%% of coordinate variance; refuge at a = %.0f km\n",
              100 * pc$sdev[1]^2 / sum( pc$sdev^2 ), a_ref ) )
for( h in names( results ) ) {
  r <- results[[ h ]]
  cat( sprintf( "\n==== %s expansion: component %s, %d sites, %.1f-%.1f N, %d edges, mean degree %.1f; d %.0f-%.0f km\n",
                h, expansions[ h ], vcount( r$graph ), min( V( r$graph )$Latitude ),
                max( V( r$graph )$Latitude ), ecount( r$graph ), mean( degree( r$graph ) ),
                min( d[ V( r$graph )$name ] ), max( d[ V( r$graph )$name ] ) ) )
  print( r$field$diagnostics )
  print( gravity_field( r$graph, gamma = 0.5 ) )
  print( r$ss_half )
  print( r$ss_one )
  print( r$fit_cgd )
  print( r$fit_pgd )
}

save( peninsula, graph, a, a_ref, d, expansions, results,
      file = "data/case_study_arapat.rda" )
