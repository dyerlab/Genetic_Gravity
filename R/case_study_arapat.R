#' Case study: directional gene flow in the peninsular lineage of the
#' Sonoran Desert bark beetle, Araptus attenuatus (Garrick et al. 2013).
#'
#' Hypothesis (a priori, from Garrick et al. 2013): expansion in two
#' directions from a single mid-peninsular refuge near 27.5 N.  Covariate for
#' the source-sink gradient test is distance from the refuge along the long
#' axis of the peninsula; directional IBGD treats a pair (i, j) as forward
#' when j lies farther from the refuge than i.
#'
#' Output: data/case_study_arapat.rda and the numbers printed below.

rm(list=ls())
library( tidyverse )
library( gstudio )
library( igraph )

set.seed( 2013 )
refuge_lat <- 27.5

# Peninsular lineage only, sites with at least four individuals.
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

# The refuge is the point on the axis at 27.5 N.
a_ref <- ( refuge_lat * 110.57 - pc$center["y"] ) / v["y"]
d     <- abs( a - a_ref )
V( graph )$axis   <- a
V( graph )$refuge <- d

sep     <- abs( outer( a, a, "-" ) )
forward <- outer( d, d, function( di, dj ) dj > di )


# Gravity at the degree-neutral bandwidth and at the local mean.
field <- c( `gamma = 1/2` = gravity_field( graph, gamma = 0.5 ),
            `gamma = 1` = gravity_field( graph, gamma = 1 ) )
field$diagnostics

# Source-sink gradient test, S against distance from the refuge.
ss_half <- source_sink_test( graph, x = d, gamma = 0.5, nperm = 999 )
ss_one  <- source_sink_test( graph, x = d, gamma = 1, nperm = 999 )

# Isolation by graph distance: symmetric (cGD) vs directional (pGD).
fit_cgd <- ibgd( graph, mode = "cgd", distance = sep, nperm = 999 )
fit_pgd <- ibgd( graph, mode = "pgd", distance = sep, forward = forward, gamma = 0.5 )


# Summary
cat( sprintf( "%d individuals, %d sites, %.1f-%.1f N\n", nrow( peninsula ), length( nodes ),
              min( V( graph )$Latitude ), max( V( graph )$Latitude ) ) )
cat( sprintf( "%d edges; %d components of size %s\n", ecount( graph ), comp$no,
              paste( sort( comp$csize, decreasing = TRUE ), collapse = ", " ) ) )
cat( sprintf( "axis: PC1 explains %.1f%% of coordinate variance; refuge at a = %.0f km; d range %.0f-%.0f km\n",
              100 * pc$sdev[1]^2 / sum( pc$sdev^2 ), a_ref, min( d ), max( d ) ) )
ss_half
ss_one
fit_cgd
fit_pgd

save( peninsula, graph, a, a_ref, d, field, ss_half, ss_one, fit_cgd, fit_pgd,
      file = "data/case_study_arapat.rda" )
