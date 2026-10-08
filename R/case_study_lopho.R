#' Case study: directional gene flow in the senita cactus, Lophocereus
#' schottii (Nason et al. 2002; Dyer & Nason 2004).
#'
#' Graph: the Population Graph of Dyer & Nason (2004), as distributed with
#' gstudio >= 1.15.0 (the 50 edges of their Table 1), with site coordinates
#' from the gstudio 'baja' metadata.
#'
#' Hypothesis (a priori, from Dyer & Nason 2004): range expansion out of
#' peninsular Baja California into SenBas (southern Arizona) by long-distance
#' dispersal, then southward from SenBas through the remaining continental
#' (Sonoran) populations.  Each site is given a stage along this route:
#' peninsula (1, unordered within), SenBas (2), then the other continental
#' sites in order of decreasing latitude (3, 4, ...).
#'
#' Source-sink gradient test: S against route stage; the hypothesis predicts
#' r(S, stage) < 0, sources early on the route.  Directional IBGD: a pair
#' (i, j) is forward when j lies at a later stage than i; pairs of peninsular
#' sites carry no direction and are left out of the directional fit.
#' Separation is great-circle distance.
#'
#' Output: data/case_study_lopho.rda and the numbers printed below.

rm(list=ls())
library( tidyverse )
library( gstudio )
library( igraph )

stopifnot( packageVersion( "gstudio" ) >= "1.15.0" )
set.seed( 2004 )

data( lopho )
data( baja )
graph <- lopho
stopifnot( ecount( graph ) == 50 )          # the Dyer & Nason (2004) topology
nodes <- V( graph )$name
coords <- baja[ match( nodes, baja$Population ), ]
stopifnot( all( coords$Population == nodes ) )
V( graph )$Latitude  <- coords$Latitude
V( graph )$Longitude <- coords$Longitude
V( graph )$Region    <- factor( V( graph )$region, levels = c( "Baja", "Sonora" ),
                                labels = c( "Peninsula", "Continent" ) )
continent <- V( graph )$Region == "Continent"

# Stage along the hypothesized route.
lat   <- setNames( V( graph )$Latitude, nodes )
stage <- setNames( rep( 1, length( nodes ) ), nodes )
stage[ "SenBas" ] <- 2
southward <- nodes[ continent & nodes != "SenBas" ]
stage[ southward ] <- 2 + rank( -lat[ southward ] )
V( graph )$stage <- stage
forward <- outer( stage, stage, function( si, sj ) sj > si )

sep <- data.frame( Stratum = nodes, Longitude = V( graph )$Longitude, Latitude = lat ) |>
  strata_distance( mode = "Circle" ) |>
  as.matrix()
sep <- sep[ nodes, nodes ]


# Gravity at the degree-neutral bandwidth and at the local mean.
field <- c( `gamma = 1/2` = gravity_field( graph, gamma = 0.5 ),
            `gamma = 1` = gravity_field( graph, gamma = 1 ) )
field$diagnostics

# Source-sink gradient test, S against route stage.
ss_half <- source_sink_test( graph, x = stage, gamma = 0.5, nperm = 999 )
ss_one  <- source_sink_test( graph, x = stage, gamma = 1, nperm = 999 )

# Isolation by graph distance: symmetric (cGD) vs directional (pGD).
fit_cgd <- ibgd( graph, mode = "cgd", distance = sep, nperm = 999 )
fit_pgd <- ibgd( graph, mode = "pgd", distance = sep, forward = forward, gamma = 0.5 )


# Summary
cat( sprintf( "%d sites (%d peninsular, %d continental), %.1f-%.1f N; %d edges, %d component(s), mean degree %.1f\n",
              length( nodes ), sum( !continent ), sum( continent ), min( lat ), max( lat ),
              ecount( graph ), components( graph )$no, mean( degree( graph ) ) ) )
cross <- apply( as_edgelist( graph ), 1, function( e ) continent[ match( e[1], nodes ) ] != continent[ match( e[2], nodes ) ] )
cat( "Edges between regions:", apply( as_edgelist( graph )[ cross, ], 1, paste, collapse = "-" ), "\n" )
cat( "Bridges:", apply( as_edgelist( graph )[ bridges( graph ), , drop = FALSE ], 1, paste, collapse = "-" ), "\n" )
ss_half
ss_one
fit_cgd
fit_pgd

S <- field$nodes |> filter( panel == "gamma = 1/2" ) |> select( Stratum, degree, S )
S$Region   <- V( graph )$Region[ match( S$Stratum, nodes ) ]
S$Latitude <- lat[ S$Stratum ]
S$stage    <- stage[ S$Stratum ]
print( arrange( S, desc( S ) ) )
print( S |> group_by( Region ) |> summarize( mean_S = mean( S ), n = n() ) )
with( subset( S, Region == "Continent" & Stratum != "SenBas" ),
      cat( sprintf( "Continent below SenBas: rho(S, latitude) = %.2f\n", cor( S, Latitude, method = "spearman" ) ) ) )
with( subset( S, Region == "Peninsula" ),
      cat( sprintf( "Peninsula: rho(S, latitude) = %.2f\n", cor( S, Latitude, method = "spearman" ) ) ) )

# Edges at SenBas, the putative long-distance colonist.
field$edges |> filter( panel == "gamma = 1/2", from == "SenBas" | to == "SenBas" ) |> print()

# Directional fit on the continental pairs alone (SenBas, then southward).
fw_cont <- forward & outer( continent, continent, "&" )
fit_pgd_cont <- ibgd( graph, mode = "pgd", distance = sep, forward = fw_cont, gamma = 0.5 )
fit_pgd_cont

save( graph, stage, sep, forward, field, ss_half, ss_one, fit_cgd, fit_pgd, fit_pgd_cont,
      file = "data/case_study_lopho.rda" )
