# R/case_study/araptus_gamma_sweep.R
#
# Bandwidth sensitivity for the Araptus case study (2026-11-23).
# Same geometric sqrt(2) grid as R/analysis/bandwidth_gamma_sweep.R,
# gamma = 0.25 ... 2, plus gamma = Inf (topological floor), applied to the
# central (South) and northern (North) components saved by
# R/case_study/case_study_arapat.R.  Per gamma: degree share, cor(S, k), cor(S, S0),
# r(S, d) with its spectral-null p, and the directional IBGD comparison.
# Also reports, descriptively, the gamma* at which cor(S, -1/k) = 0
# (the per-graph rule rejected in Supplement S8.2).

suppressMessages( { library( gstudio ); library( igraph ) } )
load( "data/case_study_arapat.rda" )
set.seed( 2013 )

grid <- c( 2^seq( -2, 1, by = 0.5 ), Inf )

sweep_one <- function( k ) {
  sub <- induced_subgraph( graph, which( V( graph )$Component == k ) )
  dk  <- d[ V( sub )$name ]
  sep <- abs( outer( dk, dk, "-" ) )
  fwd <- outer( dk, dk, function( di, dj ) dj > di )
  rows <- lapply( grid, function( g ) {
    ss <- source_sink_test( sub, x = dk, gamma = g, nperm = 999 )
    sc <- ss$scores
    fp <- tryCatch( ibgd( sub, mode = "dcGD", distance = sep, forward = fwd, gamma = g ),
                    error = function( e ) NULL )
    data.frame( gamma        = g,
                degree_share = attr( sc, "degree_share" ),
                cor_S_k      = attr( sc, "cor_S_degree" ),
                cor_S_S0     = attr( sc, "cor_S_S0" ),
                r_S_d        = unname( ss$statistic ),
                p_S_d        = ss$p.value,
                dR2          = if( is.null( fp ) ) NA else unname( fp$statistic ),
                fwd_rev      = if( is.null( fp ) ) NA else
                                 unname( fp$estimate["forward slope"] / fp$estimate["reverse slope"] ) )
  } )
  do.call( rbind, rows )
}

gamma_star <- function( k ) {
  sub <- induced_subgraph( graph, which( V( graph )$Component == k ) )
  kk  <- degree( sub )
  f   <- function( lg ) cor( source_sink_scores( sub, 2^lg )$S, -1 / kk )
  lg  <- seq( -4, 2, by = 0.05 )
  v   <- sapply( lg, f )
  i   <- which( diff( sign( v ) ) != 0 )
  if( !length( i ) ) return( c( gamma_star = NA, min_cor = min( v ), at = 2^lg[ which.min( abs( v ) ) ] ) )
  c( gamma_star = 2^uniroot( f, lg[ c( i[1], i[1] + 1 ) ] )$root, min_cor = NA, at = NA )
}

out <- lapply( expansions, sweep_one )
for( h in names( out ) ) {
  cat( sprintf( "\n==== %s (component %s)\n", h, expansions[ h ] ) )
  print( format( out[[ h ]], digits = 3 ), row.names = FALSE )
  cat( "gamma* (cor(S, -1/k) = 0): " ); print( gamma_star( expansions[ h ] ) )
}
araptus_gamma_sweep <- out
save( araptus_gamma_sweep, file = "data/derived/araptus_gamma_sweep.rda" )
