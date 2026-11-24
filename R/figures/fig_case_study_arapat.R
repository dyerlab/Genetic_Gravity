# Fig_CaseStudyArapat: gravity field of the peninsular Araptus attenuatus
# Population Graph at the degree-neutral bandwidth (gamma = 1/2), drawn in
# graph space (Kamada-Kawai layout) and over the Baja California peninsula,
# with the hypothesized refuge (27.5 N; Garrick et al. 2013) marked on the
# map.  Nodes are coloured by connected component.  Reads the graph saved by
# R/case_study/case_study_arapat.R.
# v2: node labels placed with ggrepel (gstudio's fixed labels overlapped in the
# dense central cluster); full-opacity surface on the map panel.
# v3: the two two-site components (S = 0 and Delta = 0 by symmetry) are left
# out of the map panel but kept in the graph panel.
# v4: CBP only, its graph estimated from CBP sites alone (superseded).
# v5: back to the whole peninsular graph (components 1 and 3 carry the two
# expansions); refuge latitude drawn on the map.
# v6: v5 with the two-site components left out of the map panel, as in v3,
# and nodes coloured by component again (current gstudio); no component legend.
# v7: Natural Earth 1:10m coastline (rnaturalearthhires), surface masked to each component's hull and
# clipped to land (plot.gravity_field(mask = "hull", clip_to_base = TRUE)).
# Three versions of the map panel: with the graph (-graph), the graph faded
# (-faded: nodes, edges and arrows faded, labels at full strength) and the
# surface alone (-surface).

library( tidyverse )
library( gstudio )
library( igraph )
library( patchwork )
library( sf )

FIG <- c( graph   = "media/fig-case-study-arapat-v7-graph.png",
          faded   = "media/fig-case-study-arapat-v7-faded.png",
          surface = "media/fig-case-study-arapat-v7-surface.png" )
# Opacity of each part of the map panel's graph; labels are drawn here (ggrepel).
graph_alpha <- list( graph   = c( nodes = 1,   edges = 1,    arrows = 1,    labels = 1 ),
                     faded   = c( nodes = 0.5, edges = 0.35, arrows = 0.35, labels = 1 ),
                     surface = c( nodes = 0,   edges = 0,    arrows = 0,    labels = 0 ) )
refuge_lat <- 27.5

load( "data/case_study_arapat.rda" )
field <- gravity_field( graph, gamma = 0.5 )

# Map panel: drop the two-site components, which carry no asymmetry.  S is
# computed within each component, so the remaining scores are unchanged.
keep  <- V( graph )$name[ components( graph )$csize[ components( graph )$membership ] > 2 ]
graph_map <- induced_subgraph( graph, keep )
V( graph_map )$Component <- droplevels( V( graph_map )$Component )
field_map <- gravity_field( graph_map, gamma = 0.5 )

# Same component colours in both panels.
component_colours <- setNames( scales::viridis_pal()( nlevels( V( graph )$Component ) ),
                               levels( V( graph )$Component ) )

# Natural Earth 1:10m countries (rnaturalearthhires), cropped to the region.
sf_use_s2( FALSE )
box <- c( xmin = -116, xmax = -109, ymin = 23.5, ymax = 31 )
land <- rnaturalearth::ne_countries( scale = 10, country = c( "Mexico", "United States of America" ),
                                     returnclass = "sf" ) |>
  st_geometry() |>
  st_crop( box ) |>
  suppressWarnings()

ggplot() +
  geom_sf( data = land, fill = "grey85", colour = "grey60", linewidth = 0.15 ) +
  xlab( "Longitude" ) +
  ylab( "Latitude" ) +
  theme_bw( base_size = 12 ) +
  coord_sf( xlim = c(-115, -110.3), ylim = c(24.4, 29.8), expand = FALSE ) -> baja

# gstudio draws labels as fixed text above each node, which collides in the
# dense central cluster; repel them instead, at the plot's node positions.
repel_labels <- function( p, alpha = 1 ) {
  nd <- plot_nodes( p )
  p + ggrepel::geom_text_repel( data = nd, aes( x, y, label = node ), size = 2.7,
                                fontface = "bold", colour = "grey10", point.size = 4,
                                box.padding = 0.4, point.padding = 0.15,
                                min.segment.length = 0.3, segment.colour = "grey30",
                                segment.size = 0.3, max.overlaps = Inf, seed = 1,
                                bg.colour = "white", bg.r = 0.12, alpha = alpha, inherit.aes = FALSE )
}

p_graph <- plot( field, layout = "kk", node_fill = "Component", node_labels = "none",
                 subtitle_stats = FALSE ) |>
  repel_labels() +
  scale_colour_manual( values = component_colours ) +
  labs( title = "A" ) +
  theme( legend.position = "none" )
map_panel <- function( a ) {
  p_map <- plot( field_map, base = baja, node_fill = "Component", node_labels = "none",
                 subtitle_stats = FALSE, mask = "hull", clip_to_base = TRUE,
                 graph_alpha = a[ c( "nodes", "edges", "arrows" ) ] )
  if( a[[ "labels" ]] > 0 )
    p_map <- repel_labels( p_map, alpha = a[[ "labels" ]] )
  if( a[[ "nodes" ]] > 0 ) {
    p_map <- p_map +
      scale_colour_manual( values = component_colours, guide = "none" )
    # Keep node sizes on the graph panel's scale (the map has fewer sites).
    size_scale <- p_map$scales$get_scales( "size" )$clone()
    size_scale$limits <- range( V( graph )$size )
    p_map <- p_map + size_scale
  }
  p_map +
    geom_hline( yintercept = refuge_lat, linetype = "dashed", colour = "grey20" ) +
    annotate( "text", x = -110.4, y = refuge_lat + 0.12, label = "refuge, 27.5°N",
              hjust = 1, size = 3.2, colour = "grey20" ) +
    labs( title = "B" )
}

for( v in names( FIG ) ) {
  p <- p_graph + map_panel( graph_alpha[[ v ]] ) + plot_layout( widths = c(1, 1) )
  ggsave( FIG[[ v ]], p, width = 13, height = 6.5, dpi = 200, bg = "white" )
  cat( "Wrote", FIG[[ v ]], "\n" )
}
