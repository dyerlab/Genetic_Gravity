#' Testing gravity fields on sonoran desert bark beetles

# Load in the araptus site data
library( tidyverse ) 
library( gstudio ) 
library(igraph)

data( arapat ) 
arapat |> 
  population_graph( decorate=TRUE ) -> graph
comp <- components(graph)
V( graph )$Subgraph <- factor( comp$membership )

graph |>
  gravity_field() -> gf.arapat
plot( gf.arapat, layout="kk" )

# let's find the graph components
gf.arapat$nodes |>
  filter( Subgraph == 2 ) |>
  ggplot( aes(Latitude, S, color=Subgraph) ) + 
  stat_smooth(se=FALSE, method="lm", formula = "y~x") + 
  geom_point() 


