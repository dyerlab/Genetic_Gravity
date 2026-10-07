rm(list=ls())
library( tidyverse )
library( gstudio )
library( igraph )


read_population( "data/PedimaFinal.csv", 
                    type="zyme", 
                    locus.columns = 10:15, 
                    sep="\t" ) |>
    select( -SoC, -IoL, -MSP, -Geneland, -ID) |>
    mutate( Population = factor( Pop) ) |>
    select( Population, 
            Latitude = Lat, Longitude = Lon, 
            everything() ) -> pedima


summary( pedima )

pedima |>
    population_graph() -> ped.graph 

plot( ped.graph )

gravity_field( ped.graph ) -> gf.pedima 
plot( gf.pedima )



