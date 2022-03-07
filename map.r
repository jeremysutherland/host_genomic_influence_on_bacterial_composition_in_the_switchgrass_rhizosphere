library(tidyverse)
library(ggmap)
library(sf)
library(mapview)
library(leaflet)
library(wesanderson)
  

mapviewOptions(basemaps = "Esri.WorldShadedRelief",
               raster.palette = grey.colors,
              vector.palette = colorRampPalette(colors = c("#000000","#C93312","#046C9A","#D69C4E","#ABDDDE")),
              na.color = "magenta",
              layers.control.pos = "topright")

locations_df <- read.csv("df.csv", header= TRUE, row.names = 1)
locations_df

locations <- as_tibble(locations_df)
Populations <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)
m = mapview(Populations, zcol = "Source", legend = TRUE)
m


## create .html and .png
mapshot(m, url = paste0(getwd(), "/map.html"))
mapshot(m, file = paste0(getwd(), "/Geo_map.pdf"),
        remove_controls = c("homeButton", "layersControl", "zoomControl"))

