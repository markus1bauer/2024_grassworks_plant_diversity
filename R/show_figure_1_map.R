#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Show figure map ####
# Map of the study sites
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Markus Bauer
# 2024-08-07



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(sf)

### Start ###
rm(list = ls())

### Functions ###
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_line(colour = NA),
    text = element_text(size = 10, color = "black"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(.5, 0, 0, 0, "cm")
  )
}

### Load data ###
sites <- st_read(
  here("data", "raw", "data_processed_sites_epsg4326.gpkg")
  )
germany <- geodata::gadm(
  country = "DEU", level = 0, version = "latest", resolution = 2, path = here()
) %>%
  st_as_sf() %>%
  st_set_crs(4326)
ecoregions <- st_read(here("data", "raw", "ecoregions2017.shp")) %>%
  st_transform(crs = 4326) %>%
  filter(
    ECO_NAME == "Western European broadleaf forests" |
      ECO_NAME == "Central European mixed forests" |
      ECO_NAME == "European Atlantic mixed forests" |
      ECO_NAME == "Baltic mixed forests" |
      ECO_NAME == "Alps conifer and mixed forests"
  )
elevation <- elevatr::get_elev_raster(
  locations = data.frame(x = c(5.5, 15.5), y = c(47, 55.5)),
  prj = "+proj=longlat +datum=WGS84 +no_defs",
  clip = "bbox", z = 3
) %>%
  terra::rast()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot #######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



clip <- st_intersection(ecoregions, germany)
grad <- tidyterra::hypso.colors2(10, "dem_poster")



## 1 Map with ecoregions #######################################################


graph_sites <- ggplot() +
  geom_sf(
    data = clip, colour = "black", fill = "transparent", size = 1,
    linetype = "dashed"
  ) +
  geom_sf(data = germany, colour = "black", fill = "transparent", size = 1) +
  geom_sf(data = sites, colour = "black", size = 1) +
  annotate(
    geom = "label", x = c(9.9, 12.85, 11), y = c(52.8, 51.55, 48.7),
    label = c("North", "Centre", "South"), fill = "transparent", size = 3.5
  ) +
  annotate(
    geom = "text", x = c(8, 13, 13.8, 10), y = c(53.1, 53.85, 52.2, 50),
    lineheight = .8, size = 3.5,
    label = c(
      "European Atlantic\nmixed forests",
      "Baltic\nmixed forests",
      "Central European\nmixed forests",
      "Western European\nbroadleaf forests"
      )
  ) +
  ggspatial::annotation_north_arrow(
    which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(),
    height = unit(2, "cm"),
    width = unit(2, "cm"),
    pad_y = unit(.5, "cm"),
    pad_x = unit(0.1, "cm")
  ) +
  ggspatial::annotation_scale(
    pad_y = unit(.7, "cm"),
    pad_x = unit(2, "cm")
  ) +
  theme_mb();graph_sites

### Save -----------------------------------------------------------------------

ggsave(
  "figure_1_map_ecoregions_300dpi_12x15cm.tiff",
  dpi = 300, width = 12, height = 15, units = "cm",
  path = here("outputs", "figures")
)



## 2 Map with elevation ########################################################


graph_sites <- ggplot() +
  tidyterra::geom_spatraster(data = elevation) +
  geom_sf(
    data = clip, colour = "black", fill = "transparent", size = 1,
    linetype = "dashed"
    ) +
  geom_sf(data = germany, colour = "black", fill = "transparent", size = 1) +
  geom_sf(data = sites, colour = "black", size = 1) +
  annotate(
    geom = "label", x = c(9.9, 12.85, 11), y = c(52.8, 51.55, 48.7),
    label = c("North", "Centre", "South"), fill = "transparent", size = 3.5,
  ) +
  annotate(
    geom = "text", x = c(8, 13, 13.8, 10), y = c(53.1, 53.85, 52.2, 50),
    lineheight = .8, size = 3.5,
    label = c(
      "European Atlantic\nmixed forests",
      "Baltic\nmixed forests",
      "Central European\nmixed forests",
      "Western European\nbroadleaf forests"
    )
  ) +
  ggspatial::annotation_north_arrow(
    which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(),
    height = unit(2, "cm"),
    width = unit(2, "cm"),
    pad_y = unit(.5, "cm"),
    pad_x = unit(0.1, "cm")
  ) +
  ggspatial::annotation_scale(
    pad_y = unit(.7, "cm"),
    pad_x = unit(2, "cm")
  ) +
  scale_fill_gradientn(colours = grad, na.value = NA, name = "Elevation [m]") +
  theme_mb();graph_sites


### Save -----------------------------------------------------------------------

ggsave(
  "figure_1_map_elevation_300dpi_15x15cm.tiff",
  dpi = 300, width = 15, height = 15, units = "cm",
  path = here("outputs", "figures")
)
