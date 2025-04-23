########################################################################################

# Code used to generate weighted overlay map of threats of T. matacus.

########################################################################################

# Packages needed ------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)
library(tmap)
library(tmaptools)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(ggcorrplot)

# Base maps for visual inspection -----------------------------------------
# Load South America country limits

# Country limits
Sudam_raw <- ne_countries(scale = 10, 
                          continent = "south america", 
                          returnclass = "sf")

French_guiana <- ne_countries(scale = 10,
                              type = "map_units",
                              geounit = "French Guiana",
                              returnclass = "sf")
Sudam <- st_union(Sudam_raw, French_guiana)

# Load Chaco region by Morrone 2022.
Chaco_biome <- st_read("C:/Capas_GIS/Bioregiones/BioRegionesNeotrop_Morrone2022/Geographic_vector_all/NeotropicMap_Geo.shp") |> 
  filter(Provincias == "Chaco province") |> 
  st_make_valid()

qtm(Chaco_biome)

# Crop Sudam using Proj Area for further plotting
Sudam_Proj <- st_intersection(Sudam,Chaco_biome)

# Base map
bb <- c(-72, # longitud mínima 
        -37, # latitud mínima
        -52, # longitud máxima
        -14) # latitud máxima
bb_sf <- st_sfc(st_polygon(list(matrix(c(bb[1], bb[2], 
                                         bb[1], bb[4], 
                                         bb[3], bb[4], 
                                         bb[3], bb[2], 
                                         bb[1], bb[2]), 
                                       ncol = 2, byrow = TRUE))))

mapa_base <- tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgray") +
  tmap_options(check.and.fix = TRUE)
mapa_base



# 1. Habitat loss footprint ---------------------------------------------------------

# Load habitat change raster

habitat_change <- rast("02_Habitat_loss/habitat_change.tif")
qtm(habitat_change)

# Reclassify values: Convert -1 (habitat loss) and stable unsuitable habitat (0)
# to 1, and all others (1, 2) to 0
reclass_matrix <- matrix(c(-1, 1,   # Habitat loss → 1 (threat)
                           0, 1,   # Stable unsuitable → 0
                           1, 0,   # Habitat gained → 0
                           2, 0),  # Stable suitable → 0
                         ncol = 2, byrow = TRUE)
lost_areas <- classify(habitat_change, reclass_matrix)
qtm(lost_areas)

# Apply focal mean filter (10 km buffer ~57x57 pixel window)
focal_window <- matrix(1, nrow = 57, ncol = 57)  # Define window
smoothed_loss <- focal(lost_areas, w = focal_window, fun = mean, na.policy = "all", na.rm = TRUE) |> 
  crop(Chaco_biome) |> 
  mask(Chaco_biome)
qtm(smoothed_loss)

# Map with smoothes habitat loss
mapa_smoothed_HL <- tm_shape(Chaco_biome) + tm_polygons(col = "lightgray") +
  tm_shape(smoothed_loss) + tm_raster(style = "fixed",
                                             breaks = seq(0, 1, by = 0.05),
                                             palette = get_brewer_pal("YlOrRd", n = 20,
                                                                      plot = F)) +
  tm_shape(Sudam_Proj) + tm_borders() +
  tm_shape(Chaco_biome) + tm_borders(col = "black", lwd = 2) +
  tm_layout(title = "",
            legend.outside = T,
            frame = F)
mapa_smoothed_HL

# Save image
tmap_save(mapa_smoothed_HL,
          "plots/smoothed_HL.pdf",
          dpi = 300)

# Weighted overlay with hunting pressure -----------------------------------

# Load hunting pressure raster
hunting_risk_raw <- rast("data/HP_tolymat.tif")

# Ensure same resolution and extent of loss
hunting_risk <- resample(hunting_risk_raw, smoothed_loss, method = "bilinear") |> 
  crop(Chaco_biome) |> 
  mask(Chaco_biome)
qtm(hunting_risk)

# Define weights
w_habitat_loss <- 0.6
w_hunting <- 0.4

# Weighted sum of threats
threat_overlay <- (smoothed_loss * w_habitat_loss) + (hunting_risk * w_hunting)
qtm(threat_overlay)

# Save output
writeRaster(threat_overlay, "03_Overlap/threat_overlay_weighted.tif", overwrite=TRUE)

threat_overlay <- rast("03_Overlap/threat_overlay_weighted.tif")

# Count the total number of non-NA pixels
total_non_na_pixels <- sum(!is.na(values(threat_overlay)))

# Count the number of pixels with values above 0.5
pixels_above_0.5 <- sum(values(threat_overlay) > 0.5, na.rm = TRUE)

# Calculate the proportion
proportion_above_0.5 <- pixels_above_0.5 / total_non_na_pixels

# Print the result
print(proportion_above_0.5)


# Map of hunting
mapa_hunting <- tm_shape(Chaco_biome) + tm_polygons(col = "lightgray") +
  tm_shape(hunting_risk) + tm_raster(style = "fixed",
                                      breaks = seq(0, 1, by = 0.05),
                                      palette = get_brewer_pal("YlOrRd", n = 20,
                                                               plot = F)) +
  tm_shape(Sudam_Proj) + tm_borders(col = "darkgrey", lty = "dashed") +
  tm_shape(Chaco_biome) + tm_borders(col = "black", lwd = 2) +
  tm_layout(title = "",
            legend.outside = T,
            frame = F)
mapa_hunting

# Save image
tmap_save(mapa_hunting,
          "plots/hunting_footprint.pdf",
          dpi = 300)

# Map of overlayed threats
mapa_overlap <- tm_shape(Chaco_biome) + tm_polygons(col = "lightgray") +
  tm_shape(threat_overlay) + tm_raster(style = "cont",
                                     #breaks = seq(0, 1, by = 0.05),
                                     palette = get_brewer_pal("YlOrRd", n = 20,
                                                              plot = F)) +
  tm_shape(Sudam_Proj) + tm_borders() +
  tm_shape(Chaco_biome) + tm_borders(col = "black", lwd = 2) +
  tm_layout(title = "",
            legend.outside = T,
            frame = F)
mapa_overlap

# Save image
tmap_save(mapa_overlap,
          "plots/overlap.pdf",
          dpi = 300)

