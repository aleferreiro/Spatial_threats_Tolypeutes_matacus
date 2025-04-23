########################################################################################

# Code used to model the habitat suitability of T. matacus for 1985.

# This model was then used to estimate habitat loss when compared with 2015 model.

########################################################################################

# Packages needed -----------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(tmap)
library(tmaptools)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(ggcorrplot)
library(flexsdm)

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


# 1. Ocurrence data -----------------------------------------------------------------------

# From Ferreiro et al. 2022 (doi: 10.1007/s10914-022-09627-3)
Occs_Ferreiro <- read.csv("data/Occs_Ferreiro_etal_2022.csv",
                          sep = ";" )
Occs_Ferreiro_reduced <- Occs_Ferreiro |>  
  dplyr::filter(Temporalidad == "Actual") |> 
  dplyr::filter(Ano_colecta < 2015 & Ano_colecta > 1980) |>  # Filter rows to retain samples from 70s and 00s
  dplyr::filter(!Origen == "Bordigon, 2006") # to remove an outlier presence

# Convert to sf
Occs_Ferreiro_reduced_sf <- Occs_Ferreiro_reduced |> 
  st_as_sf(coords = c("Longitud", "Latitud"),
           crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# Plot map with presences
map_Ferreiro <- mapa_base +
  tm_shape(Occs_Ferreiro_reduced_sf) + tm_symbols(col = "yellow", size = 0.5, shape = 21, alpha = 0.5) +
  tm_shape(Chaco_biome) + tm_borders(col = "black", lwd = 2)
  
map_Ferreiro

# Reduce to species, lat, long data to model
Occs_Toly_mat <- Occs_Ferreiro_reduced |>
  select(c("Longitud", "Latitud")) |> 
  mutate(sp = c(rep("Tolypeutes_matacus", nrow(Occs_Ferreiro_reduced))), .before = 1)

# Save file
write.csv(Occs_Toly_mat,
          file = "data/Occs_ENM_habitat_loss.csv",
          row.names = F, quote = F)



# 2. Accesible area -----------------------------------------------------------------------------

# As detailed land cover data was only available for the Chaco biome, we used it extent to calibrate models.

# 3. Predictors -----------------------------------------------------------
## 3.1. Soil variables -----------------------------------------------------
# Layers used for calibration of the models
# We used four physical soil variables to assess habitat suitability.
# Load layers cropped and masked from the Chaco_biome
path_Soil <- list.files(path = "C:/Capas_GIS/SoilGrids",
                        pattern = "*.tif$",full.names = TRUE)
Soil <- rast(path_Soil) 

# Reproject Chaco_biome to match raster CRS
Chaco_biome_rasterCRS <- st_transform(Chaco_biome, crs = crs(Soil))

# Convert to SpatVector (required for cropping/masking)
Chaco_biome_vect <- vect(Chaco_biome_rasterCRS)

# Crop and mask the raster
Soil_cropped <- Soil |> 
  crop(Chaco_biome_rasterCRS) |> 
  mask(Chaco_biome_rasterCRS)

# Reproject cropped raster to WGS 84 (EPSG:4326)
Soil_cropped_WGS <- project(Soil_cropped, "EPSG:4326", method = "bilinear")
names(Soil_cropped_WGS) <- c("cfvo", "clay", "sand", "silt")

# Save soil predictors
lapply(names(Soil_cropped_WGS), 
       function(x){
         writeRaster(Soil_cropped_WGS[[x]], 
                     paste0("02_Habitat_loss/01_1985_model/predictors/",x,".tif"),
                     overwrite=TRUE)})

sand <- rast("02_Habitat_loss/01_1985_model/predictors/sand.tif")

qtm(sand)

## 3.2. Land cover/land use data -------------------------------------------

# We downloaded a the lulc layers from the mapbiomas Chaco project, 
# and use it to estimate three land cover layers:

# Load LULC raster (original classification)
lulc <- rast("C:/Capas_GIS/Land_cover/MapBiomas/Chaco/mapbiomas-chaco-collection-50-none-1985.tif")


### 3.2.1. Woodland --------------------------------------------------
# Define woodland class (adjust class ID as needed)
woodland_classes <- c(3, 4, 45, 6)  
# Reclassify LULC to binary (1 = closed woodland, 0 = others)
woodland_binary <- classify(Soil_cropped_WGS, rcl = cbind(woodland_classes, 1), others = 0)
# Resample to match soil raster resolution and extent
woodland_resampled <- resample(woodland_binary, Soil_cropped_WGS, method = "bilinear")
# Aggregate LULC to calculate percentage per soil grid cell
# It Compute percentage of closed woodland within 3 km buffer (within a 25x25 pixel window)
woodland <- focal(woodland_resampled, 
                  w = matrix(1, 25, 25), fun = mean) |> 
  crop(Chaco_biome) |> 
  mask(Chaco_biome) 

qtm(woodland)
names(woodland) <- "woodland"

writeRaster(woodland,
            "02_Habitat_loss/01_1985_model/predictors/woodland.tif",
            overwrite = T)

# Reload layer
woodland <- rast("02_Habitat_loss/01_1985_model/predictors/woodland.tif")

### 3.2.2. Grasslands --------------------------------------------------
# Define grassland classes (adjust class IDs as needed)
grassland_classes <- c(12, 42, 43, 44, 11)  # General, open, closed, sparse, flooded grassland

# Reclassify LULC to binary (1 = grassland, 0 = others)
grassland_binary <- classify(Soil_cropped_WGS, rcl = cbind(grassland_classes, 1), others = 0)

# Resample to match soil raster resolution and extent
grassland_resampled <- resample(grassland_binary, Soil_cropped_WGS, method = "bilinear")

# Compute percentage of grasslands within 1 km buffer (9x9 window)
grassland <- focal(grassland_resampled, 
                   w = matrix(1, 25, 25), fun = mean) |> 
  crop(Chaco_biome) |> 
  mask(Chaco_biome) 
names(grassland) <- "grassland"

# Plot result
qtm(grassland)

# Save output raster
writeRaster(grassland,
            "02_Habitat_loss/01_1985_model/predictors/grasslands.tif",
            overwrite = TRUE)

grassland <- rast("02_Habitat_loss/01_1985_model/predictors/grasslands.tif")
       
# % of grasslands (all types)
# % of pastures
# % of cropland
# % of urban areas

### 3.2.3. Pastures --------------------------------------------------
# Define pasture class (adjust class ID as needed)
pasture_class <- c(15)  # Pasture category

# Reclassify LULC to binary (1 = pasture, 0 = others)
pasture_binary <- classify(Soil_cropped_WGS, rcl = cbind(pasture_class, 1), others = 0)

# Resample to match soil raster resolution and extent
pasture_resampled <- resample(pasture_binary, Soil_cropped_WGS, method = "bilinear")

# Compute percentage of pastures within 1 km buffer (9x9 window)
pasture <- focal(pasture_resampled, 
                 w = matrix(1, 25, 25), fun = mean) |> 
  crop(Chaco_biome) |> 
  mask(Chaco_biome) 

names(pasture) <- "pastures"

# Plot result
qtm(pasture)

# Save output raster
writeRaster(pasture,
            "02_Habitat_loss/01_1985_model/predictors/pastures.tif",
            overwrite = TRUE)

pasture <- rast("02_Habitat_loss/01_1985_model/predictors/pastures.tif")

### 3.2.4. Croplands --------------------------------------------------
# Define cropland classes (adjust class IDs as needed)
cropland_classes <- c(9, 36, 57, 58)  # Single crop (57) and multiple crop (58)

# Reclassify LULC to binary (1 = cropland, 0 = others)
cropland_binary <- classify(Soil_cropped_WGS, rcl = cbind(cropland_classes, 1), others = 0)

# Resample to match soil raster resolution and extent
cropland_resampled <- resample(cropland_binary, Soil_cropped_WGS, method = "bilinear")

# Compute percentage of cropland within 1 km buffer (9x9 window)
cropland <- focal(cropland_resampled, 
                  w = matrix(1, 25, 25), fun = mean) |> 
  crop(Chaco_biome) |> 
  mask(Chaco_biome) 

names(cropland) <- "cropland"

# Plot result
qtm(cropland)

# Save output raster
writeRaster(cropland,
            "02_Habitat_loss/01_1985_model/predictors/cropland.tif",
            overwrite = TRUE)

cropland <- rast("02_Habitat_loss/01_1985_model/predictors/cropland.tif")

## 3.3. Predictors stack --------------------------------------------------

# We combine predictor within one raster
Predictors <- c(woodland, grassland, pasture, cropland, sand)

# 4. Check correlation --------------------------------------------------

# Predictors correlation 
Predictors_Corr_result <- correct_colinvar(env_layer = Predictors,
                                           method = c('pearson', th='0.7'),
                                           maxcell = 50000)
Cor_table <- Predictors_Corr_result$cor_table # extract correlation matrix
Corplot <- ggcorrplot(Cor_table,
                      hc.order = F,
                      lab = T,
                      lab_size = 6) +             # plot correlation heatmap
  theme(axis.text.x = element_text(size = 20),  # Cambia el tamaño del texto en el eje X
        axis.text.y = element_text(size = 20))  # Cambia el tamaño del texto en el eje Y

Corplot 
ggsave("plots/S3a_Corplot.pdf", Corplot)    # save the plot
ggsave("plots/S3a_Corplot.png", Corplot) 

# Based on the correlation we decided to keep 7 predictors to model:


Predictors_uncorr <- Predictors |> 
  subset(c("woodland", "grassland", "pastures", "cropland", "sand"))

# 5. Thinning occurrences -----------------------------------------------

# We used flexsdm package to filter presence data retaining records separated
# by at least 10km
Occs_raw <- Occs_Toly_mat |> 
  select(2:3) |> 
  mutate(id = 1:length(Occs_Toly_mat$sp), .before = 1)

set.seed(10)
Occs_thin <- occfilt_geo(
  data = Occs_raw,
  x = "Longitud",
  y = "Latitud",
  method = c('defined', d = 10),
  env_layer = Predictors_uncorr,
  prj = crs(Predictors_uncorr))

# Convert to sf
Occs_thin_sf <- Occs_thin |> 
  st_as_sf(coords = c("Longitud", "Latitud"),
           crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# Map occs thinned
map_Ferreiro_thin <- map_Ferreiro +
  tm_shape(Occs_thin_sf) + tm_symbols(col = "red", size = 0.5, shape = 21, alpha = 0.5)
map_Ferreiro_thin

# 6. Data partition ----------------------------------------------------------
# To evaluate the models, I need independent data that will allow me
# to assess how well they perform. Since I don't have independent sources,
# I will randomly partition the data so that I calibrate the models with one part
# and evaluate them with the other.

# First, I need to add a column to identify the presences
Occs_pres <- Occs_thin |> 
  mutate(pr_ab = rep("1", length(Occs_thin$id)))

# Run the function to select the best partition
set.seed(10)

occ_part <- part_random(
  data = Occs_pres,
  pr_ab = "pr_ab",
  method = c(method = "rep_kfold", folds = 4, replicates = 5)
)

# Extract environmental values from ocurrences
Occs_pa <- occ_part %>%
  sdm_extract(
    data = .,
    x = "Longitud",
    y = "Latitud",
    env_layer = Predictors_uncorr,
    filter_na = TRUE
  )

Occs_pa$pr_ab <- as.numeric(Occs_pa$pr_ab)

# 7. Backgroung points sampling -------------------------------------------

# Sample background points throughout study area with random method, 
# allocating x10 background POINTS
set.seed(10)
bg <- sample_background(
    data = Occs_thin,
    x = "Longitud",
    y = "Latitud",
    n = 10000,
    #n = nrow(Occs_thin) * 10,
    method = "random",
    rlayer = Predictors_uncorr)

# Convert to sf to plot
bg_sf <- st_as_sf(bg, coords = c("Longitud","Latitud"), crs = 4326)

# Plot to see bg distributions
mapa5_bg <- mapa_base +
  tm_shape(bg_sf) + tm_symbols(size = 0.2) +
  tm_shape(Sudam) + tm_borders() +
  tm_layout(legend.show = F)
mapa5_bg

bg_envs <- bg %>%
  sdm_extract(
    data = .,
    x = "Longitud",
    y = "Latitud",
    env_layer = Predictors_uncorr,
    filter_na = TRUE
  )

# Assign background points to one of the 4 partition groups at random
set.seed(10)  # Set seed for reproducibility
for (i in 1:5) {
  bg_envs[[paste0(".part", i)]] <- sample(1:4, nrow(bg_envs), replace = TRUE)
}


# 8. Model ---------------------------------------------------------------

# Tunned maxent 
t_max <- tune_max(
  data = Occs_pa,
  response = "pr_ab",
  predictors = names(Predictors_uncorr),
  background = bg_envs,
  partition = ".part",
  grid = expand.grid(
    regmult = c(0.1, 0.25, 0.5, 1, 2, 5, 10),
    classes = c("l", "q", "lq", "lp", "lqp")
  ),
  thr = c('sensitivity', sens='0.9') ,
  metric = "TSS",
  clamp = TRUE,
  pred_type = "cloglog",
  n_cores = 2
)

# Evaluation metrics of the best model
best_model_metrics <- t_max$performance
best_model_metrics

# Table with metrics of all models analyzed
model_metrics <- t_max$hyper_performance
model_metrics_paper <- model_metrics |> 
  select(c("regmult", "classes", "TSS_mean", "TSS_sd", "AUC_mean", "AUC_sd", "OR_mean", "OR_sd"))

library(gt)

model_metrics_paper_table <- model_metrics_paper |> 
  gt() |> 
  fmt_number(columns = everything(),
             decimals = 3) |> 
  fmt_number(columns = 1,
             decimals = 2)

model_metrics_paper_table

gtsave(model_metrics_paper_table,
       "plots/HL_model_metrics.docx")



# 9. Projection ----------------------------------------------------------


## 9.1. Project to 1985 ---------------------------------------------------

# Models were projected using the same variables as used for training.

Projection_layers_1985 <- Predictors_uncorr


SDM_Current <- sdm_predict(
  models = t_max,
  pred = Projection_layers_1985,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL,
  clamp = T
)

SDM_Current_bin <- SDM_Current$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))
qtm(SDM_Current$max$sensitivity)

writeRaster(SDM_Current$max$sensitivity,
            "02_Habitat_loss/sdm_1985.tif")

sdm_1985 <- rast("02_Habitat_loss/sdm_1985.tif")
sdm_1985_binary <- sdm_1985 |> 
  classify(rcl = matrix(c(-10, 0, NA,
                          0, 1, 1), 
                        ncol = 3, byrow = TRUE))
qtm(sdm_1985_binary)
# Map with threshold
mapa_Current_bin <- tm_shape(Chaco_biome) + tm_polygons(col = "lightgray") +
  tm_shape(sdm_1985_binary) + tm_raster(style = "cat",
                                        palette = "darkgreen",
                                        labels = "Suitable") +
  tm_shape(Sudam_Proj) + tm_borders() +
  tm_shape(Chaco_biome) + tm_borders(col = "black", lwd = 2) +
  tm_layout(title = "",
            legend.outside = T,
            frame = F)
mapa_Current_bin

tmap_save(mapa_Current_bin,
          "plots/1985_HabitatModel_binary.pdf",
          dpi = 300)


## 9.2. Project to 2015 ---------------------------------------------------------

woodland_2015 <- rast("02_Habitat_loss/02_2015_model/predictors/woodland.tif")
grassland_2015 <- rast("02_Habitat_loss/02_2015_model/predictors/grasslands.tif")
cropland_2015 <- rast("02_Habitat_loss/02_2015_model/predictors/cropland.tif")
pastures_2015 <- rast("02_Habitat_loss/02_2015_model/predictors/pastures.tif")

Projection_layers_2015 <- c(woodland_2015, grassland_2015, 
                            cropland_2015, pastures_2015, sand)


SDM_Current_2015 <- sdm_predict(
  models = t_max,
  pred = Projection_layers_2015,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL,
  clamp = T
)

SDM_Current_bin_2015 <- SDM_Current_2015$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

writeRaster(SDM_Current_2015$max$sensitivity,
            "02_Habitat_loss/sdm_2015.tif")
sdm_2015 <- rast("02_Habitat_loss/sdm_2015.tif")
sdm_2015_binary <- sdm_1985 |> 
  classify(rcl = matrix(c(-10, 0, NA,
                          0, 1, 1), 
                        ncol = 3, byrow = TRUE))


# Map with threshold
mapa_Current_bin_2015 <- tm_shape(Chaco_biome) + tm_polygons(col = "lightgray") +
  tm_shape(sdm_2015_binary) +tm_raster(style = "cat",
                                            palette = "darkgreen",
                                            labels = "Suitable") +
  tm_shape(Sudam_Proj) + tm_borders() +
  tm_shape(Chaco_biome) + tm_borders(col = "black", lwd = 2) +
  tm_layout(title = "",
            legend.outside = T,
            frame = F)
mapa_Current_bin_2015


tmap_save(mapa_Current_bin_2015,
          "plots/2015_HabitatModel_binary.pdf",
          dpi = 300)


# 10. Habitat change map --------------------------------------------------

SDM_binary_1985 <-  ifel(sdm_1985 > 0, 1, 0)
qtm(SDM_binary_1985)

SDM_binary_2015 <- ifel(sdm_2015 > 0, 1, 0)
qtm(SDM_binary_2015)

# Sum both rasters
change_map <- SDM_binary_1985 + SDM_binary_2015
qtm(change_map)

habitat_loss <- (change_map == 1) & (SDM_binary_1985 == 1)  # Suitable in 1985, lost in 2015
habitat_gain <- (change_map == 1) & (SDM_binary_2015 == 1)  # Suitable in 2015, gained from 1985

# Assign categories:
final_map <- change_map
final_map[habitat_loss] <- -1  # Habitat loss coded as -1
final_map[habitat_gain] <- 1   # Habitat gain coded as 1

writeRaster(final_map,
            "02_Habitat_loss/habitat_change.tif")

final_map2 <- final_map |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))
qtm(final_map2)

# Define the color palette for the three categories
change_colors <- c("Habitat Loss" = "#E41A1C",    # Red for Habitat Loss
                   "Habitat Gain" = "blue",    # Blue for Habitat Gain
                   "Stable Suitable" = "chartreuse3") # Dark Green for Stable Suitable

# Assign labels for the categories after reclassification
category_labels <- c("Habitat Loss", "Habitat Gain", "Stable Suitable")

# Plot the habitat change map after reclassification
mapa_habitat_change <- tm_shape(Chaco_biome) + 
  tm_polygons(col = "lightgray") + 
  
  tm_shape(final_map2) + 
  tm_raster(style = "cat", 
            palette = change_colors, 
            labels = category_labels, 
            title = "Habitat Change (1985-2015)") +
  
  tm_shape(Sudam_Proj) + tm_borders() +
  tm_shape(Chaco_biome) + tm_borders(col = "black", lwd = 2) +
  
  tm_layout(title = "",
            legend.outside = TRUE,
            frame = FALSE)

# Show the updated map
mapa_habitat_change

# Save the figure
tmap_save(mapa_habitat_change,
          "plots/HabitatChange.pdf",
          dpi = 300)


