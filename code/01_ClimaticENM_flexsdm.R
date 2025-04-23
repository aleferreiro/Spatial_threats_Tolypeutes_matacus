########################################################################################

# Code used to model the ENM of T. matacus, using flexsdm package and assess current 
# climatic suitable areas and their change under future scenarios of global warming.

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
  filter(Provincias == "Chaco province")
qtm(Chaco_biome)

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
  tm_shape(Occs_Ferreiro_reduced_sf) + tm_symbols(col = "yellow", size = 0.5, shape = 21, alpha = 0.5)
map_Ferreiro

# Reduce to species, lat, long data to model
Occs_Toly_mat <- Occs_Ferreiro_reduced |>
  select(c("Longitud", "Latitud")) |> 
  mutate(sp = c(rep("Tolypeutes_matacus", 144)), .before = 1)

# Save file
write.csv(Occs_Toly_mat,
          file = "data/Occs_ENM.csv",
          row.names = F, quote = F)



# 2. Accesible area -----------------------------------------------------------------------------
# We obtained the calibration/background area using the flexsdm package

#Previous to obtain the calibration area we should load the Current climatic data
path_CurrentBios <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Current_2.1/CHELSA_Bioclim",
                               pattern = "*.tif$",full.names = TRUE)
CurrentBios <- rast(path_CurrentBios) |> 
  subset(c(1:9,12:17)) |> 
  terra::crop(Sudam) |> 
  terra::mask(Sudam)

names(CurrentBios) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                        "Bio_14", "Bio_15", "Bio_16", "Bio_17",
                        "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                        "Bio_06", "Bio_07")

# Create a calibration area using a buffer of 150km to the minimum convex polygon of the presences
Area_M_Toly_mat <-  calib_area(
  data = Occs_Toly_mat,
  x = 'Longitud',
  y = 'Latitud',
  method =  c('bmcp', width = 150000),
  crs = crs(CurrentBios)
) # create a calibration area with 100 km buffer around MCP

# Convert to sf
Area_M_Toly_mat_sf <- st_as_sf(Area_M_Toly_mat)  

# Plot Area M
map_Area_M <- map_Ferreiro +
  tm_shape(Area_M_Toly_mat_sf) + tm_borders(col = "green", 
                                            lwd = 2) 
map_Area_M

# Save the polygon vector
st_write(Area_M_Toly_mat_sf,
         "data/Area_M.gpkg")
Area_M_Toly_mat_sf <- st_read("data/Area_M.gpkg")
Area_M_Toly_mat <- vect(Area_M_Toly_mat_sf)

# 3. Predictors --------------------------------------------------------------

## 3.1. Bioclimatic variables ----------------------------------------------

# Recorto Capas Bioclimaticas por Area M.
Predictor_bios <- CurrentBios |> 
  terra::crop(Area_M_Toly_mat) |> 
  terra::mask(Area_M_Toly_mat)

names(Predictor_bios) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                           "Bio_14", "Bio_15", "Bio_16", "Bio_17",
                           "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                           "Bio_06", "Bio_07")
Predictor_bios_2 <- subset(Predictor_bios,
                           c("Bio_01", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                             "Bio_06", "Bio_07", "Bio_10", "Bio_11", "Bio_12", 
                             "Bio_13", "Bio_14", "Bio_15", "Bio_16", "Bio_17"))

## 3.2. Predictors selection -----------------------------------------------

# First, we estimate correlations among the 15 Bioclimatic variables, to 
# generate groups of correlated prdictors. Of those groups, we chose one 
# predictor according to its relevance for the species (Table S1)

# Predictors correlation 
Predictors_Corr_result <- correct_colinvar(env_layer = Predictor_bios_2,
                                           method = c('pearson', th='0.7'),
                                           maxcell = 50000)
Cor_table <- Predictors_Corr_result$cor_table # extract correlation matrix
Corplot <- ggcorrplot(Cor_table,
                      hc.order = T,
                      lab = T,
                      lab_size = 3)              # plot correlation heatmap
Corplot 
ggsave("plots/S1a_Corplot_Climatic.pdf", Corplot)    # save the plot
ggsave("plots/S1a_Corplot_Climatic.png", Corplot) 

# Based on correlation clusters and Biological meaningful we select 4 variables
# See Table S1
Predictor_Bios_Uncorr <- Predictor_bios_2 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14"))

# 4. Projection layers ----------------------------------------------------

# As area of projection we used a higher buffer of 400 km to evaluate how suitability
# may change in potential new habitats for the species.
Area_Proj_Toly_mat <-  calib_area(
  data = Occs_Toly_mat,
  x = 'Longitud',
  y = 'Latitud',
  method =  c('bmcp', width = 400000),
  crs = crs(CurrentBios)
) # create a projection area with 400 km buffer around MCP

# Convert to sf
Area_Proj_Toly_mat_sf <- st_as_sf(Area_Proj_Toly_mat)  

# Plot Area Projection area
map_Area_Proj <- map_Ferreiro +
  tm_shape(Area_Proj_Toly_mat_sf) + tm_borders(col = "green",
                                               lwd = 2) 
map_Area_Proj

# Crop Sudam using Proj Area for further plotting
Sudam_Proj <- st_intersection(Sudam, Area_Proj_Toly_mat_sf)
qtm(Sudam_Proj)

# We used projection layers from CHELSA v. 2.1 for:
# Future scenarios: SSP 1-2.6, SSP 3-7.0, and SSP 5-8.5
# GCMs: GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL

# Current -----------------------------------------------------------------

CurrentBios_proj <- CurrentBios |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)


## GDFL_126 ----------------------------------------------------------------
path_GDFL_126 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/GFDL-ESM4_ssp126/",
                            pattern = "*.tif$",full.names = TRUE)
GDFL_126 <- rast(path_GDFL_126)
names(GDFL_126) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                     "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                     "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                     "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_GDFL_126 <- GDFL_126 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## GDFL_370 ----------------------------------------------------------------
path_GDFL_370 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/GFDL-ESM4_ssp370/",
                            pattern = "*.tif$",full.names = TRUE)
GDFL_370 <- rast(path_GDFL_370)
names(GDFL_370) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                     "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                     "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                     "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_GDFL_370 <- GDFL_370 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)
# plot(Proj_GDFL_370)

## GDFL_585 ----------------------------------------------------------------
path_GDFL_585 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/GFDL-ESM4_ssp585/",
                            pattern = "*.tif$",full.names = TRUE)
GDFL_585 <- rast(path_GDFL_585)
names(GDFL_585) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                     "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                     "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                     "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_GDFL_585 <- GDFL_585 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)
# plot(Proj_GDFL_585)


## IPSL_126 ----------------------------------------------------------------
path_IPSL_126 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/IPSL-CM6A-LR_ssp126/",
                            pattern = "*.tif$",full.names = TRUE)
IPSL_126 <- rast(path_IPSL_126)
names(IPSL_126) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                     "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                     "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                     "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_IPSL_126 <- IPSL_126 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## IPSL_370 ----------------------------------------------------------------
path_IPSL_370 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/IPSL-CM6A-LR_ssp370/",
                            pattern = "*.tif$",full.names = TRUE)
IPSL_370 <- rast(path_IPSL_370)
names(IPSL_370) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                     "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                     "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                     "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_IPSL_370 <- IPSL_370 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## IPSL_585 ----------------------------------------------------------------
path_IPSL_585 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/IPSL-CM6A-LR_ssp585/",
                            pattern = "*.tif$",full.names = TRUE)
IPSL_585 <- rast(path_IPSL_585)
names(IPSL_585) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                     "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                     "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                     "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_IPSL_585 <- IPSL_585 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## MPI_126 ----------------------------------------------------------------
path_MPI_126 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/MPI-ESM1-2-HR_ssp126/",
                           pattern = "*.tif$",full.names = TRUE)
MPI_126 <- rast(path_MPI_126)
names(MPI_126) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                    "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                    "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                    "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_MPI_126 <- MPI_126 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## MPI_370 ----------------------------------------------------------------
path_MPI_370 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/MPI-ESM1-2-HR_ssp370/",
                           pattern = "*.tif$",full.names = TRUE)
MPI_370 <- rast(path_MPI_370)
names(MPI_370) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                    "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                    "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                    "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_MPI_370 <- MPI_370 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## MPI_585 ----------------------------------------------------------------
path_MPI_585 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/MPI-ESM1-2-HR_ssp585/",
                           pattern = "*.tif$",full.names = TRUE)
MPI_585 <- rast(path_MPI_585)
names(MPI_585) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                    "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                    "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                    "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_MPI_585 <- MPI_585 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## MRI_126 ----------------------------------------------------------------
path_MRI_126 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/MRI-ESM2-0_ssp126/",
                           pattern = "*.tif$",full.names = TRUE)
MRI_126 <- rast(path_MRI_126)
names(MRI_126) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                    "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                    "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                    "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_MRI_126 <- MRI_126 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## MRI_370 ----------------------------------------------------------------
path_MRI_370 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/MRI-ESM2-0_ssp370/",
                           pattern = "*.tif$",full.names = TRUE)
MRI_370 <- rast(path_MRI_370)
names(MRI_370) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                    "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                    "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                    "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_MRI_370 <- MRI_370 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## MRI_585 ----------------------------------------------------------------
path_MRI_585 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/MRI-ESM2-0_ssp585/",
                           pattern = "*.tif$",full.names = TRUE)
MRI_585 <- rast(path_MRI_585)
names(MRI_585) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                    "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                    "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                    "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_MRI_585 <- MRI_585 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## UKESM1_126 ----------------------------------------------------------------
path_UKESM1_126 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/UKESM1-0-LL_ssp126/",
                              pattern = "*.tif$",full.names = TRUE)
UKESM1_126 <- rast(path_UKESM1_126)
names(UKESM1_126) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                       "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                       "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                       "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_UKESM1_126 <- UKESM1_126 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## UKESM1_370 ----------------------------------------------------------------
path_UKESM1_370 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/UKESM1-0-LL_ssp370/",
                              pattern = "*.tif$",full.names = TRUE)
UKESM1_370 <- rast(path_UKESM1_370)
names(UKESM1_370) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                       "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                       "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                       "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_UKESM1_370 <- UKESM1_370 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

## UKESM1_585 ----------------------------------------------------------------
path_UKESM1_585 <- list.files(path = "C:/Capas_GIS/ClimVars/CHELSA/Future_2.1/2071_2100/UKESM1-0-LL_ssp585/",
                              pattern = "*.tif$",full.names = TRUE)
UKESM1_585 <- rast(path_UKESM1_585)
names(UKESM1_585) <- c("Bio_01", "Bio_10", "Bio_11", "Bio_12", "Bio_13", 
                       "Bio_14", "Bio_15", "Bio_16", "Bio_17", "Bio_18",
                       "Bio_19", "Bio_02", "Bio_03", "Bio_04", "Bio_05", 
                       "Bio_06", "Bio_07", "Bio_08", "Bio_09")
Proj_UKESM1_585 <- UKESM1_585 |> 
  subset(c("Bio_01", "Bio_04", "Bio_12", "Bio_14")) |> 
  terra::crop(Area_Proj_Toly_mat) |> 
  terra::mask(Area_Proj_Toly_mat)

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
  env_layer = Predictor_Bios_Uncorr,
  prj = crs(Predictor_Bios_Uncorr))

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
# I will partition the data so that I calibrate the models with one part
# and evaluate them with the other. I will use the function `part_sblock()` 
# from the `flexsdm` package to partition the data into spatial blocks.

# First, I need to add a column to identify the presences
Occs_pres <- Occs_thin |> 
  mutate(pr_ab = rep("1", length(Occs_thin$id)))

# Run the function to select the best partition
set.seed(10)

occ_part <- part_sblock(
  data = Occs_pres,
  env_layer = Predictor_Bios_Uncorr,
  pr_ab = "pr_ab",
  x = "Longitud",
  y = "Latitud",
  n_part = 4,
  min_res_mult = 3,
  max_res_mult = 200,
  num_grids = 30,
  min_occ = 10,
  prop = 0.5
)
Occs_p <- occ_part$part

Occs_p_sf <- st_as_sf(Occs_p,
                      coords = c("x","y"),
                      crs = 4326)


# Transform best block partition to a raster layer with same resolution and extent than
# predictor variables
block_layer <- get_block(env_layer = Predictor_Bios_Uncorr, best_grid = occ_part$grid)

# Plot best spatial blocking partitions
# Ploteo para ver como quedó la partición en el espacio
mapa4 <- mapa_base +
  tm_shape(block_layer) + tm_raster(palette = get_brewer_pal("Accent", n = 4, plot = F)) +
  tm_shape(Occs_p_sf) + tm_symbols(size = 0.2, col = ".part", palette = get_brewer_pal("Accent", n = 4, plot = F)) +
  tm_shape(Sudam) + tm_borders() +
  tm_shape(Area_M_Toly_mat_sf) + tm_borders(col = "red", 
                                            lwd = 2) +
  tm_layout(legend.show = F)
mapa4 

tmap_save(tm = mapa4,
          filename = "plots/data_partition.pdf")

tmap_save(tm = mapa4,
          filename = "plots/data_partition.png")


# Number of presences per block
Occs_p %>%
  dplyr::group_by(.part) %>%
  dplyr::count()

# Extract Bioclimatic variables from ocurrences
Occs_pa <- Occs_p %>%
  sdm_extract(
    data = .,
    x = "x",
    y = "y",
    env_layer = Predictor_Bios_Uncorr,
    filter_na = TRUE
  )

Occs_pa$pr_ab <- as.numeric(Occs_pa$pr_ab)
# 7. Backgroung points sampling -------------------------------------------

# Spatial blocks where species occurs
# Sample background points throughout study area with random method, 
# allocating 10X the number of presences a background
set.seed(10)
bg <- lapply(1:4, function(x) {
  sample_background(
    data = Occs_p,
    x = "x",
    y = "y",
    n = 10000,
    method = "random",
    rlayer = block_layer,
    maskval = x,
    calibarea = Area_M_Toly_mat
  )
}) |> 
  bind_rows()
bg <- sdm_extract(data = bg, x = "x", y = "y", env_layer = block_layer)

# Convert to sf to plot
bg_sf <- st_as_sf(bg, coords = c("x","y"), crs = 4326)

# Plot to see bg distributions
mapa5_bg <- mapa_base +
  tm_shape(block_layer) + tm_raster(palette = get_brewer_pal("Accent", n = 4, plot = F)) +
  tm_shape(bg_sf) + tm_symbols(size = 0.2, col = ".part", palette = get_brewer_pal("Accent", n = 4, plot = F)) +
  tm_shape(Sudam) + tm_borders() +
  tm_layout(legend.show = F)
mapa5_bg

bg_envs <- bg %>%
  sdm_extract(
    data = .,
    x = "x",
    y = "y",
    env_layer = Predictor_Bios_Uncorr,
    filter_na = TRUE
  )


# 8. Models ---------------------------------------------------------------
t_max <- tune_max(
  data = Occs_pa,
  response = "pr_ab",
  predictors = names(Predictor_Bios_Uncorr),
  background = bg_envs,
  partition = ".part",
  grid = expand.grid(
    regmult = c(0.1, 0.25, 0.5, 1, 2, 5, 10),
    classes = c("l", "q", "lq", "lqp")
  ),
  thr = c('sensitivity', sens='0.9') ,
  metric = "TSS",
  clamp = TRUE,
  pred_type = "cloglog",
  n_cores = 2
)

# Evaluation metrics of the best model
best_model_metrics <- t_max$performance

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
       "plots/model_metrics.docx")

# 9. Geographic projection --------------------------------------------------
# 9.1. Current --------------------------------------------------------------
SDM_Current <- sdm_predict(
  models = t_max,
  pred = CurrentBios_proj,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL,
  clamp = T
)

SDM_Current_bin <- SDM_Current$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

writeRaster(SDM_Current$max$sensitivity,
            "01_ClimaticENM/ENM_Current_suitability.tif")

# Map with threshold
mapa_Current_bin <- tm_shape(Area_Proj_Toly_mat_sf) + tm_polygons(col = "lightgray") +
  tm_shape(SDM_Current_bin) + tm_raster(style = "fixed",
                                        breaks = seq(0, 1, by = 0.05),
                                        palette = viridisLite::viridis(20)
                                        ) +
  tm_shape(Sudam_Proj) + tm_borders() +
  tm_shape(Chaco_biome) + tm_borders(col = "black", lwd = 2) +
  tm_layout(title = "Current GCM mean",
            legend.outside = T,
            frame = F)
mapa_Current_bin

tmap_save(mapa_Current_bin,
          "plots/flexsdm/SDM_Current_mean.pdf",
          dpi = 300)

# Extent of suitable area
Suitable_area <- expanse(x = SDM_Current_bin)

# 9.2. SSP_126 --------------------------------------------------------------
SDM_GDFL_126 <- sdm_predict(
  models = t_max,
  pred = Proj_GDFL_126,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL,
  clamp = T
)
SDM_GDFL_126_bin <- SDM_GDFL_126$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_IPSL_126 <- sdm_predict(
  models = t_max,
  pred = Proj_IPSL_126,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_IPSL_126_bin <- SDM_IPSL_126$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_MPI_126 <- sdm_predict(
  models = t_max,
  pred = Proj_MPI_126,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_MPI_126_bin <- SDM_MPI_126$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_MRI_126 <- sdm_predict(
  models = t_max,
  pred = Proj_MRI_126,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_MRI_126_bin <- SDM_MRI_126$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_UKESM1_126 <- sdm_predict(
  models = t_max,
  pred = Proj_UKESM1_126,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_UKESM1_126_bin <- SDM_UKESM1_126$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_126_mean <- terra::mean(c(SDM_GDFL_126_bin, SDM_IPSL_126_bin,
                              SDM_MPI_126_bin, SDM_MRI_126_bin,
                              SDM_UKESM1_126_bin),
                            na.rm = T) 

SDM_126_sd <- terra::stdev(c(SDM_GDFL_126_bin, SDM_IPSL_126_bin,
                              SDM_MPI_126_bin, SDM_MRI_126_bin,
                              SDM_UKESM1_126_bin),
                           na.rm = T)

writeRaster(SDM_126_mean,
            "01_ClimaticENM/ENM_126_suitability.tif")
writeRaster(SDM_126_sd,
            "01_ClimaticENM/ENM_126_suitability_sd.tif")


# Map with threshold
# Mean
mapa_SDM_126_mean <- tm_shape(Area_Proj_Toly_mat_sf) + tm_polygons(col = "lightgray") +
  tm_shape(SDM_126_mean) + tm_raster(style = "fixed",
                                     breaks = seq(0, 1, by = 0.05),
                                     palette = viridisLite::viridis(20)
  ) +
  tm_shape(Sudam_Proj) + tm_borders() +
  tm_shape(Chaco_biome) + tm_borders(col = "black", lwd = 2) +
  tm_layout(title = "SSP 1-2.6 GCM mean",
            legend.show = F,
            frame = F)
mapa_SDM_126_mean

tmap_save(mapa_SDM_126_mean,
          "plots/flexsdm/mapa_SDM_126_mean.pdf",
          dpi = 300)

# SD
mapa_SDM_126_sd <- mapa_base +
  tm_shape(SDM_126_sd) + tm_raster(style = "cont",
                                   palette = get_brewer_pal("Reds", n = 7, plot = F)
  ) +
  tm_shape(Sudam) + tm_borders() +
  tm_layout(title = "SSP 1-2.6 GCM sd",
            legend.show = T)
mapa_SDM_126_sd

tmap_save(mapa_SDM_126_sd,
          "plots/flexsdm/mapa_SDM_126_sd.pdf",
          dpi = 300)



# 9.3. SSP_370 --------------------------------------------------------------
SDM_GDFL_370 <- sdm_predict(
  models = t_max,
  pred = Proj_GDFL_370,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL,
  clamp = T
)
SDM_GDFL_370_bin <- SDM_GDFL_370$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_IPSL_370 <- sdm_predict(
  models = t_max,
  pred = Proj_IPSL_370,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_IPSL_370_bin <- SDM_IPSL_370$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_MPI_370 <- sdm_predict(
  models = t_max,
  pred = Proj_MPI_370,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_MPI_370_bin <- SDM_MPI_370$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_MRI_370 <- sdm_predict(
  models = t_max,
  pred = Proj_MRI_370,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_MRI_370_bin <- SDM_MRI_370$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_UKESM1_370 <- sdm_predict(
  models = t_max,
  pred = Proj_UKESM1_370,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_UKESM1_370_bin <- SDM_UKESM1_370$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_370_mean <- terra::mean(c(SDM_GDFL_370_bin, SDM_IPSL_370_bin,
                              SDM_MPI_370_bin, SDM_MRI_370_bin,
                              SDM_UKESM1_370_bin),
                            na.rm = T) 

SDM_370_sd <- terra::stdev(c(SDM_GDFL_370_bin, SDM_IPSL_370_bin,
                             SDM_MPI_370_bin, SDM_MRI_370_bin,
                             SDM_UKESM1_370_bin),
                           na.rm = T)

writeRaster(SDM_370_mean,
            "01_ClimaticENM/ENM_370_suitability.tif")
writeRaster(SDM_370_sd,
            "01_ClimaticENM/ENM_370_suitability_sd.tif")

# Map with threshold
# Mean
mapa_SDM_370_mean <- tm_shape(Area_Proj_Toly_mat_sf) + tm_polygons(col = "lightgray") +
  tm_shape(SDM_370_mean) + tm_raster(style = "fixed",
                                     breaks = seq(0, 1, by = 0.05),
                                     palette = viridisLite::viridis(20)
  ) +
  tm_shape(Sudam_Proj) + tm_borders() +
  tm_shape(Chaco_biome) + tm_borders(col = "black", lwd = 2) +
  tm_layout(title = "SSP 3-7.0 GCM mean",
            legend.show = F)
mapa_SDM_370_mean

tmap_save(mapa_SDM_370_mean,
          "plots/flexsdm/mapa_SDM_370_mean.pdf",
          dpi = 300)

# SD
mapa_SDM_370_sd <- mapa_base +
  tm_shape(SDM_370_sd) + tm_raster(style = "cont",
                                   palette = get_brewer_pal("Reds", n = 7, plot = F)
  ) +
  tm_shape(Sudam) + tm_borders() +
  tm_layout(title = "SSP 3-7.0 GCM sd",
            legend.show = T)
mapa_SDM_370_sd

tmap_save(mapa_SDM_370_sd,
          "plots/flexsdm/mapa_SDM_370_sd.pdf",
          dpi = 300)
# 9.4. SSP_585 --------------------------------------------------------------
SDM_GDFL_585 <- sdm_predict(
  models = t_max,
  pred = Proj_GDFL_585,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL,
  clamp = T
)
SDM_GDFL_585_bin <- SDM_GDFL_585$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_IPSL_585 <- sdm_predict(
  models = t_max,
  pred = Proj_IPSL_585,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_IPSL_585_bin <- SDM_IPSL_585$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_MPI_585 <- sdm_predict(
  models = t_max,
  pred = Proj_MPI_585,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_MPI_585_bin <- SDM_MPI_585$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_MRI_585 <- sdm_predict(
  models = t_max,
  pred = Proj_MRI_585,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_MRI_585_bin <- SDM_MRI_585$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_UKESM1_585 <- sdm_predict(
  models = t_max,
  pred = Proj_UKESM1_585,
  thr = c('sensitivity', sens='0.9'),
  con_thr = TRUE,
  predict_area = NULL
)
SDM_UKESM1_585_bin <- SDM_UKESM1_585$max[[2]] |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = TRUE))

SDM_585_mean <- terra::mean(c(SDM_GDFL_585_bin, SDM_IPSL_585_bin,
                              SDM_MPI_585_bin, SDM_MRI_585_bin,
                              SDM_UKESM1_585_bin),
                            na.rm = T) 

SDM_585_sd <- terra::stdev(c(SDM_GDFL_585_bin, SDM_IPSL_585_bin,
                             SDM_MPI_585_bin, SDM_MRI_585_bin,
                             SDM_UKESM1_585_bin),
                           na.rm = T)
writeRaster(SDM_585_mean,
            "01_ClimaticENM/ENM_585_suitability.tif")
writeRaster(SDM_585_sd,
            "01_ClimaticENM/ENM_585_suitability_sd.tif")


# Map with threshold
# Mean
mapa_SDM_585_mean <- tm_shape(Area_Proj_Toly_mat_sf) + tm_polygons(col = "lightgray") +
  tm_shape(SDM_585_mean) + tm_raster(style = "fixed",
                                     breaks = seq(0, 1, by = 0.05),
                                     palette = viridisLite::viridis(20)
  ) +
  tm_shape(Sudam_Proj) + tm_borders() +
  tm_shape(Chaco_biome) + tm_borders(col = "black", lwd = 2) +
  tm_layout(title = "SSP 5-8.5 GCM mean",
            legend.show = F, frame = F)
mapa_SDM_585_mean

tmap_save(mapa_SDM_585_mean,
          "plots/flexsdm/mapa_SDM_585_mean.pdf",
          dpi = 300)

# SD
mapa_SDM_585_sd <- mapa_base +
  tm_shape(SDM_585_sd) + tm_raster(style = "cont",
                                   palette = get_brewer_pal("Reds", n = 7, plot = F)
  ) +
  tm_shape(Sudam) + tm_borders() +
  tm_layout(title = "SSP 5-8.5 GCM sd",
            legend.show = T)
mapa_SDM_585_sd

tmap_save(mapa_SDM_585_sd,
          "plots/flexsdm/mapa_SDM_585_sd.pdf",
          dpi = 300)

# 10. Inset map --------------------------------------------------
bb_Sudam <- c(-81, # longitud mínima 
              -55, # latitud mínima
              -33, # longitud máxima
              12) # latitud máxima


Inset_map_ENS <- tm_shape(Sudam, bbox = bb_Sudam) + tm_borders(col = "darkgrey") +
  tm_shape(Area_Proj_Toly_mat_sf) + tm_polygons(col = "lightgrey") +
  tm_shape(Sudam) + tm_borders(col = "darkgrey") +
  tm_shape(Chaco_biome) + tm_borders(lwd = 1, col = "black")
Inset_map_ENS

Inset_map_SD <- tm_shape(Sudam, bbox = bb_Sudam) + tm_borders(col = "darkgrey") +
  tm_shape(bb_sf) + tm_borders(col = "black") +
  tm_shape(Sudam) + tm_borders(col = "darkgrey") 
Inset_map_SD



tmap_save(Inset_map_ENS, "plots/flexsdm/Inset_map_ENS.pdf")

tmap_save(Inset_map_SD, "plots/flexsdm/Inset_map_SD.pdf")

