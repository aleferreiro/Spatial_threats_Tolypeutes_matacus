#################################################################################

# Code used to generate weighted overlay map of threats of T. matacus.

########################################################################################

# Packages and functions needed ------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)
library(tmap)
library(tmaptools)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(patchwork)

# SRI index
sri_index <- function(x) {
  Prod_RiSi <- Sudam_PA_rast*x # Primero multiplicar el Si raster con el raster de PAs.
  Sum_RixSi <- global(Prod_RiSi, fun = "sum", na.rm=T) # Sumatoria de los valores de celdas de Producto Si x Ri
  Sum_Si <- global(x, fun = "sum", na.rm=T) # Sumatoria de valores de pixeles del raster Si
  SRI <-  (Sum_RixSi[1,1]/Sum_Si[1,1])*100
  SRI                       # la última línea es lo que devolverá la función
}

# PSA index
area_idonea_protegida <- function(x) {
  Prod_RiSi <- Sudam_PA_rast*x # Primero multiplicar el Si raster con el raster de PAs.
  AIP <-  expanse(Prod_RiSi, unit = "km", transform = T)
  AIP                       # la última línea es lo que devolverá la función
} 

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





## 0.1. Load PAs data ------------------------------------------------------------


# Download Protected Areas (PA) data
# Create a folder to store the layers
# dir.create("WDPA")
# 
# # Download the file: Go to the WDPA website, copy the download link instead of clicking "Download", and paste it below
# download.file(url = "https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_Apr2021_Public_shp.zip",
#               destfile = "WDPA/wdpa_Apr2021_shp.zip")
# 
# # Unzip the downloaded file
# unzip(zipfile = "WDPA/wdpa_Apr2021_shp.zip", exdir = "WDPA")
# 
# # The database includes both polygons and points. Points represent areas without proper boundary data,
# # so they will not be considered in this analysis.
# 
# # The global shapefile is split into three parts, each compressed.
# # Unzip each of the three polygon shapefiles
# dir.create("WDPA/wdpa_Apr2021_Public_shp_0")
# unzip(zipfile = "WDPA/wdpa_Apr2021_Public_shp_0.zip", exdir = "WDPA/wdpa_Apr2021_Public_shp_0")

# Load the shapefile
# Sudam_pa_data <-  st_read(dsn = "WDPA//wdpa_Apr2021_Public_shp_0/WDPA_Apr2021_Public_shp-polygons.shp") |> 
#   select(c("NAME", "IUCN_CAT", "MARINE", "GIS_AREA", "ISO3" )) |> # Reduce dataset to be easiest to handle
#   filter(MARINE == 0) |>   # Remove marine PA
#   st_intersection(Sudam) # Keep only South American
# Save it
# st_write(obj = Sudam_pa_data, dsn = "WDPA_Sudam_terr.gpkg")

# Load PA vector
Sudam_pa_data <- st_read("data/04_WDPA_Sudam_terr.gpkg") |> 
  st_intersection(Chaco_biome)

## 0.2. Load raster habitat raster to be used in rescaling all rast --------

habitat_suitability_change <- rast("02_Habitat_loss/habitat_change.tif")

# 1. Representativity in Current climatic suitable areas ----------------------

# Load cureent climatic suitability raster umbralized
# umbralized means that above the threshold used, the pixel retains the 
# suitability values
current_clim_suit <- rast("01_ClimaticENM/ENM_Current_suitability.tif") |> 
  resample(habitat_suitability_change, method = "bilinear") |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = T))


qtm(current_clim_suit)

# Rasterize PA vector
Sudam_PA_rast <- rasterize(Sudam_pa_data, habitat_suitability_change, touches = F)
qtm(Sudam_PA_rast)

# Estimate SRI
SRI_current_clim <- sri_index(current_clim_suit)

# Estimate total and protected suitable area (TSA and PSA), and the percentage
# of suitable areas protected (PSAP)
TSA_current_clim <- expanse(current_clim_suit, unit = "km", transform = T)
PSA_current_clim <- area_idonea_protegida(current_clim_suit)
PSAP_current_clim <- (PSA_current_clim / TSA_current_clim) * 100

# Create a data frame to store representativeness metrics results:
Representativity_df <- data.frame(SRI = c(SRI_current_clim),
                                  TSA = c(TSA_current_clim[1,2]),
                                  PSA = c(PSA_current_clim[1,2]),
                                  PSAP = c(PSAP_current_clim[1,2]),
                                  row.names = c("Current_clim"))

# 2. Representativity in Current climatic suitable areas corrected by habitat loss -------

# Load habitat loss intensity layer
habitat_loss_intensity <- rast("03_Overlap/habitat_loss_intensity.tif")

# Correct suitability by habitat loss
current_clim_HL_suit <- current_clim_suit * (1 - habitat_loss_intensity)

# Estimate SRI
SRI_current_clim_HL <- sri_index(current_clim_HL_suit)

# Estimate total and protected suitable area (TSA and PSA), and the percentage
# of suitable areas protected (PSAP)
TSA_current_clim_HL <- expanse(current_clim_HL_suit, unit = "km", transform = T)
PSA_current_clim_HL <- area_idonea_protegida(current_clim_HL_suit)
PSAP_current_clim_HL <- (PSA_current_clim_HL / TSA_current_clim_HL) * 100

# Create a data frame to store representativeness metrics results:
Representativity_df_HL <- data.frame(SRI = c(SRI_current_clim_HL),
                                  TSA = c(TSA_current_clim_HL[1,2]),
                                  PSA = c(PSA_current_clim_HL[1,2]),
                                  PSAP = c(PSAP_current_clim_HL[1,2]),
                                  row.names = c("current_clim_HL"))
Representativity_df <- rbind(Representativity_df, Representativity_df_HL)

# 3. Representativity in Current climatic suitable areas corrected by habitat loss -------

# Load habitat loss intensity layer
hunting_intensity <- rast("data/HP_tolymat.tif") |> 
  resample(habitat_suitability_change, method = "bilinear")

# Correct suitability by habitat loss
current_clim_HP_suit <- current_clim_suit * (1 - hunting_intensity)

# Estimate SRI
SRI_current_clim_HP <- sri_index(current_clim_HP_suit)

# Estimate total and protected suitable area (TSA and PSA), and the percentage
# of suitable areas protected (PSAP)
TSA_current_clim_HP <- expanse(current_clim_HP_suit, unit = "km", transform = T)
PSA_current_clim_HP <- area_idonea_protegida(current_clim_HP_suit)
PSAP_current_clim_HP <- (PSA_current_clim_HP / TSA_current_clim_HP) * 100

# Create a data frame to store representativeness metrics results:
Representativity_df_HP <- data.frame(SRI = c(SRI_current_clim_HP),
                                  TSA = c(TSA_current_clim_HP[1,2]),
                                  PSA = c(PSA_current_clim_HP[1,2]),
                                  PSAP = c(PSAP_current_clim_HP[1,2]),
                                  row.names = c("current_clim_HP"))
Representativity_df <- rbind(Representativity_df, Representativity_df_HP)


# 4. Representativity in Current climatic suitable areas corrected by habitat loss and hunting pressure-------

# Correct suitability by habitat loss and hunting pressure
current_clim_HL_HP_suit <- current_clim_suit * (1 - habitat_loss_intensity) * (1 - hunting_intensity)

# Estimate SRI
SRI_current_clim_HL_HP <- sri_index(current_clim_HL_HP_suit)

# Estimate total and protected suitable area (TSA and PSA), and the percentage
# of suitable areas protected (PSAP)
TSA_current_clim_HL_HP <- expanse(current_clim_HL_HP_suit, unit = "km", transform = T)
PSA_current_clim_HL_HP <- area_idonea_protegida(current_clim_HL_HP_suit)
PSAP_current_clim_HL_HP <- (PSA_current_clim_HL_HP / TSA_current_clim_HL_HP) * 100

# Create a data frame to store representativeness metrics results:
Representativity_df_HL_HP <- data.frame(SRI = c(SRI_current_clim_HL_HP),
                                  TSA = c(TSA_current_clim_HL_HP[1,2]),
                                  PSA = c(PSA_current_clim_HL_HP[1,2]),
                                  PSAP = c(PSAP_current_clim_HL_HP[1,2]),
                                  row.names = c("current_clim_HL_HP"))

Representativity_df <- rbind(Representativity_df, Representativity_df_HL_HP)


# 5. Representativity in Future climatic suitabile areas -----------

## 5.1. SSP 1-2.6 ----------------------------------------------------------

# Load future climatic suitability raster umbralized

ssp126_clim_suit <- rast("01_ClimaticENM/ENM_126_suitability.tif") |> 
  resample(habitat_suitability_change, method = "bilinear")  |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = T))
qtm(ssp126_clim_suit)

# Estimate SRI
SRI_ssp126_clim <- sri_index(ssp126_clim_suit)

# Estimate total and protected suitable area (TSA and PSA), and the percentage
# of suitable areas protected (PSAP)

TSA_ssp126_clim <- expanse(ssp126_clim_suit, unit = "km", transform = T)
PSA_ssp126_clim <- area_idonea_protegida(ssp126_clim_suit)
PSAP_ssp126_clim <- (PSA_ssp126_clim / TSA_ssp126_clim) * 100

# Create a data frame to store representativeness metrics results:
Representativity_df_ssp126 <- data.frame(SRI = c(SRI_ssp126_clim),
                                  TSA = c(TSA_ssp126_clim[1,2]),
                                  PSA = c(PSA_ssp126_clim[1,2]),
                                  PSAP = c(PSAP_ssp126_clim[1,2]),
                                  row.names = c("ssp126_clim"))

Representativity_df <- rbind(Representativity_df, Representativity_df_ssp126)

## 5.2. SSP 3-7.0 ----------------------------------------------------------
# Load future climatic suitability raster umbralized
ssp370_clim_suit <- rast("01_ClimaticENM/ENM_370_suitability.tif") |> 
  resample(habitat_suitability_change, method = "bilinear") |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = T))
qtm(ssp370_clim_suit)
# Estimate SRI
SRI_ssp370_clim <- sri_index(ssp370_clim_suit)

# Estimate total and protected suitable area (TSA and PSA), and the percentage
# of suitable areas protected (PSAP)
TSA_ssp370_clim <- expanse(ssp370_clim_suit, unit = "km", transform = T)
PSA_ssp370_clim <- area_idonea_protegida(ssp370_clim_suit)
PSAP_ssp370_clim <- (PSA_ssp370_clim / TSA_ssp370_clim) * 100

# Create a data frame to store representativeness metrics results:
Representativity_df_ssp370 <- data.frame(SRI = c(SRI_ssp370_clim),
                                         TSA = c(TSA_ssp370_clim[1,2]),
                                         PSA = c(PSA_ssp370_clim[1,2]),
                                         PSAP = c(PSAP_ssp370_clim[1,2]),
                                         row.names = c("ssp370_clim"))

Representativity_df <- rbind(Representativity_df, Representativity_df_ssp370)


## 5.3. SSP 5-8.5 ----------------------------------------------------------
# Load future climatic suitability raster umbralized
ssp585_clim_suit <- rast("01_ClimaticENM/ENM_585_suitability.tif") |> 
  resample(habitat_suitability_change, method = "bilinear") |> 
  classify(rcl = matrix(c(0, NA), ncol = 2, byrow = T))

# Estimate SRI
SRI_ssp585_clim <- sri_index(ssp585_clim_suit)

# Estimate total and protected suitable area (TSA and PSA), and the percentage
# of suitable areas protected (PSAP)
TSA_ssp585_clim <- expanse(ssp585_clim_suit, unit = "km", transform = T)
PSA_ssp585_clim <- area_idonea_protegida(ssp585_clim_suit)
PSAP_ssp585_clim <- (PSA_ssp585_clim / TSA_ssp585_clim) * 100

# Create a data frame to store representativeness metrics results:
Representativity_df_ssp585 <- data.frame(SRI = c(SRI_ssp585_clim),
                                         TSA = c(TSA_ssp585_clim[1,2]),
                                         PSA = c(PSA_ssp585_clim[1,2]),
                                         PSAP = c(PSAP_ssp585_clim[1,2]),
                                         row.names = c("ssp585_clim"))

Representativity_df <- rbind(Representativity_df, Representativity_df_ssp585)


# Save data frame with metrics
write.csv(Representativity_df, file = "04_Representativity_PA/Representativity_metrics.csv", row.names = T)
Representativity_df <- read.csv("04_Representativity_PA/Representativity_metrics.csv")


# 6. Plot metrics -----------------------------

# Add SSP column to plot according to temporan scenarios
Representativity_df$SSP <- c("Current", "Current", "Current", "Current",
                             "SSP 1.2-6", "SSP 3.7-0", "SSP 5.8-5")
Representativity_df$Threat <- c("None", "HL", "HP", "HL + HP",
                             "None", "None", "None")
Representativity_df$Threat <- factor(Representativity_df$Threat,
                                     levels = c("None", "HL", "HP", "HL + HP"))

# 6.1. SRI ----------------------------------------------------------------
plot_SRI <- ggplot(Representativity_df, aes(x = SSP, y = SRI)) + 
                                            #group = GCM)) + 
  geom_point(aes(colour = Threat), size = 4) +
  # geom_line() +
  labs(
    x = "Scenarios", 
    y = 'Species Representation Index') + 
    # colour = "GCM") +
  scale_x_discrete(labels = c("Present" , "SSP 1.2-6", "SSP 3.7-0", "SSP 5.8-5")) +
  # scale_y_continuous(limits = (7, 10)) +
  theme_replace()
plot_SRI
# Lo exporto como pdf para editarlo en inkscape
ggsave("plots/SRI_plot.pdf", plot = plot_SRI)

# 6.2. TSA ----------------------------------------------------------------
plot_TSA <- ggplot(Representativity_df, aes(x = SSP, y = TSA)) + 
  #group = GCM)) + 
  geom_point(aes(colour = Threat), size = 4) +
  # geom_line() +
  labs(
    x = "Scenarios", 
    y = 'Total suitable area'~(km^2)) + 
  # colour = "GCM") +
  scale_x_discrete(labels = c("Present" , "SSP 1.2-6", "SSP 3.7-0", "SSP 5.8-5")) +
  # scale_y_continuous(limits = (7, 10)) +
  theme_replace()
plot_TSA
# Lo exporto como pdf para editarlo en inkscape
ggsave("plots/TSA_plot.pdf", plot = plot_TSA)

# 6.3. PSA ----------------------------------------------------------------
plot_PSA <- ggplot(Representativity_df, aes(x = SSP, y = PSA)) + 
  #group = GCM)) + 
  geom_point(aes(colour = Threat), size = 4) +
  # geom_line() +
  labs(
    x = "Scenarios", 
    y = 'Protected suitable area'~(km^2)) + 
  # colour = "GCM") +
  scale_x_discrete(labels = c("Present" , "SSP 1.2-6", "SSP 3.7-0", "SSP 5.8-5")) +
  # scale_y_continuous(limits = (7, 10)) +
  theme_replace()
plot_PSA
# Lo exporto como pdf para editarlo en inkscape
ggsave("plots/PSA_plot.pdf", plot = plot_PSA)

# 6.4. PSAP ----------------------------------------------------------------
plot_PSAP <- ggplot(Representativity_df, aes(x = SSP, y = PSAP)) + 
  #group = GCM)) + 
  geom_point(aes(colour = Threat), size = 4) +
  # geom_line() +
  labs(
    x = "Scenarios", 
    y = 'Protected sutable suitable area (%)') + 
  # colour = "GCM") +
  scale_x_discrete(labels = c("Present" , "SSP 1.2-6", "SSP 3.7-0", "SSP 5.8-5")) +
  # scale_y_continuous(limits = (7, 10)) +
  theme_replace()
plot_PSAP
# Lo exporto como pdf para editarlo en inkscape
ggsave("plots/PSAP_plot.pdf", plot = plot_PSAP)

qtm(current_clim_suit)

# Rasterize PA vector
Sudam_PA_rast <- rasterize(Sudam_pa_data, habitat_suitability_change, touches = F)
qtm(Sudam_PA_rast)

# Estimate SRI
SRI_current_clim <- sri_index(current_clim_suit)

# Estimate total and protected suitable area (TSA and PSA), and the percentage
# of suitable areas protected (PSAP)
TSA_current_clim <- expanse(current_clim_suit, unit = "km", transform = T)
PSA_current_clim <- area_idonea_protegida(current_clim_suit)
PSAP_current_clim <- (PSA_current_clim / TSA_current_clim) * 100

# Create a data frame to store representativeness metrics results:
Representativity_df <- data.frame(SRI = c(SRI_current_clim),
                                  TSA = c(TSA_current_clim[1,2]),
                                  PSA = c(PSA_current_clim[1,2]),
                                  PSAP = c(PSAP_current_clim[1,2]),
                                  row.names = c("Current_clim"))
