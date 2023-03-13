# Paquetes ----------------------------------------------------------------

library(raster)
library(sf)
library(tmap)
library(tmaptools)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(terra)

# Uso funcion que robe de "https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/"
# pero la modifico para primero calcular el valor de thrshold que quiero:
# Esta funcion tiene para calcular tres tipos de thrEsholds: MTP, P5 Y P10.
thresh = function(sdm, occs, type = "mtp"){
  occPredVals <- raster::extract(sdm, occs) # Estraigo valores de idoneidad de presencias
  if(type == "mtp"){
    thresh <- min(na.omit(occPredVals))
  } else if(type == "p5"){
    if(length(occPredVals) < 10){
      p5 <- floor(length(occPredVals) * 0.95)
    } else {
      p5 <- ceiling(length(occPredVals) * 0.95)
    }
    thresh <- rev(sort(occPredVals))[p5]
  } else if(type == "p10"){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  return(thresh)
}
# Funcion para binarizar los mapas de acuerdo al umbral calculado previamente
sdm_threshold <- function(sdm, Umbral, binary = FALSE){
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < Umbral] <- 0
  if(binary){
    sdm_thresh[sdm_thresh >= Umbral] <- 1
  }
  return(sdm_thresh)
}

# Capas geograficas -------------------------------------------------------
# 1.1. Países
Sudam = ne_countries(scale = 10, continent = "south america", returnclass = "sf")
Argentina = ne_countries(scale = 10, country = "argentina", returnclass = "sf")

# 1.2. Provincias
ProvSudam = ne_states(iso_a2 = c("AR", "BO", "PY", "BR", "CL", "UY", "PE", "EC", "VE", "CO", "GY", "SR", "FR"),
                      returnclass = "sf")
ProvArg = ne_states(iso_a2 = c("AR"),
                    returnclass = "sf")

## Ecorregiones (Olson et al., 2001)
# Descargado de https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world 

# Cargo mapa de bioregiones de sudamerica de Morrone 2014
Bioregions = st_read("C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/Capas_biologicas/BioRegions_Morrone2014/Lowenberg_Neto_2014.shp")
Bioregions_geo = st_transform(Bioregions,crs = 4326)
Bioregion_geo = Bioregions %>%  
  filter(Province_1=="Chacoan province" | Province_1=="Monte province" | Province_1=="Pampean province" ) %>%
  st_transform(crs = 4326)

Chaco_Bioregion = Bioregions %>%  
  filter(Province_1=="Chacoan province") %>%
  st_transform(crs = 4326)

ChacoArg_Bioregion = st_intersection(Chaco_Bioregion, Argentina)

# Mapa base
mapa_base = tm_shape(ProvSudam, bbox = Chaco_Bioregion) + tm_borders(col = "#252525" , lty = "dotted") + tmap_options(check.and.fix = TRUE)+
  tm_shape(Sudam, bbox = Chaco_Bioregion) +  tm_borders(col = "#252525", lwd = 2)

# 1. Amenazas --------------------------------------------------

# 1.1. Indice de caza --------------------------------------------------------

# Índice de presión de caza en el Chaco estimado por Romero-Muñoz, 2020.
HP_Romero = raster("Paper_Hunting/Tolypeutes_matacus_2015_pu.tif")
HP_Romero_geo <- projectRaster(HP_Romero, crs = 4326, method = "bilinear")
HP_Romero_geo2 <- abs(HP_Romero_geo -1)
HP_Romero_geo_resampl <- resample(HP_Romero_geo2, Current_Si_raster, method= "ngb")
HP_Romero_geo_resampl1 <- crop(HP_Romero_geo_resampl, Chaco_Bioregion)
HP_Romero_geo_resampl2 <- mask(HP_Romero_geo_resampl1, Chaco_Bioregion)

# Plot
HP_Mapa =  tm_shape(Sudam, bbox = Chaco_Bioregion) + tm_polygons(col = "lightgrey") + 
  tm_shape(HP_Romero_geo_resampl2) + tm_raster(style = "cont",
                                palette =viridisLite::viridis(20),
                                title= "HP",
                                legend.format =list(text.separator="-"),
                                legend.reverse = T) +
  tm_shape(Sudam) + tm_borders(lwd = 2) +
  tm_shape (ProvSudam) + tm_borders(col = "black") + 
  tm_layout(legend.position = c("right","bottom")) + 
  tmap_options(check.and.fix = TRUE)
HP_Mapa
tmap_save(HP_Mapa, "HP_mapa.pdf")

# 1.2. Destrucción de habitat ---------------------------------------------

## Descargo mapa de land cover de Copernicus ("https://zenodo.org/record/3939050/files/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif?download=1")
# Cargo el tif
Land_cover_global = raster::raster("Land_cover/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif")
# Lo corto por sudam
Land_cover_chacoCropped = raster::crop(Land_cover_global, ChacoEcoregion)#genero un raster con las mismas dimensiones q la mascara
Land_cover_chaco = raster::mask(Land_cover_chacoCropped, ChacoEcoregion)#convierte lo que esta fuera del shape a NA u otro valor (updatevalue =)
writeRaster(Land_cover_chaco, "Land_cover/LC_Chaco.tif", overwrite = T)
Land_cover_chaco = raster("Land_cover/LC_Chaco.tif")

# Pensar en valores de asignacion a cada cobertura de suelo, que afecten a la especie
# Por ejemplo, Forests = 0, Shrubs, Cropland, Herbaceous vegetation, Herbaceous wetland,
# urban, snow and ice.
rcl = matrix(c(0, NA, #No input data available
               111, 0, #Closed forest, evergreen needle leaf
               112, 0, #Closed forest, deciduous needle leaf
               113, 0, #Closed forest, evergreen, broad leaf
               114, 0, #Closed forest, deciduous broad leaf
               115, 0, #Closed forest, mixed
               116, 0, #Closed forest, unknown
               121, 0, #Open forest, evergreen needle leaf
               122, 0, #Open forest, deciduous needle leaf
               123, 0, #Open forest, evergreen broad leaf
               124, 0, #Open forest, deciduous broad leaf
               125, 0, #Open forest, mixed Open
               126, 0, #Open forest, unknown
               20,  0, #Shrubs
               30, 0, #Herbaceous vegetation
               90, 0, #Herbaceous wetland
               100, 0, #Moss and lichen
               60, 0, #Bare / Sparse vegetation
               40, 0.75, #Cropland
               50, 1, #Built-up or urban
               70, 0, #Snow & ice
               80, 0, #Permanent water bodies
               200, NA ), #Open sea
              ncol = 2, byrow = TRUE)

LC_Chaco_Reclas = reclassify(Land_cover_chaco, rcl)
# Cambiar resolucion para que matchee la de los modelos.
Current_Si_raster = raster::raster("C:/Users/ale_f/OneDrive/Documentos/Doctorado/ENM/TolyMat_Futuro/Future_WC2.1_Soil/2_Final_Models/M_0.1_F_lq_Set_1_EC/Tolypeutes_matacus_Actual_median.asc")
PerdidaHab <- resample(LC_Chaco_Reclas, Current_Si_raster, method="ngb")# Ver si uso metodo bilinear o sin metodo
crs(PerdidaHab) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
writeRaster(PerdidaHab, "PerdidaHab_Chaco.tif", overwrite = T)
PerdidaHab = raster("4_PerdidaHab_Chaco.tif")
plot(PerdidaHab)

# Uso mapa de desmontes obtenido de http://monitoreodesmonte.com.ar/
# Colección 8.0
# Enlace de descarga: "http://monitoreodesmonte.com.ar/descargar.php?id=19"
DesmonteChacoArg = st_read("C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/Capas_biologicas/Desmontes_Chaco/Coleccion_8.0_Argentina-1976-2019.shp")
DesmonteChacoArge_geo = st_transform(DesmonteChacoArg, crs = 4326)

help("rasterize")
DesmonteChacoArg_rast = fasterize::fasterize(DesmonteChacoArge_geo, 
                                  PerdidaHab, 
                                  field = NULL,
                                  background = 0)
plot(DesmonteChacoArg_rast)

# Uso tercera fuente de datos de Global Forest Watch
# Raster que muestra pixel donde se perdieron cobertura arborea clasificados segun la causa

GFW_tree_cover_loss <- raster("C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/Capas_biologicas/GlobalForestWatch/Tree_cover_loss_dominant_driver/Goode_FinalClassification_19_05pcnt_prj.tif")
GFW_tree_cover_loss_Chaco <- crop(GFW_tree_cover_loss, Chaco_Bioregion)
GFW_HL_Chaco <- mask(GFW_tree_cover_loss_Chaco, Chaco_Bioregion)
qtm(GFW_HL_Chaco)
rclGFW = matrix(c(0, 0,
               1, 1, 
               2, 1, 
               3, 0, 
               4, 1, 
               5, 1, 
               NA, NA),
             ncol = 2, byrow = TRUE)

GFW_HL_Chaco_reclassified <- reclassify(GFW_HL_Chaco, rclGFW)
qtm(GFW_HL_Chaco_reclassified)
GFW_HL_Chaco_resmpled <- resample(GFW_HL_Chaco_reclassified, Current_Si_raster, method="ngb")# Ver si uso metodo bilinear o sin metodo
crs(GFW_HL_Chaco_resmpled) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
qtm(GFW_HL_Chaco_resmpled)

# Sumo raster de desmontes + raster de perdida de hab obtenido de copernicus
PH_para_reclasificar = PerdidaHab + DesmonteChacoArg_rast + GFW_HL_Chaco_resmpled
m = c(-1, 0.1, 0, 0.1, 3.5, 1) 
rcl = matrix(m, ncol = 3, byrow = T)
HL = reclassify(PH_para_reclasificar, rcl)
qtm(HL)

# Corto usando la bioregión del Chaco de Argentina
PH_crop = crop(PH, ChacoArg_Bioregion) 
HL_rast = mask(PH_crop, ChacoArg_Bioregion)

# Plot
HL_Mapa =  tm_shape(Sudam, bbox = Chaco_Bioregion) + tm_polygons(col = "lightgrey") + 
  tm_shape(HL) + tm_raster(style = "fixed",
                                breaks = c(-0,1,2),
                                labels = c("Unmodified sites","Human modified sites"),
                                palette = viridisLite::viridis(20),
                                title= "HL",
                                legend.format =list(text.separator="-"),
                                legend.reverse = T) +
  tm_shape(Sudam) + tm_borders(lwd = 2) +
  tm_shape (ProvSudam) + tm_borders(col = "black") + 
  tm_layout(legend.position = c("right","bottom")) + 
  tmap_options(check.and.fix = TRUE)
HL_Mapa
tmap_save(HL_Mapa, "HL_mapa.pdf")

# 1.3. Cambio climático -------------------------------------------------
## CambioClim = 1 - Prob en el pixel i promediada entre las diferentes escenarios y GCMs
## Indice de proteccion futura = Proteccion en el pixel en N GCMs
# Cargo mapas de idoneidad futurA
Suitab_future_paths = list.files(path = "Future_WC2.1_Soil/2_Final_Model_Stats/Statistics_EC",
                                 pattern = "*med.tif$",full.names = TRUE)
Suitab_future_stack <- raster::stack(Suitab_future_paths)
Suitab_future_stack = dropLayer(Suitab_future_stack, 1)
Suitab_future_mean = overlay(Suitab_future_stack, fun = "mean")
CambioClim = 1 - Suitab_future_mean
crs(CambioClim) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
writeRaster(CambioClim, "CambioClim.tif", overwrite = T)
CambioClim = raster("CambioClim.tif")
CambioClim_chacoCropped = raster::crop(CambioClim, ChacoEcoregion)#genero un raster con las mismas dimensiones q la mascara
CambioClim_chaco = raster::mask(CambioClim_chacoCropped, ChacoEcoregion)#convierte lo que esta fuera del shape a NA u otro valor (updatevalue =)

# Corto usando la bioregión del Chaco de Argentina
DF_crop = crop(CambioClim_chaco, ChacoArg_Bioregion) 
DF_rast = mask(DF_crop, ChacoArg_Bioregion)

# Plot
UF_Mapa =  tm_shape(Sudam, bbox = Chaco_Bioregion) + tm_polygons(col = "lightgrey") + 
  tm_shape(CambioClim_chaco) + tm_raster(style = "cont",
                                palette =viridisLite::viridis(20),
                                title= "UF",
                                legend.format =list(text.separator="-"),
                                legend.reverse = T) +
  tm_shape(Sudam) + tm_borders(lwd = 2) +
  tm_shape (ProvSudam) + tm_borders(col = "black") + 
  tm_layout(legend.position = c("right","bottom")) + 
  tmap_options(check.and.fix = TRUE)
UF_Mapa
tmap_save(UF_Mapa, "UF_mapa.pdf")


# 1.4. Indice de amenaza --------------------------------------------------
IA_rast = (CambioClim_chaco + HP_Romero_geo_resampl2 + HL )/3  

# Plot
IA_Mapa =  tm_shape(Sudam, bbox = Chaco_Bioregion) + tm_polygons(col = "lightgrey") + 
  tm_shape(IA_rast) + tm_raster(style = "cont",
                                palette =viridisLite::viridis(20),
                                title= "TI",
                                legend.format =list(text.separator="-"),
                                legend.reverse = T) +
  tm_shape(Sudam) + tm_borders(lwd = 2) +
  tm_shape (ProvSudam) + tm_borders(col = "black") + 
  tm_layout(legend.position = c("right","bottom")) + 
  tmap_options(check.and.fix = TRUE)
IA_Mapa
tmap_save(IA_Mapa, "TI_mapa.pdf")

# Mapa final
Mapafinal = tmap_arrange(Caza_Map, PerdidaHab_Map, CambioClim_Map, Amenaza_Map, 
                         ncol = ,
                         nrow = 1)
Mapafinal