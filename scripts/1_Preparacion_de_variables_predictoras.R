### Future SDM of Tolypeutes matacus
### Paquetes usados y directorio -----------------------------------------------------------
library(dismo)
library(ENMeval)
library(maxnet)
library(raster)
library(rworldxtra)
library(sf)
library(tidyverse)
library(kuenm)

setwd("D:/MEGAsync/Modelado/TolyMat_Futuro")

### Ocurrencias -----------------------------------------------------------------------
OccTolyMat <- read.csv("D:/MEGAsync/Modelado/Datos_Presencia/Final_Paper_TolyMat/Occs_QGIS.csv",
                       sep = ";" )

# Genero un set de datos solo Con datos Actuales
OccTolyMat_Actual <- OccTolyMat %>% dplyr::filter(Temporalidad == "Actual")

#Transformamos en SF
OccTolyMat_SF <- OccTolyMat %>% st_as_sf(coords = c(7, 6), 
                                       crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

### Area M -----------------------------------------------------------------------------
AreaM <- shapefile("D:/MEGAsync/Modelado/Areas M/AreaM_paper.shp")
AreaM_SF <- read_sf("D:/MEGAsync/Modelado/Areas M/AreaM_paper.shp")
st_transform(AreaM_SF, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

### Mapas base
# Politico
data("countriesHigh")
Politico <- countriesHigh %>% st_as_sf() %>% st_make_valid() %>% st_crop(AreaM_SF)

# Bioregiones Morrone 2014
Bioregion_SF <- read_sf("D:/BackUp-Ale/Capas/Capas_Biologicos/BioRegions_Morrone2014/Lowenberg_Neto_2014.shp", 
                        crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

st_transform(Bioregion_SF, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

Bioregion_SF <- st_crop(Bioregion_SF, AreaM_SF)

# Ecoregiones Olson 2001
Ecoregion_SF <- read_sf("D:/BackUp-Ale/Capas/Capas_Biologicos/Ecorregiones_Olson/wwf_terr_ecos.shp") 
st_transform(Ecoregion_SF, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
EcoregionCrop_SF <- st_crop(Ecoregion_SF, AreaM_SF)

# Ploteo lo puntos de acuerdo a su temporalidad junto con el area M.
ggplot() + 
  geom_sf(data = Bioregion_SF, aes(fill = Province_1)) +# mapa base politico, se puede cambiar por el de bioregiones
  scale_fill_grey(na.translate = FALSE) +
  geom_sf(data = AreaM_SF, fill = NA, color = "blue") +
  geom_point(data = OccTolyMat, 
             aes(x = Longitud,
                 y = Latitud, 
                 shape = factor(Temporalidad)), 
             size = 2) + 
  facet_wrap(~Temporalidad) + # Si dejo este termino me hace un grafico para cada tiempo
  labs(shape = "Temporalidad") +
  ylab("") + xlab("") +
  theme_bw()
### Variables Ambientales --------------------------------------------------------------
# Genero archivo con links a las capas cortadas por sudam
# Climaticas
paths_Actual <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Current/Actual_2.5min/Sudam",
                           pattern = "*.asc$",full.names = TRUE)
paths_RCP26_CCSM <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Future_2.5m/cc26bi70/Sudam",
                               pattern = "*.asc$",full.names = TRUE)
paths_RCP45_CCSM <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Future_2.5m/cc45bi70/Sudam",
                               pattern = "*.asc$",full.names = TRUE)
paths_RCP60_CCSM <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Future_2.5m/cc60bi70/Sudam",
                               pattern = "*.asc$",full.names = TRUE)
paths_RCP85_CCSM <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Future_2.5m/cc85bi70/Sudam",
                               pattern = "*.asc$",full.names = TRUE)
paths_RCP26_MIROC <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Future_2.5m/mr26bi70/Sudam",
                                pattern = "*.asc$",full.names = TRUE)
paths_RCP45_MIROC <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Future_2.5m/mr45bi70/Sudam",
                                pattern = "*.asc$",full.names = TRUE)
paths_RCP60_MIROC <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Future_2.5m/mr60bi70/Sudam",
                                pattern = "*.asc$",full.names = TRUE)
paths_RCP85_MIROC <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Future_2.5m/mr85bi70/Sudam",
                                pattern = "*.asc$",full.names = TRUE)
paths_RCP26_MPI <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Future_2.5m/mp26bi70/Sudam",
                              pattern = "*.asc$",full.names = TRUE)
paths_RCP45_MPI <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Future_2.5m/mp45bi70/Sudam",
                              pattern = "*.asc$",full.names = TRUE)
paths_RCP85_MPI <- list.files(path = "D:/BackUp-Ale/Capas/Variables_climaticas/WorldClim_1.4/Future_2.5m/mp85bi70/Sudam",
                              pattern = "*.asc$",full.names = TRUE)
###### Genero raster stack de cada set de variables
Actual <- raster::stack(paths_Actual)

LIG <- raster::stack(paths_LIG)

RCP26_cc <- raster::stack(paths_RCP26_CCSM)
RCP45_cc <- raster::stack(paths_RCP45_CCSM)
RCP60_cc <- raster::stack(paths_RCP60_CCSM)
RCP85_cc <- raster::stack(paths_RCP85_CCSM)
RCP26_mr <- raster::stack(paths_RCP26_MIROC)
RCP45_mr <- raster::stack(paths_RCP45_MIROC)
RCP60_mr <- raster::stack(paths_RCP60_MIROC)
RCP85_mr <- raster::stack(paths_RCP85_MIROC)
RCP26_mp <- raster::stack(paths_RCP26_MPI)
RCP45_mp <- raster::stack(paths_RCP45_MPI)
RCP85_mp <- raster::stack(paths_RCP85_MPI)

Soil <- raster::stack(paths_Soil)
SoilLGM <- raster::stack(paths_SoilLGM)

###Resamplear capas diferentes ####
SoilResampled <- resample(Soil, Actual, method="bilinear" )#ver si anda la resolucion
plot(SoilResampled)
?dim
dim(HTMccsm) == dim(SoilResampled)
res(Actual)
res(SoilResampled)

SoilLGMResampled <- resample(SoilLGM, LGMccsm, method="bilinear" )#ver si anda la resolucion
plot(SoilLGMResampled)
res(HTMccsm) == res(SoilResampled)
res(Actual)
res(SoilResampled)


LIGResampled <- resample(LIG, Actual, method="bilinear" )#ver si anda la resolucion
plot(LIGResampled)
res(LIGResampled) == res(SoilResampled)
res(Actual)
res(SoilResampled)

### Genero un stack con todas las variables?? Darles CRS? #####
crs(Actual) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(HTMccsm) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(HTMmiroc) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(HTMmpi) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(LGMccsm) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(LGMmiroc) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(LGMmpi) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(LIGResampled) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(SoilResampled) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(SoilLGMResampled) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(RCP26_cc) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(RCP45_cc) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(RCP60_cc) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(RCP85_cc) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(RCP26_mr) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(RCP45_mr) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(RCP60_mr) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(RCP85_mr) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(RCP26_mp) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(RCP45_mp) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(RCP85_mp) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

nlayers(SoilResampled) == nlayers(RCP26_cc)

CasiTodos <- c(Actual, 
               HTMccsm, HTMmiroc, HTMmpi, 
               LGMccsm, LGMmiroc, LGMmpi, LIGResampled, 
               SoilResampled, SoilLGMResampled)

CasiTodosFuturo <- c(Actual, 
                     RCP26_cc, RCP45_cc, RCP60_cc, RCP85_cc,
                     SoilResampled, SoilLGMResampled)


Todos <- c(Actual, 
           HTMccsm, HTMmiroc, HTMmpi, 
           RCP26_cc, RCP45_cc, RCP60_cc, RCP85_cc,
           RCP26_mr, RCP45_mr, RCP60_mr, RCP85_mr,
           RCP26_mp, RCP45_mp, RCP85_mp, 
           LGMccsm, LGMmiroc, LGMmpi, LIGResampled, 
           SoilResampled, SoilLGMResampled)

StackCasiAll <- raster::stack(CasiTodos)
StackCasiAllFuture <- raster::stack(CasiTodosFuturo)

### Cargar mascara de Area M #######
Mask <- shapefile("D:/MEGAsync/Modelado/Areas M/AreaM_paper.shp")
plot(Mask)

### Cropear ####
StackCasiAllCropped <- raster::crop(x = StackCasiAllFuture,y = Mask)

### Maskear ####
StackAllMasked <- raster::mask(StackCasiAllCropped,Mask)

### Guardar capas de salida #####

dir.create("TodasFuturo")
#Creo los rasters en mi carpeta
lapply(names(StackAllMasked), 
       function(x){
         writeRaster(StackAllMasked[[x]], 
                     paste0("TodasFuturo/",x,".asc"),
                     overwrite=TRUE)})