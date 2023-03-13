# Cargo paquetes a usar
library(sf)
library(tidyverse)
library(terra)
library(tmap) #paquete para hacer mapas

# Capas geográficas -------------------------------------------------------

# 1.1. Países
Sudam = rnaturalearth::ne_countries(scale = 10, continent = "south america", returnclass = "sf")
Argentina = rnaturalearth::ne_countries(scale = 10, country = "argentina", returnclass = "sf")


# Explicación del SRI index -----------------------------------------------

# SRI= sum (valor de proteccion del pixel * valor de idoneidad)/sum (valor de idoneidad)
# Entonces para su calculo necesito:
### 1. Raster con valor de protección del pixel
### 2. Raster con valor de idoneidad del pixel
### 3. Calculo de SRI index

# 1. Raster con valor de proteccion (Ri) ---------------------------------------

## 1a. Descargo info espacial de PAs -------------------------------------------

# Descargo archivo areas protegidas
# Creo la carpeta para guardar las capas
dir.create("WDPA")
# Las descargo, buscar en la pagina de la WDPA, y al elegir la descarga en vez 
# de hacer click en download, pongo copiar enlace de descarga y lo pego aqui
download.file(url = "https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_Apr2021_Public_shp.zip",
              destfile = "WDPA/wdpa_Apr2021_shp.zip")
unzip(zipfile = "WDPA/wdpa_Apr2021_shp.zip", exdir = "WDPA") # descomprimo archivo

# Esta base de datos trae poligonos y tb puntos. Los puntos representan areas en 
# donde no hay informacion sobre su demarcacion. Por eso no seran tenidas en 
# cuenta en este analisis.

# El shape global viene dividido en tres partes, y todas comprimidas, entonces
# descomprimo cada uno de los tres shapes de poligonos
dir.create("WDPA/wdpa_Apr2021_Public_shp_0")
unzip(zipfile = "WDPA/wdpa_Apr2021_Public_shp_0.zip", exdir = "WDPA/wdpa_Apr2021_Public_shp_0")

# Descargo PAs de Argentina
dir.create("Arg_PA")
download.file(url = "https://dnsg.ign.gob.ar/apps/api/v1/capas-sig/Geodesia+y+demarcaci%C3%B3n/L%C3%ADmites/linea_de_limite_070114/shp",
              destfile = "Arg_PA/Arg_PA.zip")
unzip(zipfile = "Arg_PA/Arg_PA.zip")

## 1b. Cargo info espacial de PAs en R -------------------------------------
global_pa_data_0 = st_read(dsn = "WDPA//wdpa_Apr2021_Public_shp_0/WDPA_Apr2021_Public_shp-polygons.shp") # cargo en R
# Achico las columnas del vector para que se mas facil trabajarlo
head(global_pa_data_0)
global_pa_data_0b = global_pa_data_0[, c("NAME", "IUCN_CAT", "MARINE", "GIS_AREA", "ISO3" )] # subset columns by name

global_pa_data_0c = global_pa_data_0b %>% 
  filter(MARINE == 0) # Saco las PAs marinas

# Cargo un vector de sudamerica para usarlo para cortar los datos de PAs
library(spData) # paquete con vector con limites politicos mundiales
Sudam = spData::world %>% filter(continent == "South America")
Sudam_pa_data = global_pa_data_0c[Sudam, ]
# Guardo el archivo
st_write(obj = Sudam_pa_data, dsn = "WDPA_Sudam_terr.gpkg")
# En el futuro puedo cargar el archivo directamente
Sudam_pa_data = st_read("3_WDPA_Sudam_terr.gpkg")

# Mapeo PAs
map_sudam = tm_shape(Sudam) + tm_borders()
map_sudam
map_sudam_pa = map_sudam + tm_shape(Sudam_pa_data) + tm_fill(col = "ISO3")
map_sudam_pa

## 1c. Rasterizar vector de PAs usando el raster de idoneidad --------
# Cargo el raster de idoneidad
library(raster)
Current_Si_raster = raster::raster("C:/Users/ale_f/OneDrive/Documentos/Doctorado/ENM/TolyMat_Futuro/Future_WC2.1_Soil/2_Final_Models/M_0.1_F_lq_Set_1_EC/Tolypeutes_matacus_Actual_median.asc")
# Rasterizo Sudam_pa_data
## Raster binario con presencia (=1) y ausencia (=0) de PAs
Ri_bin_raster = fasterize::fasterize(Sudam_pa_data, Current_Si_raster, field = NULL, background = 0)
values(Ri_bin_raster)
plot(Ri_bin_raster)

# 2. Raster con valor de idoneidad (Si) ---------------------------------------
# Cargo las raster de idoneidad de Maxent
paths_Si_rasters <- list.files(path = "C:/Users/ale_f/OneDrive/Documentos/Doctorado/ENM/TolyMat_Futuro/Future_WC2.1_Soil/2_Final_Models/M_0.1_F_lq_Set_1_EC",
                           pattern = "*_median.asc$",full.names = TRUE)
Si_rasters <- raster::stack(paths_Si_rasters)

# Remuevo el raster de suitability actual que esta repetido
names(Si_rasters)
help("dropLayer")
Si_rasters = raster::dropLayer(Si_rasters, 14)

# 3. Calculo SRI index ----------------------------------------------------

# One threshold SRI, Raster de proteccion binarizado, y Raster de idoneidad continuo
# SRI= sum (valor de proteccion del pixel * valor de idoneidad)/sum (valor de idoneidad)
# Current one-thresholded SRI

# Funcion que calcula el SRI index
sri_index = function(x) {
  Prod_RiSi <- Ri_bin_raster*x # Primero multiplicar el Si raster con el Ri raster.
  Sum_RixSi <- cellStats(Prod_RiSi, stat = sum) # Sumatoria de los valores de celdas de pRODUCTO sI_rI
  Sum_Si = cellStats(x, stat = sum) # Sumatoria de valores de pixeles del raster Si
  SRI = Sum_RixSi/Sum_Si 
  SRI                       # la última línea es lo que devolverá la función
}   

#Calculo SRI en diferentes escenarios climaticos
scenarios_SRI = lapply(names(Si_rasters), 
       function(x){
         Prod_RiSi <- Ri_bin_raster*Si_rasters[[x]] # Primero multiplicar el Si raster con el Ri raster.
         Sum_RixSi <- cellStats(Prod_RiSi, stat = sum) # Sumatoria de los valores de celdas de pRODUCTO sI_rI
         Sum_Si = cellStats(Si_rasters[[x]], stat = sum) # Sumatoria de valores de pixeles del raster Si
         SRI = Sum_RixSi/Sum_Si 
         SRI      
         })


# 4. Ploteo de los valores de SRI -----------------------------------------
library(ggplot2)
# Genero el data frame, para dsps graficar el cambio de SRI index, segun ssp, y GCM
SRI_values = as.vector(unlist(scenarios_SRI))

# Remover "Tolypeutes_matacus" y "_median" del vector Scenarios
Scenarios = names(Si_rasters) %>% 
  str_remove("Tolypeutes_matacus_") %>% 
  str_remove("_median")

SRI_df <- data.frame(
  "Scenario" = Scenarios, 
  "SRI_index" = SRI_values)

SRI_df1 = SRI_df %>% 
  separate(Scenario, c("GCM", "SSP"), sep = "_ssp")
SRI_df1[1,1]= "Present"

# Como me aparece el eje X desordenado (Current al final), reordeno
SRI_df1$SSP <- factor(SRI_df1$SSP,levels = c("Present", "126", "245", "370", "585"))

# Guardo el archivo con los valores de SRI obtenidos
write.csv(SRI_df1, file = "SRIindex.csv")
SRI_df1 = read.csv(file = "SRIindex.csv")
# Como me aparece el eje X desordenado (Current al final), reordeno
SRI_df1$SSP <- factor(SRI_df1$SSP,levels = c("Present", "126", "245", "370", "585"))

base <- ggplot(SRI_df1, aes(x = SSP, y = SRI_index, group = GCM))
base + 
  geom_point(aes(colour = GCM), size = 3) +
  geom_line(aes(colour = GCM)) +
  labs(
    x = "Climatic scenario", 
    y = "SRI", 
    colour = "GCM") +
  scale_x_discrete(labels = c("Present" , "SSP 1.2-6", "SSP 2.4-5", "SSP 3.7-0", "SSP 5.8-5")) +
  scale_fill_brewer(palette = "Set1") +
  theme_replace()
