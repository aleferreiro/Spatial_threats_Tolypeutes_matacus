# Paquetes y wd ------

setwd("/Future_WC2.1_Soil")

library(kuenm)
#Para graficar
library(ggplot2)
library(ggspatial)
library(ggforce)
library(ggrepel)
library(scales)
library(raster)
library(rworldxtra)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

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
# Cargar archivo en R
Ecoregion = st_read("C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/Capas_biologicas/Ecorregiones_Olson/wwf_terr_ecos.shp")
ChacoEcoregion = Ecoregion %>% filter(ECO_NAME == "Dry Chaco" |ECO_NAME == "Humid Chaco") 

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

AreaM <- shapefile("D:/MEGAsync/Modelado/Areas M/AreaM_paper.shp")
AreaM_SF <- read_sf("C:/Users/ale_f/OneDrive/Documentos/Doctorado/ENM/TolyMat_Past/Areas M/AreaM_paper.shp")


# Modelado en kuenm
# 1. Modelos candidatos ---------------------------------------------------------------------------------------------------------
occ_joint <- "Sp_joint.csv"
occ_tra <- "Sp_train.csv"
M_var_dir <- "2_M_variables"
batch_cal <- "2_Candidate_models"
out_dir <- "2_Candidate_Models"
reg_mult <- c(0.1,0.5,1,2,5,8,10)
f_clas <- c("no.t.h") # alternativamente "all"
args <- "maximumbackground=10000"  # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
# note that some arguments are fixed in the function and should not be changed
maxent_path <- "D:/MEGAsync/Modelado/KUENM/." # path de maxent
wait <- FALSE
run <- TRUE

#Inicio el proceso de calibracion de los modelos
help(kuenm_cal)
kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
          maxent.path = maxent_path, wait = wait, run = run)


# 2. Evaluacion de modelos candidatos --------------------------------------------------------------------
occ_test <- "Sp_test.csv"
out_eval <- "2_Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- "OR_AICc"
paral_proc <- FALSE 

cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                        out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                        kept = kept, selection = selection, parallel.proc = paral_proc)

# 3. Modelos finales ---------------------------------------------------------------------------------
batch_fin <- "2_Final_models"
mod_dir <- "2_Final_Models"
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- TRUE
out_format <- "logistic"#buscar paper que hable de que salida es mas correcta
project <- TRUE
G_var_dir <- "2_G_variables"
ext_type <- "all" # para extrapolacion "all" PORQUE ME LO PIDE AHORA 
write_mess <- TRUE
write_clamp <- TRUE
wait1 <- FALSE
run1 <- TRUE
args <- "maximumbackground=10000"  

kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir =G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, 
          maxent.path = maxent_path, args = args, wait = wait1, run = run1)


# 4. MOP -------------------------------------------------------------------------------------------------------------
sets_var <- "Set_3" # a vector of various sets can be used
out_mop <- "MOP_results"
percent <- 5
paral <- F # make this true to perform MOP calculations in parallel, recommended
# only if a powerfull computer is used (see function's help)
# Two of the variables used here as arguments were already created for previous functions
kuenm_mmop(G.var.dir = G_var_dir, M.var.dir = M_var_dir, sets.var = sets_var, out.mop = out_mop,
           percent = percent, parallel = F, is.swd = FALSE)


# 5. Descriptive statistics of results --------------------------------------------------------------------------
help(kuenm_modstats)

spname <- "Tolypeutes_matacus"
modstats <- "2_Final_Model_Stats"
proj_scenarios <- c("Actual", 
                    "BCC_CSM2_MR_ssp126", "BCC_CSM2_MR_ssp245", "BCC_CSM2_MR_ssp370", "BCC_CSM2_MR_ssp585",
                    "CanESM5_ssp126", "CanESM5_ssp245", "CanESM5_ssp370", "CanESM5_ssp585", 
                    "IPSL_CM6A_LR_ssp126", "IPSL_CM6A_LR_ssp245", "IPSL_CM6A_LR_ssp370", "IPSL_CM6A_LR_ssp585",
                    "MIROC_ES2L_ssp126", "MIROC_ES2L_ssp245", "MIROC_ES2L_ssp370", "MIROC_ES2L_ssp585",
                    "MRI_ESM2_ssp126", "MRI_ESM2_ssp245", "MRI_ESM2_ssp370", "MRI_ESM2_ssp585")

kuenm_modstats(sp.name = spname, fmod.dir = mod_dir, format = "asc", 
               project = TRUE, 
               statistics = c("med", "range"), 
               replicated = TRUE,
               proj.scenarios = proj_scenarios, 
               ext.type = "EC", out.dir = modstats)

# 6. Cambios en las proyecciones-----------------------------------------------------------------------
help(kuenm_projchanges)

curpatt <- "Actual"
emi_scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")
clim_models <- c("BCC_CSM2_MR", "CanESM5", "IPSL_CM6A_LR", "MIROC_ES2L", "MRI_ESM2")

kuenm_projchanges(occ = occ_joint, fmod.stats = modstats, 
                  threshold = 5, 
                  emi.scenarios =  emi_scenarios,  
                  clim.models = clim_models,
                  ext.type = "EC", out.dir = "2_Projection_changes")

## 6a. Grafico comparacion binaria en cada escenario -----------------------------------------------------------

scenario126 <- rep("TSC 1-2.6", 323140)
ProjChanges_ssp126 <- raster("Future_WC2.1_Soil/2_Projection_changes/Changes_EC/Period_1/Scenario_ssp126/binary_comparison.tif") 
ProjChanges_ssp126b <- ProjChanges_ssp126 %>% as("SpatialPixelsDataFrame") %>% 
  as.data.frame()
ProjChanges_ssp126b$scenario <- scenario126

scenario245 <- rep("TSC 2-4.5", 323140)
ProjChanges_ssp245 <- raster("Future_WC2.1_Soil/2_Projection_changes/Changes_EC/Period_1/Scenario_ssp245/binary_comparison.tif") 
ProjChanges_ssp245b <- ProjChanges_ssp245 %>% as("SpatialPixelsDataFrame") %>% 
  as.data.frame()
ProjChanges_ssp245b$scenario <- scenario245  

scenario370 <- rep("TSC 3-7.0", 323140)
ProjChanges_ssp370 <- raster("Future_WC2.1_Soil/2_Projection_changes/Changes_EC/Period_1/Scenario_ssp370/binary_comparison.tif") 
ProjChanges_ssp370b <- ProjChanges_ssp370 %>% as("SpatialPixelsDataFrame") %>% 
  as.data.frame()
ProjChanges_ssp370b$scenario <- scenario370

scenario585 <- rep("TSC 5-8.5", 323140)
ProjChanges_ssp585 <- raster("Future_WC2.1_Soil/2_Projection_changes/Changes_EC/Period_1/Scenario_ssp585/binary_comparison.tif") 
ProjChanges_ssp585b <- ProjChanges_ssp585 %>% as("SpatialPixelsDataFrame") %>% 
  as.data.frame()
ProjChanges_ssp585b$scenario <- scenario585

ProjChangesA <-rbind(ProjChanges_ssp126b, ProjChanges_ssp245b)
ProjChangesB <-rbind(ProjChanges_ssp370b, ProjChanges_ssp585b)
ProjChanges <-rbind(ProjChangesA, ProjChangesB)

ProjChanges <- ProjChanges %>%
  mutate(binary_comparison = case_when(binary_comparison == 1 ~ "Area ganada en 1 MCG",
                            binary_comparison == 2 ~ "Area ganada en 2 MCGs",
                            binary_comparison == 3 ~ "Area ganada en 3 MCGs",
                            binary_comparison == 4 ~ "Area ganada en 4 MCGs",
                            binary_comparison == 5 ~ "Area ganada en 5 MCGs",
                            binary_comparison == 6 ~ "Area perdida en 5 MCGs",
                            binary_comparison == 7 ~ "Area perdida en 4 MCGs",
                            binary_comparison == 8 ~ "Area perdida en 3 MCGs",
                            binary_comparison == 9 ~ "Area perdida en 2 MCGs",
                            binary_comparison == 10 ~ "Area perdida en 1 MCG",
                            binary_comparison == 11 ~ "Estable")) %>% 
  dplyr::filter(!is.na(binary_comparison)) %>% 
  mutate(binary_comparison = fct_relevel(binary_comparison, "Area ganada en 5 MCGs", "Area ganada en 4 MCGs", "Area ganada en 3 MCGs", "Area ganada en 2 MCGs", "Area ganada en 1 MCG",
                              "Estable", "Area perdida en 1 MCG", "Area perdida en 2 MCGs", "Area perdida en 3 MCGs", "Area perdida en 4 MCGs", "Area perdida en 5 MCGs"))

P <- ggplot() +  
  geom_sf(data = Sudam, size = 0.1) +
  geom_tile(data = ProjChanges, aes(x = x, y = y, fill = binary_comparison)) + 
  geom_sf(data = Sudam, fill = NA) +
  facet_wrap(~scenario,nrow = 1, ncol = 4) +
  scale_fill_manual(name = NULL, values = c('#bd0026','#f03b20','#fd8d3c','#fecc5c','#ffffb2',
                                            '#00441b','#3182bd','#92c5de','#4393c3','#2166ac','#053061')) + 
  ylab("") + xlab("") +
  coord_sf(xlim = sf::st_bbox(AreaM_SF)[c(1,3)],
           ylim = sf::st_bbox(AreaM_SF)[c(2,4)],
           expand = FALSE) +
  theme_no_axes() +
  theme(plot.background = element_rect(fill = "NA"),
        plot.margin = margin(t = -0.6, r = 0.1, b = -0.6, l = 0, unit = "cm"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "grey50", size = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.25, 'cm'),
        legend.margin=margin(t = -0.3, r = 0, b = 0, l = 0, unit = "cm")) +
  guides(fill = guide_legend(ncol = , nrow = 2))
P