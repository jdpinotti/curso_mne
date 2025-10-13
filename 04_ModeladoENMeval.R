##############################################
# Modelado con ENMeval 2.0
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################

# Cargar las librerías necesarias
library(terra)
library(dplyr)
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")
getwd()

# Definir la ruta de los archivos
csv_path <- "cvenusus_clean_thin1.csv"

datos <- read.csv(csv_path)

tif_path <- list.files("/capas_chelsa", pattern = "*.tif", full.names = TRUE)


#################Generacion bg ##############################
library(sp)
library(dismo)
library(sf)
mask <- st_read("/home/dargwind/Modelado de Nicho/Tmatacus/Area_M.gpkg")

mask_sp <- as_Spatial(mask)

library(terra)

# Ruta al directorio donde están los rasters
ruta <- "/home/dargwind/Modelado de Nicho/Tmatacus/modelo_2km/promedios_ponderados"

# Cargar todos los .tif del directorio
envs <- rast(list.files(ruta, pattern = "\\.tif$", full.names = TRUE))
envs
names(envs)
# Generar 10000 puntos aleatorios dentro de la máscara
set.seed(42)  # Fijar la semilla para reproducibilidad
bg_points <- st_sample(mask, size = 10000, type = "random")

# Extraer valores ambientales en los puntos de background
bg_values <- terra::extract(envs, vect(bg_points))

head(bg_values)
dim(bg_values)

# Crear un DataFrame para el background y guardar en un CSV
bg_full <- data.frame(st_coordinates(bg_points), bg_values) %>% 
  dplyr::select(!3)

write.csv(bg_full, "background_tmatacus.csv", row.names = FALSE)


###### Modelo con las variables selecionadas por PCA ######

library(dplyr)

occs <- read.csv("tmatacus_val_ambientales.csv")
names(occs)
occs <- dplyr::select(occs, 
                      longitude, 
                      latitude,
                      bio01 = X2km_bios_bio01,
                      bio04 = X2km_bios_bio04,
                      bio14 = X2km_bios_bio14,
                      bio18 = X2km_bios_bio18,
                      close_woodland = X2km_porc_close_woodland,
                      crop = X2km_porc_crop)



names(occs)
names(bg_full)
bg <- dplyr::select(bg_full,
                    longitude = X,
                    latitude = Y,
                    bio01 = promedio_ponderado_2km_bios_bio01,
                    bio04 = promedio_ponderado_2km_bios_bio04,
                    bio14 = promedio_ponderado_2km_bios_bio14,
                    bio18 = promedio_ponderado_2km_bios_bio18,
                    close_woodland = promedio_ponderado_2km_porc_close_woodland,
                    crop = promedio_ponderado_2km_porc_crop)

names(bg)
names(occs)
###################### MODELADO #########################

#Agregar pROC

library(kuenm)
proc <- function(vars) {
  proc <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
  out <- data.frame(proc_auc_ratio = proc$pROC_summary[1], 
                    proc_pval = proc$pROC_summary[2], row.names = NULL)
  return(out)
}

library(ENMeval)

# Dividir las ocurrencias en entrenamiento y prueba
#set.seed(42)
#train_fraction <- 0.8
#n_occs <- nrow(occs)
#train_indices <- sample(1:n_occs, size = floor(train_fraction * n_occs))
#occs_training <- occs[train_indices, ]
#occs_testing <- occs[-train_indices, ]


#options(java.parameters = "-Xmx2g") # Aumentar la memoria para MaxEnt si es necesario
#library(dismo)
#system.file("maxent.jar", package = "dismo")

#fc = c("L", "Q","LQHPT", "LQ", "LQH", "LQHT", "LQHP", "LQP")

tmatacus_model <- ENMevaluate(
  occs = occs,
  envs = NULL,            # No se usan capas en RasterStack para el enfoque SWD
  bg = bg,                # Background con los mismos valores ambientales
  algorithm = "maxent.jar", 
  partitions = 'block',
  user.eval = proc,
  doClamp = TRUE,
  tune.args = list(fc = c("L", "Q", "LQ", "LQH", "LQHT"), rm = c(0.5, 1, 1.5, 2, 5)),
  taxon.name = "Tolypeutes matacus"
)

#c("L", "Q","LQHPT", "LQ", "LQH", "LQHT", "LQHP", "LQP")
#para ir viendo resultados

show(tmatacus_model)

eval.results(tmatacus_model)
write.csv(tmatacus_model@results, "modelos_tmatacus_PCA.csv")



####model 8 was selected (FC lqp, rm, o.5)####

mod.seq <- eval.models(tmatacus_model)[[8]]
mod.seq@results
pr <- predicts::partialResponse(tmatacus_model@models[[8]])

###########guardar el modelo #########

saveRDS(tmatacus_model, file = "tmatacus_enmeval_pca.rds")

