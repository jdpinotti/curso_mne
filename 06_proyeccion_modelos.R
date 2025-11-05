##############################################
# Proyeccion en el Espacio Geo
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################
library(terra)
library(ENMeval)
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")

#### invocar el modelo####
cvenustus_model <- readRDS("cvenustus_enmeval.rds")
show(cvenustus_model)
#### invocar el stack de capas del año a proyectar ###
ruta <- "capas_chelsa/"

envs <- rast(list.files(ruta, pattern = "\\.tif$", full.names = TRUE))
envs
names(envs)

# Crea un vector con los nuevos nombres (en el mismo orden!)
nuevos_nombres <- c(
  "bio01", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17",
  "bio18", "bio19", "bio02", "bio03", "bio04", "bio05", "bio06", "bio07", "bio08", "bio09")

# Asigna los nuevos nombres
names(envs) <- nuevos_nombres
names(envs)

### Proyectar ######

suitability_map <- predict(envs, cvenustus_model@models[[3]]) ### Poner modelo selccionado aca!

writeRaster(suitability_map, filename = "outputs/modelo_cvenustus", filetype = "GTiff", overwrite = TRUE)

