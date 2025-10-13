##############################################
# Modelado con ENMeval 2.0
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################

# Cargar las librerías necesarias
library(terra)
library(dplyr)
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")

# Definir la ruta de los archivos
csv_path <- "cvenusus_clean_thin1.csv"

datos <- read.csv(csv_path)

tif_path <- list.files("/corte_m", pattern = "*.tif", full.names = TRUE)


#################Generacion bg ##############################
library(sp)
library(dismo)
library(sf)
mask <- st_read("mvenustus.shp")

mask_sp <- as_Spatial(mask)

# Ruta al directorio donde están los rasters

ruta <- "corte_m"

# Cargar todos los .tif del directorio
envs <- rast(list.files(ruta, pattern = "\\.tif$", full.names = TRUE))
envs
names(envs)

# Función para simplificar nombres
renombrar_bio <- function(nombres_originales) {
  # Extraer el número bio usando expresiones regulares
  numeros_bio <- gsub(".*bio(\\d+).*", "\\1", nombres_originales)
  # Crear nuevos nombres
  paste0("bio", numeros_bio)
}

# Aplicar la función
names(envs) <- renombrar_bio(names(envs))
print(names(envs))

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

write.csv(bg_full, "background_cvenustus.csv", row.names = FALSE)


###### Modelo con las variables selecionadas por PCA ######

occs <- read.csv("var_amb_presencias.csv")
names(occs)
# Simplificar nombres
nombres_finales <- c("sp", "long", "lat", 
                     "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9",
                     "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", 
                     "bio18", "bio19")

names(occs) <- nombres_finales
print(names(occs))

#Seleccionar las variables de acuerdo al PCA

occs <- dplyr::select(occs, 
                      long, 
                      lat,
                      bio03 = bio3,
                      bio04 = bio4,
                      bio06 = bio6,
                      bio07 = bio7,
                      bio17 = bio17)

names(occs)
names(bg_full)
bg <- dplyr::select(bg_full,
                    long = X,
                    lat = Y,
                    bio03 = bio3,
                    bio04 = bio4,
                    bio06 = bio6,
                    bio07 = bio7,
                    bio17 = bio17)

names(bg)
names(occs)
###################### MODELADO #########################

library(ENMeval)


cvenustus_model <- ENMevaluate(
  occs = occs,
  envs = NULL,            # No se usan capas en RasterStack para el enfoque SWD
  bg = bg,                # Background con los mismos valores ambientales
  algorithm = "maxent.jar", 
  partitions = 'block',
  user.eval = NULL,
  doClamp = TRUE,
  tune.args = list(fc = c("L", "Q", "LQ", "LQH", "LQHT"), rm = c(0.5, 1, 1.5, 2, 5)),
  taxon.name = "Calomys venustus"
)

#c("L", "Q","LQHPT", "LQ", "LQH", "LQHT", "LQHP", "LQP")

#para ir viendo resultados
show(cvenustus_model)
eval.results(cvenustus_model)
write.csv(cvenustus_model@results, "modelos_cvenustus_enmeval.csv")


####Elegir el mejor modelo ####


mod.seq <- eval.models(cvenustus_model)[[3]]
mod.seq@results
write.csv(mod.seq@results, "resultados_modelo_3.csv")


pr <- predicts::partialResponse(cvenustus_model@models[[3]])

###########guardar el modelo #########

saveRDS(cvenustus_model, file = "cvenustus_enmeval.rds")

