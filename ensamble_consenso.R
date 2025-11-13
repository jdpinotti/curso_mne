##############################################
#  Script: 08_ensamble_simple_desde_tif_y_rds.r
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################

# ---- 1. Librerías ----
library(terra)
library(dplyr)

# ---- 2. Configuración ----
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")

# Directorios
dir_in <- "mapas_exportados/para_consenso/" #colocar los tif de los modelos que quiero ensamblar en esta carpeta
dir_out <- "mapas_exportados/consenso/"
dir_rds <- "modelos_rds/"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# ---- 3. Cargar archivos ----
# Rasters continuos (uno por modelo)
rutas <- list.files(dir_in, pattern = "\\.tif$", full.names = TRUE)
modelos <- rast(rutas)
cat("✓ Modelos cargados:", nlyr(modelos), "bandas\n")

# Evaluaciones biomod (umbrales TSS)
evals <- readRDS(file.path(dir_rds, "evaluaciones_modelos.rds"))

# ---- 4. Seleccionar algoritmos de interés ----
# Cambiá esta línea según tus modelos disponibles
algoritmos_usar <- c("RF", "GLM", "GAM", "GBM")

# ---- 5. Extraer umbrales TSS por algoritmo ----
umbrales_TSS <- evals %>%
  filter(metric.eval == "TSS") %>%
  select(algo, cutoff) %>%
  group_by(algo) %>%
  summarise(umbral = mean(cutoff, na.rm = TRUE))

cat("\nUmbrales promedio TSS:\n")
print(umbrales_TSS)

# ---- 6. Binarizar modelos ----
modelos_bin <- list()

for (i in 1:nlyr(modelos)) {
  modelo_nombre <- names(modelos)[i]
  
  # Detectar algoritmo dentro del nombre del raster
  algo_encontrado <- NA
  for (a in algoritmos_usar) {
    if (grepl(a, modelo_nombre, ignore.case = TRUE)) {
      algo_encontrado <- a
      break
    }
  }
  
  if (is.na(algo_encontrado)) next  # omitir si no pertenece
  
  #Si quiero puedo hacer un vector con los umbrales que elija (por ej. umbral <- c(0.05, 0.25, 0.32, 0.25))
  umbral <- umbrales_TSS$umbral[umbrales_TSS$algo == algo_encontrado]
  
  # Binarizar
  bin <- modelos[[i]] >= umbral
  modelos_bin[[modelo_nombre]] <- bin
  
  # Exportar binario
  writeRaster(bin,
              filename = paste0(dir_out, modelo_nombre, "_TSSbin.tif"),
              overwrite = TRUE)
  
  cat("✅", modelo_nombre, "binarizado (", algo_encontrado, ", umbral =", round(umbral, 3), ")\n")
}

# ---- 7. Crear ensamble ----
ensamble <- sum(rast(modelos_bin), na.rm = TRUE)

# Exportar ensamble final
nombre_ensamble <- paste0("ENSAMBLE_suma_", paste(algoritmos_usar, collapse = "_"), ".tif")
writeRaster(ensamble, paste0(dir_out, nombre_ensamble), overwrite = TRUE)

cat("\n🎯 Ensamble final exportado:", nombre_ensamble, "\n")