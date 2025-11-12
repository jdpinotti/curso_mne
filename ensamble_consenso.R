##############################################
#  Script: 08_ensamble_simple_desde_rds.r
#  Curso: Modelado de Nicho y Distribución
##############################################

# ---- 1. Librerías ----
library(terra)
library(dplyr)
library(biomod2)

# ---- 2. Configuración ----
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")

# Carpeta de salida
dir.create("mapas_exportados/consenso", recursive = TRUE, showWarnings = FALSE)

# ---- 3. Cargar los objetos necesarios ----
# Asegúrate de tener estos archivos generados por scripts anteriores:
# - modelos_rds/biomod_modelos.rds (modelos ajustados)
# - modelos_rds/biomod_proyeccion_current.rds (proyecciones)
# Si usás otra especie, cambiale el nombre de carpeta o archivo.

myBiomodModelOut <- readRDS("modelos_rds/biomod_modelos.rds")
myBiomodProjection <- readRDS("modelos_rds/biomod_proyeccion_current.rds")

cat("✓ Objetos cargados correctamente.\n")

# ---- 4. Obtener predicciones continuas ----
predicciones <- get_predictions(myBiomodProjection)
names(predicciones)

# ---- 5. Seleccionar algoritmos de interés ----
# Cambiá esta línea según tus modelos disponibles:
# Ejemplo: c("RF", "GLM"), o c("RF", "GAM", "GBM"), etc.
algoritmos_usar <- c("RF", "GLM")

# Filtrar solo los modelos de esos algoritmos
modelos_sel <- predicciones[[grep(paste(algoritmos_usar, collapse = "|"), names(predicciones))]]

cat("Modelos seleccionados:\n")
print(names(modelos_sel))

# ---- 6. Obtener los umbrales TSS desde las evaluaciones ----
evals <- get_evaluations(myBiomodModelOut)

# Convertir la lista compleja a un data.frame utilizable
umbrales_TSS <- as.data.frame.table(evals["TSS", "Cutoff", , , , drop = FALSE]) %>%
  rename(algo = Var4, umbral = Freq) %>%
  group_by(algo) %>%
  summarise(umbral = mean(umbral, na.rm = TRUE))

cat("\nUmbrales promedio TSS por algoritmo:\n")
print(umbrales_TSS)

# ---- 7. Binarizar los modelos seleccionados ----
modelos_bin <- list()

for (i in 1:nlyr(modelos_sel)) {
  modelo_nombre <- names(modelos_sel)[i]
  
  # Buscar qué algoritmo corresponde (ej. RF, GLM, etc.)
  algo_encontrado <- NA
  for (a in algoritmos_usar) {
    if (grepl(a, modelo_nombre)) {
      algo_encontrado <- a
      break
    }
  }
  
  if (is.na(algo_encontrado)) next
  
  umbral <- umbrales_TSS$umbral[umbrales_TSS$algo == algo_encontrado]
  
  # Binarizar: 1 = hábitat adecuado, 0 = no adecuado
  bin <- modelos_sel[[i]] >= umbral
  modelos_bin[[modelo_nombre]] <- bin
  
  # Exportar binario individual
  writeRaster(bin,
              filename = paste0("mapas_exportados/consenso/", modelo_nombre, "_TSSbin.tif"),
              overwrite = TRUE)
  
  cat("✅ Binarizado y exportado:", modelo_nombre, " | Umbral =", round(umbral, 3), "\n")
}

# ---- 8. Crear el ensamble (suma de binarios) ----
ensamble <- sum(rast(modelos_bin), na.rm = TRUE)

# Exportar el ensamble final
writeRaster(ensamble, "mapas_exportados/consenso/ENSAMBLE_suma_RF_GLM.tif", overwrite = TRUE)
cat("\n🎯 Ensamble final exportado: 'ENSAMBLE_suma_RF_GLM.tif'\n")

# ---- 9. Verificar el rango de valores ----
cat("\nRango del ensamble (esperado: 0 a", length(modelos_bin), "):\n")
print(range(values(ensamble), na.rm = TRUE))

# ---- 10. Visualizar ----
plot(ensamble, main = "Ensamble simple (suma binaria)")


#Despues cargar los mapas binarios en QGIS y sumarlos con calculadora raster