##############################################
#  Analisis de Nicho con tipicidad de mahalanobis
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################

# ============================================================================
# ANÁLISIS DE NICHO ECOLÓGICO USANDO TIPICIDAD DE MAHALANOBIS
# ============================================================================

# ---- 1. Librerías ----
library(terra)
library(ntbox)
library(dplyr)

# ---- 2. Configuración del entorno ----
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")
rm(list = ls())
gc()

# ---- 3. Cargar datos ----
# Presencias
pres <- read.csv("cvenustus_clean_thin1.csv")
colnames(pres) <- c("species", "lon", "lat")

# Capas ambientales
env_files <- list.files(
  ruta_capas,
  pattern = "CHELSA_bio(6|17|4|7).*\\.tif$",
  full.names = TRUE
)

env_files
# Cargar como raster stack de terra
env <- rast(env_files)
names(env) <- sub(".*bio(\\d+).*", "bio\\1", names(env))

cat("Variables cargadas:", names(env), "\n")
cat("Número de presencias:", nrow(pres), "\n")

# ---- 4. Extraer valores ambientales de las presencias ----
pres_env <- terra::extract(env, pres[, c("lon", "lat")], ID = FALSE)
pres_full <- cbind(pres, pres_env)
pres_full <- na.omit(pres_full)

cat("Presencias con datos ambientales:", nrow(pres_full), "\n")

# ---- 5. Calcular centroide y matriz de covarianza ----
# Variables ambientales
vars_amb <- names(env)
X <- as.matrix(pres_full[, vars_amb])

# Centroide (media multivariada)
centroid <- colMeans(X)

# Matriz de covarianza
cov_matrix <- cov(X)

cat("\n=== CENTROIDE DEL NICHO ===\n")
print(centroid)

cat("\n=== VOLUMEN DEL NICHO ===\n")
# Volumen del elipsoide (proporcional al determinante)
niche_volume <- det(cov_matrix)^0.5
cat("Volumen:", niche_volume, "\n")

# ---- 6. Calcular distancias de Mahalanobis para las presencias ----
dists_pres <- mahalanobis(X, center = centroid, cov = cov_matrix)

# Tipicidad (probabilidad chi-cuadrado)
# Valores altos = más típico
typicality_pres <- pchisq(dists_pres, df = ncol(X), lower.tail = FALSE)

# Añadir al dataframe
pres_full$mahal_dist <- dists_pres
pres_full$typicality <- typicality_pres

cat("\n=== ESTADÍSTICAS DE TIPICIDAD ===\n")
cat("Media:", mean(typicality_pres), "\n")
cat("Min:", min(typicality_pres), "\n")
cat("Max:", max(typicality_pres), "\n")

# ---- 7. Proyectar el modelo sobre las capas ambientales ----
cat("\n=== PROYECTANDO MODELO ===\n")

# Función para calcular tipicidad pixel por pixel
calc_typicality <- function(env_values, centroid, cov_matrix) {
  # Remover NAs
  na_idx <- complete.cases(env_values)
  
  if (sum(na_idx) == 0) {
    return(rep(NA, nrow(env_values)))
  }
  
  # Calcular distancias de Mahalanobis
  dists <- mahalanobis(env_values, center = centroid, cov = cov_matrix)
  
  # Tipicidad
  typicality <- pchisq(dists, df = length(centroid), lower.tail = FALSE)
  
  return(typicality)
}

# Aplicar sobre el raster
typicality_raster <- app(env, function(x) {
  calc_typicality(matrix(x, ncol = nlyr(env)), centroid, cov_matrix)
})

names(typicality_raster) <- "typicality"

# ---- 9. Exportar resultados ----
# Crear carpeta de salida
dir.create("resultados_mahalanobis", showWarnings = FALSE)

# Guardar presencias con tipicidad
write.csv(pres_full, "resultados_mahalanobis/presencias_con_tipicidad.csv", 
          row.names = FALSE)

# Guardar rasters
writeRaster(typicality_raster, 
            "resultados_mahalanobis/typicality_map.tif", 
            overwrite = TRUE)

writeRaster(suitability_raster, 
            "resultados_mahalanobis/suitability_map.tif", 
            overwrite = TRUE)

# Guardar parámetros del modelo
sink("resultados_mahalanobis/parametros_modelo.txt")
cat("=== PARÁMETROS DEL MODELO DE NICHO ===\n\n")
cat("Especies:", unique(pres$species), "\n")
cat("Número de presencias:", nrow(pres_full), "\n")
cat("Variables ambientales:", paste(vars_amb, collapse = ", "), "\n\n")
cat("CENTROIDE:\n")
print(centroid)
cat("\nMATRIZ DE COVARIANZA:\n")
print(cov_matrix)
cat("\nVOLUMEN DEL NICHO:", niche_volume, "\n")
sink()

cat("\n=== RESULTADOS GUARDADOS ===\n")
cat("Carpeta: resultados_mahalanobis/\n")
cat("- presencias_con_tipicidad.csv\n")
cat("- typicality_map.tif\n")
cat("- suitability_map.tif\n")
cat("- parametros_modelo.txt\n")

# ---- 10. Visualización básica ----
# Plot de tipicidad
plot(typicality_raster, main = "Tipicidad de Mahalanobis")
points(pres_full[, c("lon", "lat")], pch = 20, cex = 0.5, col = "red")

# Plot de idoneidad
plot(suitability_raster, main = "Idoneidad (Suitability)")
points(pres_full[, c("lon", "lat")], pch = 20, cex = 0.5, col = "red")

# Histograma de tipicidad en presencias
hist(pres_full$typicality, 
     main = "Distribución de Tipicidad en Presencias",
     xlab = "Tipicidad", 
     col = "steelblue")
