##############################################
# Alinear una capa raster a otra de referencia
##############################################

# ---- 1. Cargar librerías ----
library(terra)

# ---- 2. Definir rutas ----
# Capa a ajustar (por ejemplo, suelo o NDVI)
ruta_suelo <- "suelo.tif"

# Capa de referencia (por ejemplo, bio1)
ruta_bio <- "bio1.tif"

# Salida
ruta_salida <- "suelo_alineado.tif"

# ---- 3. Cargar rasters ----
suelo <- rast(ruta_suelo)
bio <- rast(ruta_bio)

# ---- 4. Reproyectar si el CRS es diferente ----

suelo <- project(suelo, bio, method = "bilinear")


# ---- 5. Recortar al área de la capa bioclimática ----

suelo_crop <- crop(suelo, bio)

# ---- 6. Remuestrear para igualar resolución y alineamiento ----
cat("Alineando resolución y grilla...\n")
suelo_aligned <- resample(suelo_crop, bio, method = "bilinear")

# ---- 7. Verificar alineamiento ----
cat("Verificando geometría...\n")
print(compareGeom(suelo_aligned, bio, stopOnError = FALSE))

# ---- 8. Guardar resultado ----
cat("Guardando raster alineado...\n")
writeRaster(suelo_aligned, ruta_salida, overwrite = TRUE)

cat("\n✅ Proceso completado. Raster guardado en:", ruta_salida, "\n")
