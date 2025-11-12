library(terra)

# Ruta del archivo multibanda
tif_path <- "/home/dargwind/Modelado de Nicho/curso_modelado/wd/mapas_exportados/predicciones_individuales_continuas.tif"

# Cargar el archivo multibanda
multibanda <- rast(tif_path)

# Ruta del directorio de salida
output_dir <- "/home/dargwind/Modelado de Nicho/curso_modelado/wd/mapas_exportados/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Obtener nombres de las bandas
nombres_bandas <- names(multibanda)

# Crear un bucle para exportar cada banda con su nombre original
for (i in 1:nlyr(multibanda)) {
  # Extraer la banda
  banda <- multibanda[[i]]
  
  # Usar el nombre original de la banda
  nombre_banda <- paste0(output_dir, nombres_bandas[i], ".tif")
  
  # Exportar la banda como un archivo raster
  writeRaster(banda, filename = nombre_banda, overwrite = TRUE)
  
  # Mensaje para confirmar la exportación
  cat("Exportado:", nombre_banda, "\n")
}

