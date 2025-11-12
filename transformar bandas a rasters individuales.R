library(terra)

# Ruta del archivo multibanda
tif_path <- "/media/dargwind/Elements/CHELSA_ts/bioclima/bios_2018.tif"

# Cargar el archivo multibanda
multibanda <- rast(tif_path)

# Ruta del directorio de salida (asegúrate de que exista y termina con '/')
output_dir <- "/media/dargwind/Elements/Tmatacus_tenm/2018/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Obtener el número de bandas
num_bandas <- nlyr(multibanda)

# Crear un bucle para exportar cada banda
for (i in 1:num_bandas) {
  # Extraer la banda
  banda <- multibanda[[i]]
  
  # Crear el nombre del archivo (bio01.tif, bio02.tif, etc.)
  nombre_banda <- paste0(output_dir, "bio", sprintf("%02d", i), ".tif")
  
  # Exportar la banda como un archivo raster
  writeRaster(banda, filename = nombre_banda, overwrite = TRUE)
  
  # Mensaje para confirmar la exportación
  cat("Exportado:", nombre_banda, "\n")
}
