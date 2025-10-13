##############################################
# Corte por Area de Acceso Histórico (M)
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################

# ---- Librerías ----
library(terra)
library(sf)

# ---- Configuración ----
# 
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")

# Carpeta con las capas CHELSA
dir_capas <- "capas_chelsa"

# Carpeta de salida
dir_salida <- "corte_m"
if (!dir.exists(dir_salida)) dir.create(dir_salida)

# ---- 1. Cargar shape ----
# Descargar shapefile desde GitHub si no existe localmente
url_shape <- "https://github.com/jdpinotti/curso_mne/raw/main/mvenustus/mvenustus.shp"

# Comprobamos si ya existe el shape localmente
shape_path <- "mvenustus.shp"
if (!file.exists(shape_path)) {
  message("Descargando shapefile desde GitHub...")
  download.file(url_shape, destfile = shape_path, mode = "wb")
  # También descargamos los archivos complementarios (.dbf, .shx, .prj)
  archivos <- c("mvenustus.dbf", "mvenustus.shx", "mvenustus.prj")
  for (a in archivos) {
    url_a <- paste0("https://github.com/jdpinotti/curso_mne/raw/main/mvenustus/", a)
    download.file(url_a, destfile = a, mode = "wb")
  }
}

# Cargar el shape
m_shape <- vect("mvenustus.shp")

# ---- 2. Cargar capas bioclimáticas ----
# Busca todos los archivos .tif en la carpeta capas_chelsa
capas_files <- list.files(dir_capas, pattern = "\\.tif$", full.names = TRUE)

# Cargar como SpatRaster
capas <- rast(capas_files)

# ---- 3. Recorte ----
capas_corte <- crop(capas, m_shape)
capas_m <- mask(capas_corte, m_shape)

# ---- 4. Guardar resultados ----

for (i in 1:nlyr(capas_m)) {
  nombre <- names(capas_m)[i]
  writeRaster(capas_m[[i]], 
              filename = file.path(dir_salida, paste0(nombre, ".tif")), 
              overwrite = TRUE)
}
