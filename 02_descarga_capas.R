##############################################
#  Descarga de capas climáticas CHELSA (Sudamérica)
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################

## 1. Cargar librerías necesarias ----
library(terra)
library(utils)

## 2. Configuración inicial ----
# Definir directorio de trabajo
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")

# Carpeta para guardar capas
dir_capas <- "capas_chelsa"
if (!dir.exists(dir_capas)) dir.create(dir_capas)

## 3. URLs de las capas en tu repositorio GitHub ----
# Base URL de tu repositorio
base_url <- "https://github.com/jdpinotti/curso_mne/raw/main/capas_chelsa/"

# Lista de todas las variables bioclimáticas (bio1 a bio19)
variables <- sprintf("bio%d.tif", 1:19)

# URLs completas
urls <- paste0(base_url, variables)

# Nombres de archivos locales
archivos <- file.path(dir_capas, variables)

## 4. Descargar capas si no existen ----
for (i in seq_along(urls)) {
  if (!file.exists(archivos[i])) {
    cat("Descargando:", variables[i], "\n")
    tryCatch({
      download.file(urls[i], archivos[i], mode = "wb")
      cat("✓ Descarga completada:", variables[i], "\n")
    }, error = function(e) {
      cat("✗ Error descargando:", variables[i], "-", e$message, "\n")
    })
  } else {
    cat("Archivo ya existe:", variables[i], "\n")
  }
}

## 5. Verificar descargas ----
archivos_descargados <- list.files(dir_capas, pattern = "\\.tif$", full.names = TRUE)
cat("Se descargaron", length(archivos_descargados), "archivos de", length(urls), "\n")