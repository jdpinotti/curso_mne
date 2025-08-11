##############################################
#  Limpieza y filtrado de ocurrencias
#  Curso: Introducción a los Modelos Predictivos de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################

## 1. Cargar librerías necesarias ----
library(CoordinateCleaner)  # Limpieza de coordenadas
library(spThin)             # Filtrado espacial (thinning)
library(dplyr)              # Manipulación de datos
library(sf)                 # Datos espaciales
library(mapview)            # Mapas interactivos

## 2. Configuración inicial ----
# Definir el directorio de trabajo (ajustar si es necesario)
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")

# Nombre del archivo de ocurrencias
csv_file <- "cvenustus_compilado.csv"

# Nombres de columnas
lat_col <- "lat"
long_col <- "long"
sp_col <- "sp"

## 3. Leer el archivo de ocurrencias ----
occs_crudas <- read.csv(csv_file)

# Mostrar información inicial
cat("Número de registros originales:", nrow(occs_crudas), "\n")
head(occs_crudas)

## 4. Limpieza de coordenadas ----
# Detectar registros problemáticos: duplicados, coordenadas en el mar, ceros, etc.
occs_clean_flags <- clean_coordinates(
  x = occs_crudas,
  lon = long_col,
  lat = lat_col,
  species = sp_col,
  tests = c("centroids", "equal", "duplicates", "gbif", "seas", "zeros")
)

# Filtrar solo los registros válidos
occs_clean <- occs_crudas[occs_clean_flags$.summary, ]

cat("Registros después de limpieza:", nrow(occs_clean), "\n")

## 5. Visualizar registros limpios

occs_clean_sf <- st_as_sf(occs_clean, coords = c(long_col, lat_col), crs = 4326)
mapview(occs_clean_sf)


## 6. Filtrado espacial (thinning) ----
# Queremos un mínimo de 10 km entre registros
dist_min <- 10  # en km

filtered_list <- thin(
  loc.data = occs_clean,
  lat.col = lat_col,
  long.col = long_col,
  spec.col = sp_col,
  thin.par = dist_min,
  reps = 1,
  locs.thinned.list.return = TRUE,
  write.files = TRUE,
  out.dir = ".",
  out.base = "cvenusus_clean",
  verbose = TRUE
)

# Extraer el data.frame con las ocurrencias filtradas
filtered_data <- filtered_list[[1]]
cat("Número de registros después de thinning:", nrow(filtered_data), "\n")
