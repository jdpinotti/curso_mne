library(terra)

# 1. Cargar los modelos
rutas <- list.files("mapas_exportados/para_consenso/", pattern = "\\.tif$", full.names = TRUE)
modelos <- rast(rutas)

# 2. Revisar nombres
names(modelos)

# 3. Binarizar (mismo umbral para todos)
umbral <- 0.25
modelos_bin <- modelos > umbral

# 4. Sumar los modelos binarios ignorando NA
ensamble <- sum(modelos_bin, na.rm = TRUE)

# 5. Exportar resultado
writeRaster(ensamble, "ensamble_suma.tif", overwrite = TRUE)

# 6. Ver rango esperado
range(values(ensamble), na.rm = TRUE)
