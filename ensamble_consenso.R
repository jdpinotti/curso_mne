library(terra)

# 1. Cargar los modelos
rutas <- list.files("mapas_exportados/para_consenso/", pattern = "\\.tif$", full.names = TRUE)
modelos <- rast(rutas)

# Ver nombres
names(modelos)


umbral <- 0.025
modelos_bin <- modelos > umbral

#Para umbrales distintos para cada modelo
#umbrales <- c(0.45, 0.51, 0.47, 0.53) # uno por cada modelo
#modelos_bin <- modelos
#for (i in 1:nlyr(modelos)) {
#  modelos_bin[[i]] <- modelos[[i]] > umbrales[i]
#}

ensamble <- sum(modelos_bin)
writeRaster(ensamble, "ensamble_suma.tif", overwrite = TRUE)

