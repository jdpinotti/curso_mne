##############################################
#  Script: 07_introduccion_a_biomod2.r
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
#  Objetivo:
# ###################################

# ---- 1. Librerías ----
library(biomod2)
library(terra)

# ---- 2. Configuración del entorno ----
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")
rm(list = ls())
gc()

# ---- 3. Cargar datos ----
# Presencias (archivo CSV con columnas: species, lon, lat)
pres <- read.csv("cvenustus_clean_thin1.csv")
colnames(pres) <- c("species", "lon", "lat")

# ---- 3bis. Selección de capas ambientales ----

# RUTA general de capas
ruta_capas <- "capas_chelsa"

# OPCIÓN A: cargar todas las capas
# env_files <- list.files(ruta_capas, pattern = ".tif$", full.names = TRUE)

# OPCIÓN B: seleccionar solo algunas variables o subcarpetas

env_files <- list.files(
  ruta_capas,
  pattern = "bio(6|17|4|7)\\.tif$", #aca poner los numeros de las bios seleccionadas
  full.names = TRUE
)
env_files

# Cargar como raster stack
env <- rast(env_files)
names(env)
# Renombrar para quedarse solo con "bioX"
names(env) <- sub(".*bio(\\d+).*", "bio\\1", names(env))
# Verificar
names(env)

# Visualizar una capa y las presencias
plot(env[[1]], main = paste("Ejemplo de variable:", names(env)[1]))
points(pres$lon, pres$lat, pch = 20, col = "red")

# ---- 4. Preparar datos para biomod2 ----
myBiomodData <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(pres)),          # 1 = presencia
  expl.var = env,                         # variables ambientales seleccionadas
  resp.xy = pres[, c("lon", "lat")],      # coordenadas
  resp.name = "Calomys_venustus",         # nombre de la especie
  PA.nb.rep = 1,                          # nº de repeticiones de pseudoausencias
  PA.nb.absences = 10000,                 # nº de pseudoausencias
  PA.strategy = "random"                  # estrategia de selección
)
myBiomodData

# ---- 5. Ajuste de modelos ----
# Configurar opciones

myBiomodOptions <- biomod2::bm_ModelingOptions(
  data.type = "binary",
  strategy = "default"  # o "bigboss", "user.defined"
)

myBiomodModelOut <- biomod2::BIOMOD_Modeling(
  bm.format = myBiomodData,           # Datos formateados
  modeling.id = "intro_biomod2",      # ID del proyecto
  models = c("RF", "GBM", "GLM", "GAM"),  # Modelos a usar
  bm.options = myBiomodOptions,       # Opciones configuradas
  nb.rep = 3,                         # Repeticiones
  data.split.perc = 70,               # % entrenamiento
  var.import = 3,                     # Importancia variables
  metric.eval = c("TSS", "ROC", "KAPPA"),  # Métricas
  do.progress = TRUE,                 # Mostrar progreso
  scale.models = TRUE,                # Escalar predicciones
  nb.cpu = 1                          # Núcleos CPU (opcional)
)

# ---- 6. Evaluación ----

# 1. Gráfico principal de evaluación
biomod2::bm_PlotEvalMean(bm.out = myBiomodModelOut)

# 2. Boxplot con DOS grupos (correcto)
biomod2::bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, 
                            group.by = c("algo", "full.name"))

# 3. Ver tabla de resultados
evals <- biomod2::get_evaluations(myBiomodModelOut)
print(evals)

# ---- 7. Importancia de variables ----

var_import <- biomod2::get_variables_importance(myBiomodModelOut)

# Calcular importancia promedio por variable
library(dplyr)
var_summary <- var_import %>%
  group_by(expl.var) %>%
  summarise(mean_importance = mean(var.imp, na.rm = TRUE)) %>%
  arrange(desc(mean_importance))

print(var_summary)

# Graficar
barplot(var_summary$mean_importance, 
        names.arg = var_summary$expl.var,
        las = 2, 
        main = "Importancia promedio de variables", 
        col = "steelblue",
        ylab = "Importancia relativa",
        cex.names = 0.9)

# O por algoritmo
var_by_algo <- var_import %>%
  group_by(expl.var, algo) %>%
  summarise(mean_importance = mean(var.imp, na.rm = TRUE), .groups = "drop")

# Graficar con ggplot2
library(ggplot2)
ggplot(var_by_algo, aes(x = reorder(expl.var, mean_importance), 
                        y = mean_importance, 
                        fill = algo)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Importancia de variables por algoritmo",
       x = "Variable", y = "Importancia", fill = "Modelo") +
  theme_minimal()


# ---- 8. Proyección de modelos ----
myBiomodProjection <- biomod2::BIOMOD_Projection(
  bm.mod = myBiomodModelOut,           # Cambio: modeling.output -> bm.mod
  proj.name = "current",
  new.env = env,
  models.chosen = "all",               # Cambio: selected.models -> models.chosen
  metric.binary = "TSS",               # Cambio: binary.meth -> metric.binary
  build.clamping.mask = TRUE
)


# ---- 9. Visualización de predicciones individuales ----

# Obtener predicciones
predicciones <- get_predictions(myBiomodProjection)

# Ver qué modelos tienes
names(predicciones)

# Graficar todas las predicciones continuas
plot(predicciones, main = names(preddicciones))

# Graficar modelos individuales con presencias
par(mfrow = c(2, 2))  # 2x2 grid
for(i in 1:nlyr(predicciones)) {
  plot(pred_continuas[[i]], 
       main = paste("Continuo:", names(predicciones)[i]),
       col = terrain.colors(100))
  points(pres$lon, pres$lat, pch = 20, col = "red", cex = 0.5)
}
par(mfrow = c(1, 1))  # Resetear

# Obtener predicciones binarias (usando umbral TSS)
# Cargar proyecciones binarias
proj_dir <- paste0(getwd(), "/Calomys_venustus/proj_current")
bin_files <- list.files(proj_dir, pattern = ".*_TSSbin.tif$", full.names = TRUE)

if(length(bin_files) > 0) {
  pred_binarias <- rast(bin_files)
  
  # Graficar binarias
  par(mfrow = c(2, 2))
  for(i in 1:nlyr(pred_binarias)) {
    plot(pred_binarias[[i]], 
         main = paste("Binario:", names(pred_binarias)[i]),
         col = c("gray90", "darkgreen"),
         legend = FALSE)
    points(pres$lon, pres$lat, pch = 20, col = "red", cex = 0.5)
  }
  par(mfrow = c(1, 1))
}

# ---- 10. Ensamble  ----

# Identificar buenos algoritmos
buenos_algos <- evals %>%
  filter(metric.eval == "TSS", "cutoff" > 700) %>%
  pull(algo) %>%
  unique()
head(evals)
# Usar esos algoritmos en versión allData
models_finales <- paste0("Calomys.venustus_PA1_RUN1_", buenos_algos)

myBiomodEnsemble <- biomod2::BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut,
  models.chosen = models_finales,
  em.by = "all",
  em.algo = c("EMmean", "EMwmean"),
  metric.eval = c("TSS", "ROC"),
  var.import = 3
)

# Ver qué modelos se incluyeron en el ensemble
myBiomodEnsemble

# ---- 11. Proyección del ensemble ----
myBiomodEnsembleProj <- biomod2::BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEnsemble,
  proj.name = "current",
  new.env = env,
  models.chosen = "all",
  metric.binary = c("TSS", "ROC")
)

# ---- 12. Visualización del ensemble ----

# Obtener predicciones del ensemble
ensemble_pred <- get_predictions(myBiomodEnsembleProj)
names(ensemble_pred)

# Graficar ensemble continuo
plot(ensemble_pred, 
     main = names(ensemble_pred),
     col = terrain.colors(100))


# Comparar media simple vs media ponderada
par(mfrow = c(1, 2))
# EMmean (media simple)
plot(ensemble_pred[[grep("EMmean", names(ensemble_pred))[1]]], 
     main = "Ensemble: Media simple",
     col = terrain.colors(100))

# EMwmean (media ponderada por TSS)
plot(ensemble_pred[[grep("EMwmean", names(ensemble_pred))[1]]], 
     main = "Ensemble: Media ponderada",
     col = terrain.colors(100))

par(mfrow = c(1, 1))

# ---- 13. Exportar mapas seleccionados ----

# Crear carpeta para exportar
dir.create("mapas_exportados", showWarnings = FALSE)


# Exportar ensemble (media ponderada - RECOMENDADO)
ensemble_wmean <- ensemble_pred[[grep("EMwmean", names(ensemble_pred))[1]]]
writeRaster(ensemble_wmean, 
            "mapas_exportados/ENSEMBLE_media_ponderada.tif", 
            overwrite = TRUE)

# Exportar ensemble (media simple)
ensemble_mean <- ensemble_pred[[grep("EMmean", names(ensemble_pred))[1]]]
writeRaster(ensemble_mean, 
            "mapas_exportados/ENSEMBLE_media_simple.tif", 
            overwrite = TRUE)

# ---- 15. Resumen de evaluación del ensemble ----
# Ver evaluaciones del ensemble
ensemble_eval <- get_evaluations(myBiomodEnsemble)
print(ensemble_eval)

# Comparar con modelos individuales
cat("\n=== RESUMEN DE MODELOS ===\n")
cat("\nModelos individuales (TSS):\n")
print(evals %>% 
        filter(metric.eval == "TSS", run == "RUN1") %>%
        select(algo, calibration, validation))

cat("\nEnsemble (TSS):\n")
print(ensemble_eval %>% 
        filter(metric.eval == "TSS") %>%
        select(full.name, calibration))

# ---- 16. Guardar workspace completo ----
save(myBiomodModelOut, 
     myBiomodProjection, 
     myBiomodEnsemble, 
     myBiomodEnsembleProj,
     file = "biomod2_workspace_completo.RData")

# ---- 17. Guardar objetos individuales en RDS ----

# Crear carpeta para RDS
dir.create("modelos_rds", showWarnings = FALSE)

# Guardar modelos individuales
saveRDS(myBiomodModelOut, "modelos_rds/biomod_modelos.rds")
saveRDS(myBiomodProjection, "modelos_rds/biomod_proyeccion_current.rds")
saveRDS(myBiomodEnsemble, "modelos_rds/biomod_ensemble.rds")
saveRDS(myBiomodEnsembleProj, "modelos_rds/biomod_ensemble_proyeccion.rds")

# Guardar evaluaciones y resultados
saveRDS(evals, "modelos_rds/evaluaciones_modelos.rds")
saveRDS(var_import, "modelos_rds/importancia_variables.rds")

# Guardar predicciones como rasters en RDS
saveRDS(predicciones, "modelos_rds/predicciones_continuas.rds")
saveRDS(ensemble_pred, "modelos_rds/ensemble_predicciones.rds")

cat("\n✓ Modelos guardados en RDS en carpeta 'modelos_rds/'\n")

# ---- Para cargar después ----
# myBiomodModelOut <- readRDS("modelos_rds/biomod_modelos.rds")
# myBiomodEnsemble <- readRDS("modelos_rds/biomod_ensemble.rds")