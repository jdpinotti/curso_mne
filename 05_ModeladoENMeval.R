##############################################
# Modelado con ENMeval 2.0
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################

# Cargar las librerías necesarias
library(terra)
library(dplyr)
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")

# Definir la ruta de los archivos
csv_path <- "cvenusus_clean_thin1.csv"

datos <- read.csv(csv_path)

tif_path <- list.files("/corte_m", pattern = "*.tif", full.names = TRUE)


#################Generacion bg ##############################
library(sp)
library(dismo)
library(sf)
mask <- st_read("mvenustus.shp")

mask_sp <- as_Spatial(mask)

# Ruta al directorio donde están los rasters

ruta <- "corte_m"

# Cargar todos los .tif del directorio
envs <- rast(list.files(ruta, pattern = "\\.tif$", full.names = TRUE))
envs
names(envs)

# Función para simplificar nombres
renombrar_bio <- function(nombres_originales) {
  # Extraer el número bio usando expresiones regulares
  numeros_bio <- gsub(".*bio(\\d+).*", "\\1", nombres_originales)
  # Crear nuevos nombres
  paste0("bio", numeros_bio)
}

# Aplicar la función
names(envs) <- renombrar_bio(names(envs))
print(names(envs))

# Generar 10000 puntos aleatorios dentro de la máscara
set.seed(42)  # Fijar la semilla para reproducibilidad
bg_points <- st_sample(mask, size = 10000, type = "random")

# Extraer valores ambientales en los puntos de background
bg_values <- terra::extract(envs, vect(bg_points))

head(bg_values)
dim(bg_values)

# Crear un DataFrame para el background y guardar en un CSV
bg_full <- data.frame(st_coordinates(bg_points), bg_values) %>% 
  dplyr::select(!3)

write.csv(bg_full, "background_cvenustus.csv", row.names = FALSE)


###### Modelo con las variables selecionadas por PCA ######

occs <- read.csv("var_amb_presencias.csv")
names(occs)
# Simplificar nombres
nombres_finales <- c("sp", "long", "lat", 
                     "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9",
                     "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", 
                     "bio18", "bio19")

names(occs) <- nombres_finales
print(names(occs))

#Seleccionar las variables de acuerdo al PCA

occs <- dplyr::select(occs, 
                      long, 
                      lat,
                      bio03 = bio3,
                      bio04 = bio4,
                      bio06 = bio6,
                      bio07 = bio7,
                      bio17 = bio17)

names(occs)
names(bg_full)
bg <- dplyr::select(bg_full,
                    long = X,
                    lat = Y,
                    bio03 = bio3,
                    bio04 = bio4,
                    bio06 = bio6,
                    bio07 = bio7,
                    bio17 = bio17)

names(bg)
names(occs)
###################### MODELADO #########################

library(ENMeval)


cvenustus_model <- ENMevaluate(
  occs = occs,
  envs = NULL,            # No se usan capas en RasterStack para el enfoque SWD
  bg = bg,                # Background con los mismos valores ambientales
  algorithm = "maxnet", 
  partitions = 'block',
  user.eval = NULL,
  doClamp = TRUE,
  tune.args = list(fc = c("L", "Q", "LQ", "LQH", "LQHT"), rm = c(0.5, 1, 1.5, 2, 5)),
  taxon.name = "Calomys venustus"
)

#c("L", "Q","LQHPT", "LQ", "LQH", "LQHT", "LQHP", "LQP")

#para ir viendo resultados
show(cvenustus_model)
eval.results(cvenustus_model)
write.csv(cvenustus_model@results, "modelos_cvenustus_enmeval.csv")


####Elegir el mejor modelo ####

# Obtener el modelo correspondiente a la mejor fila según AICc
best_row <- which.min(cvenustus_model@results$delta.AICc)
best_model <- eval.models(cvenustus_model)[[best_row]]

# Ver parámetros del mejor modelo
best_params <- cvenustus_model@results[best_row, c("fc", "rm", "tune.args")]
print(best_params)

# Opción 2: Acceder directamente al modelo usando los parámetros óptimos
# En este caso: fc = LQH, rm = 0.5 (fila 4)
#best_model <- eval.models(cvenustus_model)[[4]]

# Resumir resultados del mejor modelo
cat("Mejor modelo según AICc:\n")
cat("Feature Classes:", best_params$fc, "\n")
cat("Regularization Multiplier:", best_params$rm, "\n")
cat("AUC Validación (media):", cvenustus_model@results[best_row, "auc.val.avg"], "\n")
cat("CBI Validación (media):", cvenustus_model@results[best_row, "cbi.val.avg"], "\n")
cat("Delta AICc:", cvenustus_model@results[best_row, "delta.AICc"], "\n")
cat("Número de coeficientes:", cvenustus_model@results[best_row, "ncoef"], "\n")

# Guardar el mejor modelo
saveRDS(best_model, "cvenustus_enmeval.rds")

#### EXPLORAR TODOS LOS RESULTADOS ####  ver bien desde aca. falta importancia y aporte de las variables#######

# 1. VER ESTRUCTURA DEL OBJETO ENMevaluation
slotNames(cvenustus_model)
cvenustus_model@algorithm
cvenustus_model@variable.importance


# 4. IMPORTANCIA DE VARIABLES (Permutation Importance)
# ENMeval calcula esto automáticamente
var.importance <- cvenustus_model@variable.importance
var.importance
write.csv(var.importance, "variable_importance_cvenustus.csv")

# 6. UMBRALES ÓPTIMOS
# Threshold 10% presencia de entrenamiento (10p)
thresh_10p <- quantile(pred_occs, 0.1)
cat("Umbral 10% presencias:", thresh_10p, "\n")

# Maximum Training Presence (MTP)
thresh_mtp <- min(pred_occs)
cat("Umbral MTP:", thresh_mtp, "\n")

# Equal sensitivity-specificity (ESS)
coords_roc <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
thresh_ess <- coords_roc$threshold
cat("Umbral ESS (mejor sensibilidad-especificidad):", thresh_ess, "\n")

# Sensibilidad y especificidad en umbrales
sens_spec <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity"))
head(sens_spec, 10)

# 7. OBTENER COEFICIENTES / PESOS DE LAS VARIABLES
cat("\n=== COEFICIENTES DEL MODELO ===\n")
betas <- best_model$betas
print(betas)

# Crear tabla interpretable
coef_table <- data.frame(
  Variable = names(betas[-1]),  # Excluir el intercepto
  Coeficiente = betas[-1]
)
coef_table <- coef_table[order(abs(coef_table$Coeficiente), decreasing = TRUE), ]
print(coef_table)
write.csv(coef_table, "coeficientes_modelo.csv", row.names = FALSE)

# 8. RESUMEN COMPLETO DEL MEJOR MODELO
cat("\n=== RESUMEN DEL MEJOR MODELO ===\n")
print(summary(best_model))

# 9. MÉTRICAS DE EVALUACIÓN DEL MEJOR MODELO
best_metrics <- cvenustus_model@results[best_idx, ]
cat("\n=== MÉTRICAS DE EVALUACIÓN ===\n")
print(best_metrics)

# 10. DESEMPEÑO POR FOLD (cross-validation)
# Extraer predicciones de validación
pred_table <- cvenustus_model@predictions[[best_idx]]
head(pred_table)

# Calcular AUC por cada fold
folds <- unique(pred_table$fold)
for(fold in folds) {
  fold_data <- pred_table[pred_table$fold == fold, ]
  fold_response <- c(rep(1, sum(fold_data$data_type == "presence")), 
                     rep(0, sum(fold_data$data_type == "background")))
  fold_pred <- fold_data$prediction
  fold_auc <- auc(fold_response, fold_pred)
  cat("AUC Fold", fold, ":", fold_auc, "\n")
}

# 11. TABLA COMPLETA CON TODOS LOS UMBRALES
umbrales <- data.frame(
  Umbral = c("10% Presencia", "MTP", "ESS"),
  Valor = c(thresh_10p, thresh_mtp, thresh_ess),
  Descripcion = c(
    "10% de presencias de entrenamiento más baja",
    "Presencia mínima de entrenamiento",
    "Mejor balance sensibilidad-especificidad"
  )
)
print(umbrales)
write.csv(umbrales, "umbrales_optimos.csv", row.names = FALSE)

# 12. EXPORTAR TODO A UN ARCHIVO DE TEXTO
sink("reporte_completo_cvenustus.txt")
cat("=== REPORTE COMPLETO DEL MODELO - Calomys venustus ===\n\n")
cat("Mejores parámetros:\n")
print(cvenustus_model@results[best_idx, c("fc", "rm", "tune.args")])
cat("\n\nMétricas de evaluación:\n")
print(best_metrics)
cat("\n\nCoeficientes del modelo:\n")
print(coef_table)
cat("\n\nUmbrales óptimos:\n")
print(umbrales)
cat("\n\nAUC de validación:", best_metrics$auc.val.avg, "\n")
cat("CBI de validación:", best_metrics$cbi.val.avg, "\n")
sink()

