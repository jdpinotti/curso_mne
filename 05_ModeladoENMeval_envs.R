##############################################
# Modelado con ENMeval 2.0
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################

# Cargar las librerías necesarias
library(terra)
library(dplyr)
library(ENMeval)
library(pROC)
library(ggplot2)
library(maxnet)
setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")

# Definir la ruta de los archivos
csv_path <- "cvenustus_clean_thin1.csv"
datos <- read.csv(csv_path)

library(sp)
library(dismo)
library(sf)
mask <- st_read("mvenustus.shp")
mask_sp <- as_Spatial(mask)

# Ruta al directorio donde están los rasters
ruta <- "corte_m"

# Cargar todos los .tif del directorio
envs <- rast(list.files(ruta, pattern = "\\.tif$", full.names = TRUE))
names(envs)
# Función para simplificar nombres
renombrar_bio <- function(nombres_originales) {
  numeros_bio <- gsub(".*bio(\\d+).*", "\\1", nombres_originales)
  paste0("bio", numeros_bio)
}

# Aplicar la función
names(envs) <- renombrar_bio(names(envs))
print(names(envs))

###### Modelo con las variables seleccionadas por PCA ######

occs <- read.csv("outputs/valores_presencias.csv")
names(occs)
# Simplificar nombres
nombres_finales <- c("sp", "long", "lat", 
                     "bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17",
                     "bio18", "bio19", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                     "bio8", "bio9")
names(occs) <- nombres_finales

#Seleccionar las variables de acuerdo al PCA
occs_coords <- occs[, c("long", "lat")]
names(occs_coords) <- c("x", "y")

# Seleccionamos las mismas variables en el raster
envs_sub <- envs[[c("bio3", "bio4", "bio6", "bio7", "bio17")]]


###################### MODELADO #########################

cvenustus_model <- ENMevaluate(
  occs = occs_coords,
  envs = envs_sub,
  algorithm = "maxnet", 
  partitions = 'block',
  doClamp = TRUE,
  raster.preds = FALSE,
  tune.args = list(fc = c("L", "Q", "LQ", "LQH", "LQHT"),
                   rm = c(0.5, 1, 1.5, 2, 5)),
  taxon.name = "Calomys venustus"
)


#Para ir viendo resultados
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

# Resumir resultados del mejor modelo
cat("Mejor modelo según AICc:", best_row, "\n")
cat("Feature Classes:", best_params$fc, "\n")
cat("Regularization Multiplier:", best_params$rm, "\n")
cat("AUC Validación (media):", cvenustus_model@results[best_row, "auc.val.avg"], "\n")
cat("CBI Validación (media):", cvenustus_model@results[best_row, "cbi.val.avg"], "\n")
cat("Delta AICc:", cvenustus_model@results[best_row, "delta.AICc"], "\n")
cat("Número de coeficientes:", cvenustus_model@results[best_row, "ncoef"], "\n")

# Guardar el mejor modelo
saveRDS(best_model, "cvenustus_enmeval.rds")


###################### Curvas de Respuesta del Modelo seleccionado #########################

mod.seq <- eval.models(cvenustus_model)[[4]]  #colocar el modelo a evaluar aca

plot(mod.seq, type = "cloglog")


###################### IMPORTANCIA DE VARIABLES #########################

# Extraer las variables base
extraer_var_base <- function(nombres) {
  var_base <- nombres
  
  # Para funciones hinge: extraer solo el nombre de la variable
  var_base <- gsub("hinge\\(([^)]+)\\):.*", "\\1", var_base)
  
  # Para transformaciones cuadráticas I(var^2)
  var_base <- gsub("I\\(([^)^]+).*", "\\1", var_base)
  
  # Limpiar cualquier residuo
  var_base <- trimws(var_base)
  
  return(var_base)
}

# Aplicar la función
contrib_df <- data.frame(
  termino = names(var_importance),
  coeficiente = var_importance,
  variable_base = extraer_var_base(names(var_importance))
)

# Ver cómo quedó la extracción
print(contrib_df)

# Calcular contribución por variable
contrib_por_variable <- contrib_df %>%
  group_by(variable_base) %>%
  summarise(
    contribucion_total = sum(abs(coeficiente)),
    n_terminos = n()
  ) %>%
  mutate(
    porcentaje = (contribucion_total / sum(contribucion_total)) * 100
  ) %>%
  arrange(desc(porcentaje))

print(contrib_por_variable)

# Gráfico
ggplot(contrib_por_variable, aes(x = reorder(variable_base, porcentaje), 
                                 y = porcentaje)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = paste0(round(porcentaje, 1), "%")), 
            hjust = -0.2, size = 4.5, fontface = "bold") +
  coord_flip() +
  labs(title = "Contribución de Variables al Modelo",
       subtitle = "Calomys venustus - Modelo Maxnet",
       x = "Variable Bioclimática",
       y = "Contribución Relativa (%)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_blank()
  ) +
  ylim(0, max(contrib_por_variable$porcentaje) * 1.2)

ggsave("contribucion_variables.png", width = 10, height = 6, dpi = 300)

###################### IMPORTANCIA DE PERMUTACION #########################

# Función para Permutación
calcular_importancia_permutacion <- function(modelo, envs, occs_coords, n_permutaciones = 10) {
  
  cat("Calculando importancia mediante Permutación...\n")
  cat("Número de permutaciones por variable:", n_permutaciones, "\n\n")
  
  # Extraer valores ambientales
  occs_values <- terra::extract(envs, occs_coords, ID = FALSE)
  
  # Predicción baseline
  pred_baseline <- predict(modelo, occs_values, type = "cloglog")
  auc_baseline <- mean(pred_baseline, na.rm = TRUE)
  
  variables <- names(envs)
  resultados <- data.frame(
    variable = character(),
    auc_promedio_permutado = numeric(),
    perdida_auc = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Evaluar cada variable
  for (var in variables) {
    cat("  Permutando", var, "...\n")
    
    auc_permutados <- numeric(n_permutaciones)
    
    # Realizar múltiples permutaciones
    for (i in 1:n_permutaciones) {
      # Copiar datos y permutar la variable
      datos_perm <- occs_values
      datos_perm[[var]] <- sample(datos_perm[[var]])
      
      # Predecir con variable permutada
      pred_perm <- predict(modelo, datos_perm, type = "cloglog")
      auc_permutados[i] <- mean(pred_perm, na.rm = TRUE)
    }
    
    # Calcular promedio de AUC permutado
    auc_prom <- mean(auc_permutados)
    perdida <- auc_baseline - auc_prom
    
    resultados <- rbind(resultados, data.frame(
      variable = var,
      auc_promedio_permutado = auc_prom,
      perdida_auc = perdida
    ))
  }
  
  # Calcular porcentaje
  resultados <- resultados %>%
    mutate(
      perdida_positiva = pmax(perdida_auc, 0),
      importancia_porcentaje = (perdida_positiva / sum(perdida_positiva)) * 100
    ) %>%
    arrange(desc(importancia_porcentaje))
  
  return(resultados)
}

# EJECUTAR PERMUTACIÓN
importancia_perm <- calcular_importancia_permutacion(
  best_model, 
  envs_sub, 
  occs_coords, 
  n_permutaciones = 100  # Puedes ajustar este número
)

print(importancia_perm)

# Guardar
write.csv(importancia_perm, "importancia_permutacion.csv", row.names = FALSE)

# VISUALIZAR
ggplot(importancia_perm, aes(x = reorder(variable, importancia_porcentaje), 
                             y = importancia_porcentaje)) +
  geom_col(fill = "#D32F2F", width = 0.7) +
  geom_text(aes(label = paste0(round(importancia_porcentaje, 1), "%")), 
            hjust = -0.2, size = 5, fontface = "bold") +
  coord_flip() +
  labs(title = "Importancia de Variables - Método Permutación",
       subtitle = "Pérdida de precisión al permutar cada variable (n=20)",
       x = NULL,
       y = "Importancia (%)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("importancia_permutacion.png", width = 10, height = 6, dpi = 300)



# ============================================
# CÁLCULO DE UMBRALES
# ============================================

# ---- 1. GENERAR PUNTOS DE BACKGROUND ----
cat("Generando puntos de background...\n")
set.seed(123)  # Para reproducibilidad

# Generar puntos aleatorios en el área de estudio
bg_points <- spatSample(envs_sub, size = 10000, method = "random", 
                        na.rm = TRUE, xy = TRUE, values = TRUE)

# Convertir a dataframe
bg <- as.data.frame(bg_points)

# Verificar
cat("Puntos de background generados:", nrow(bg), "\n")
cat("Variables en background:", paste(names(bg), collapse = ", "), "\n\n")

# ---- 2. PREDICCIONES ----
cat("Calculando predicciones...\n")

# Definir las variables del modelo
vars <- c("bio3", "bio4", "bio6", "bio7", "bio17")

# Verificar que las variables existan
if(!all(vars %in% names(occs))) {
  stop("Error: No todas las variables están en 'occs'")
}
if(!all(vars %in% names(bg))) {
  stop("Error: No todas las variables están en 'bg'")
}

# Predicciones en presencias
pred_occs <- predict(best_model, occs[, vars], type = "cloglog")
cat("Predicciones en presencias:", length(pred_occs), "\n")

# Predicciones en background
pred_bg <- predict(best_model, bg[, vars], type = "cloglog")
cat("Predicciones en background:", length(pred_bg), "\n\n")

# ---- 3. CALCULAR UMBRALES ----
cat("Calculando umbrales...\n")

# ROC para calcular umbrales
library(pROC)
all_preds <- c(pred_occs, pred_bg)
all_obs <- c(rep(1, length(pred_occs)), rep(0, length(pred_bg)))
roc_obj <- roc(all_obs, all_preds, quiet = TRUE)

cat("AUC:", round(auc(roc_obj), 4), "\n\n")

# 1. Threshold 10% presencia de entrenamiento (10p)
thresh_10p <- quantile(pred_occs, 0.1, na.rm = TRUE)
cat("Umbral 10% presencias:", round(thresh_10p, 4), "\n")

# 2. Maximum Training Presence (MTP)
thresh_mtp <- min(pred_occs, na.rm = TRUE)
cat("Umbral MTP:", round(thresh_mtp, 4), "\n")

# 3. Equal sensitivity-specificity (ESS)
coords_roc <- coords(roc_obj, "best", best.method = "closest.topleft", 
                     ret = c("threshold", "sensitivity", "specificity"))
thresh_ess <- coords_roc$threshold
cat("Umbral ESS (mejor sensibilidad-especificidad):", round(thresh_ess, 4), "\n")
cat("  Sensibilidad:", round(coords_roc$sensitivity, 4), "\n")
cat("  Especificidad:", round(coords_roc$specificity, 4), "\n\n")

# 4. Maximum Training Sensitivity plus Specificity (MaxSSS)
coords_maxsss <- coords(roc_obj, "best", best.method = "youden", 
                        ret = c("threshold", "sensitivity", "specificity"))
thresh_maxsss <- coords_maxsss$threshold
cat("Umbral MaxSSS (Youden):", round(thresh_maxsss, 4), "\n")
cat("  Sensibilidad:", round(coords_maxsss$sensitivity, 4), "\n")
cat("  Especificidad:", round(coords_maxsss$specificity, 4), "\n\n")

# ---- 4. GUARDAR UMBRALES ----
umbrales <- data.frame(
  Metodo = c("10% Training Presence", 
             "Minimum Training Presence (MTP)", 
             "Equal Sensitivity-Specificity",
             "Maximum SSS (Youden)"),
  Umbral = c(thresh_10p, thresh_mtp, thresh_ess, thresh_maxsss),
  Sensibilidad = c(NA, NA, coords_roc$sensitivity, coords_maxsss$sensitivity),
  Especificidad = c(NA, NA, coords_roc$specificity, coords_maxsss$specificity)
)

print(umbrales)
write.csv(umbrales, "umbrales_cvenustus.csv", row.names = FALSE)

# ---- 5. VISUALIZAR DISTRIBUCIÓN DE PREDICCIONES ----
# Crear dataframe para ggplot
pred_df <- data.frame(
  valor = c(pred_occs, pred_bg),
  tipo = c(rep("Presencia", length(pred_occs)), 
           rep("Background", length(pred_bg)))
)


###################   Proyeccion #########

# ---- 7. APLICAR UMBRALES AL RASTER ----
cat("\nGenerando mapas binarios con umbrales...\n")

# Proyectar modelo
pred_raster <- predict(envs_sub, best_model, na.rm=TRUE, type = "cloglog")

writeRaster(pred_raster, "modelo_continuo.tif", overwrite = TRUE)

# Crear mapas binarios
binary_10p <- pred_raster >= thresh_10p
binary_mtp <- pred_raster >= thresh_mtp
binary_ess <- pred_raster >= thresh_ess
binary_maxsss <- pred_raster >= thresh_maxsss

# Guardar mapas binarios
writeRaster(binary_10p, "binario_10p.tif", overwrite = TRUE)
writeRaster(binary_mtp, "binario_mtp.tif", overwrite = TRUE)
writeRaster(binary_ess, "binario_ess.tif", overwrite = TRUE)
writeRaster(binary_maxsss, "binario_maxsss.tif", overwrite = TRUE)

# Visualizar mapas
par(mfrow = c(2, 2))

plot(binary_10p, main = "Umbral 10% Presencias",
     col = c("gray90", "darkgreen"), legend = FALSE)
points(occs_coords$x, occs_coords$y, pch = 20, cex = 0.5)

plot(binary_mtp, main = "Umbral MTP",
     col = c("gray90", "darkgreen"), legend = FALSE)
points(occs_coords$x, occs_coords$y, pch = 20, cex = 0.5)

plot(binary_ess, main = "Umbral ESS",
     col = c("gray90", "darkgreen"), legend = FALSE)
points(occs_coords$x, occs_coords$y, pch = 20, cex = 0.5)

plot(binary_maxsss, main = "Umbral MaxSSS",
     col = c("gray90", "darkgreen"), legend = FALSE)
points(occs_coords$x, occs_coords$y, pch = 20, cex = 0.5)

par(mfrow = c(1, 1))

# Guardar figura
png("mapas_binarios_umbrales.png", width = 3000, height = 3000, res = 300)

