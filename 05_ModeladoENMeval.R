##############################################
# Modelado con ENMeval 2.0 + Importancia de Variables
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################

# Cargar las librerías necesarias
library(terra)
library(dplyr)
library(ENMeval)
library(pROC)
library(ggplot2)

setwd("/home/dargwind/Modelado de Nicho/curso_modelado/wd")

# Definir la ruta de los archivos
csv_path <- "cvenustus_clean_thin1.csv"
datos <- read.csv(csv_path)

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

# Función para simplificar nombres
renombrar_bio <- function(nombres_originales) {
  numeros_bio <- gsub(".*bio(\\d+).*", "\\1", nombres_originales)
  paste0("bio", numeros_bio)
}

# Aplicar la función
names(envs) <- renombrar_bio(names(envs))
print(names(envs))

# Generar 10000 puntos aleatorios dentro de la máscara
set.seed(42)
bg_points <- st_sample(mask, size = 10000, type = "random")

# Extraer valores ambientales en los puntos de background
bg_values <- terra::extract(envs, vect(bg_points))

# Crear un DataFrame para el background y guardar en un CSV
bg_full <- data.frame(st_coordinates(bg_points), bg_values) %>% 
  dplyr::select(!3)

write.csv(bg_full, "background_cvenustus.csv", row.names = FALSE)


###### Modelo con las variables seleccionadas por PCA ######

occs <- read.csv("var_amb_presencias.csv")
# Simplificar nombres
nombres_finales <- c("sp", "long", "lat", 
                     "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9",
                     "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", 
                     "bio18", "bio19")
names(occs) <- nombres_finales

#Seleccionar las variables de acuerdo al PCA
occs <- dplyr::select(occs, 
                      long, 
                      lat,
                      bio03 = bio3,
                      bio04 = bio4,
                      bio06 = bio6,
                      bio07 = bio7,
                      bio17 = bio17)

bg <- read.csv("background_cvenustus.csv")
bg <- dplyr::select(bg,
                    long = X,
                    lat = Y,
                    bio03 = bio3,
                    bio04 = bio4,
                    bio06 = bio6,
                    bio07 = bio7,
                    bio17 = bio17)

###################### MODELADO #########################

cvenustus_model <- ENMevaluate(
  occs = occs,
  envs = NULL,
  bg = bg,
  algorithm = "maxnet", 
  partitions = 'block',
  user.eval = NULL,
  doClamp = TRUE,
  tune.args = list(fc = c("L", "Q", "LQ", "LQH", "LQHT"), rm = c(0.5, 1, 1.5, 2, 5)),
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


###################### IMPORTANCIA DE VARIABLES #########################

# Variables ambientales (sin lon/lat)
vars <- c("bio03", "bio04", "bio06", "bio07", "bio17")

# ---- 2. IMPORTANCIA POR CONTRIBUCIÓN ----
# Basada en los coeficientes del modelo maxnet
calc_contribution <- function(model, vars) {
  # Extraer coeficientes del modelo maxnet
  # En maxnet, los coeficientes están en model$betas
  coefs <- model$betas
  
  # Eliminar el intercepto si existe
  if ("(Intercept)" %in% names(coefs)) {
    coefs <- coefs[names(coefs) != "(Intercept)"]
  }
  
  # Agrupar por variable base (sin transformaciones)
  var_contrib <- sapply(vars, function(v) {
    # Buscar todas las features que contengan esta variable
    idx <- grep(v, names(coefs), fixed = TRUE)
    if (length(idx) > 0) {
      sum(abs(coefs[idx]))
    } else {
      0
    }
  })
  
  # Normalizar a porcentaje
  var_contrib_pct <- (var_contrib / sum(var_contrib)) * 100
  
  return(sort(var_contrib_pct, decreasing = TRUE))
}

contribution <- calc_contribution(best_model, vars)

cat("\n=== IMPORTANCIA POR CONTRIBUCIÓN ===\n")
print(round(contribution, 2))

# ---- 3. IMPORTANCIA POR PERMUTACIÓN ----
# Evalúa cómo cambia el AUC al permutar cada variable

calc_permutation <- function(model, occs, bg, vars, n_permutations = 10) {
  
  # Combinar datos (presencias = 1, background = 0)
  occs_data <- occs[, vars]
  bg_data <- bg[, vars]
  
  all_data <- rbind(occs_data, bg_data)
  y <- c(rep(1, nrow(occs_data)), rep(0, nrow(bg_data)))
  
  # Predicción original
  pred_original <- predict(model, all_data, type = "cloglog")
  
  # AUC original
  library(pROC)
  auc_original <- auc(roc(y, pred_original, quiet = TRUE))
  
  # Permutar cada variable
  var_importance <- numeric(length(vars))
  names(var_importance) <- vars
  
  for (i in seq_along(vars)) {
    var <- vars[i]
    
    auc_permuted <- numeric(n_permutations)
    
    for (j in 1:n_permutations) {
      # Copiar datos
      data_perm <- all_data
      
      # Permutar la variable
      data_perm[, var] <- sample(data_perm[, var])
      
      # Predecir con variable permutada
      pred_perm <- predict(model, data_perm, type = "cloglog")
      
      # Calcular AUC
      auc_permuted[j] <- auc(roc(y, pred_perm, quiet = TRUE))
    }
    
    # Importancia = reducción promedio en AUC
    var_importance[i] <- auc_original - mean(auc_permuted)
  }
  
  # Normalizar a porcentaje
  var_importance_pct <- (var_importance / sum(var_importance)) * 100
  
  return(sort(var_importance_pct, decreasing = TRUE))
}

cat("\nCalculando importancia por permutación (puede tardar)...\n")
permutation <- calc_permutation(best_model, occs, bg, vars, n_permutations = 10)

cat("\n=== IMPORTANCIA POR PERMUTACIÓN ===\n")
print(round(permutation, 2))

# ---- 4. RESUMEN COMBINADO ----
importance_df <- data.frame(
  Variable = names(contribution),
  Contribution = as.numeric(contribution),
  Permutation = permutation[names(contribution)]
)

# Ordenar por contribución
importance_df <- importance_df[order(-importance_df$Contribution), ]

cat("\n=== TABLA RESUMEN ===\n")
print(importance_df)

# Guardar resultados
write.csv(importance_df, "variable_importance_cvenustus.csv", row.names = FALSE)

# ---- 5. VISUALIZACIÓN ----

# Preparar datos para gráfico
library(tidyr)
importance_long <- pivot_longer(
  importance_df,
  cols = c(Contribution, Permutation),
  names_to = "Metric",
  values_to = "Importance"
)

# Nombres más descriptivos
var_labels <- c(
  "bio03" = "Isotermalidad (BIO3)",
  "bio04" = "Estacionalidad T° (BIO4)",
  "bio06" = "T° Mín Mes Frío (BIO6)",
  "bio07" = "Rango Anual T° (BIO7)",
  "bio17" = "Precip Trimestre Seco (BIO17)"
)

importance_long$Variable <- factor(
  importance_long$Variable,
  levels = importance_df$Variable
)

# Gráfico de barras agrupadas
p1 <- ggplot(importance_long, aes(x = Variable, y = Importance, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c("Contribution" = "#2E86AB", "Permutation" = "#A23B72"),
    labels = c("Contribución", "Permutación")
  ) +
  scale_x_discrete(labels = var_labels) +
  labs(
    title = "Importancia de Variables Ambientales",
    subtitle = "Calomys venustus - Modelo Maxnet",
    x = "",
    y = "Importancia (%)",
    fill = "Método"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

print(p1)
ggsave("importancia_variables_barplot.png", p1, width = 8, height = 6, dpi = 300)


# ---- 6. CURVAS DE RESPUESTA ----
# Mostrar cómo responde el modelo a cada variable

# Definir nombres descriptivos
var_labels <- c(
  "bio03" = "Isotermalidad (BIO3)",
  "bio04" = "Estacionalidad T° (BIO4)",
  "bio06" = "T° Mín Mes Frío (BIO6)",
  "bio07" = "Rango Anual T° (BIO7)",
  "bio17" = "Precip Trimestre Seco (BIO17)"
)

plot_response_curves <- function(model, data, vars, var_labels) {
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  
  for (var in vars) {
    # Crear secuencia de valores para la variable
    var_range <- seq(min(data[, var], na.rm = TRUE), 
                     max(data[, var], na.rm = TRUE), 
                     length.out = 100)
    
    # Mantener las otras variables en su mediana
    new_data <- data.frame(
      matrix(
        rep(apply(data[, vars], 2, median, na.rm = TRUE), each = 100),
        nrow = 100
      )
    )
    names(new_data) <- vars
    new_data[, var] <- var_range
    
    # Predecir
    pred <- predict(model, new_data, type = "cloglog")
    
    # Plotear
    plot(var_range, pred, type = "l", lwd = 2, col = "#2E86AB",
         xlab = var_labels[var], ylab = "Idoneidad",
         main = paste("Respuesta a", var_labels[var]))
    grid()
  }
  
  par(mfrow = c(1, 1))
}

# Datos combinados para las curvas
all_env_data <- rbind(occs[, vars], bg[, vars])

# Verificar que hay datos
cat("\nDimensiones de datos para curvas de respuesta:", dim(all_env_data), "\n")

png("response_curves_cvenustus.png", width = 1200, height = 800, res = 120)
plot_response_curves(best_model, all_env_data, vars, var_labels)
dev.off()

cat("Curvas de respuesta guardadas en: response_curves_cvenustus.png\n")

# ---- 6. PREDICCIÓN Y UMBRALES ----

# Predicciones en presencias
pred_occs <- predict(best_model, occs[, vars], type = "cloglog")

# Predicciones en background
pred_bg <- predict(best_model, bg[, vars], type = "cloglog")

# ROC para calcular umbrales
library(pROC)
all_preds <- c(pred_occs, pred_bg)
all_obs <- c(rep(1, length(pred_occs)), rep(0, length(pred_bg)))
roc_obj <- roc(all_obs, all_preds, quiet = TRUE)

# Threshold 10% presencia de entrenamiento (10p)
thresh_10p <- quantile(pred_occs, 0.1)
cat("\nUmbral 10% presencias:", thresh_10p, "\n")

# Maximum Training Presence (MTP)
thresh_mtp <- min(pred_occs)
cat("Umbral MTP:", thresh_mtp, "\n")

# Equal sensitivity-specificity (ESS)
coords_roc <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
thresh_ess <- coords_roc$threshold
cat("Umbral ESS (mejor sensibilidad-especificidad):", thresh_ess, "\n")

# Guardar umbrales
umbrales <- data.frame(
  Metodo = c("10% Training Presence", "MTP", "Equal Sens-Spec"),
  Umbral = c(thresh_10p, thresh_mtp, thresh_ess)
)
write.csv(umbrales, "umbrales_cvenustus.csv", row.names = FALSE)