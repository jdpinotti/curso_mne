# ============================================================
# Script: 04_corr_vars_princomp.r
# Objetivo: 
#   1. Extraer valores ambientales en presencias
#   2. Estandarizar variables
#   3. PCA con matriz de correlación (princomp)
#   4. Exportar correlaciones variable-ejes para selección manual
# ============================================================

# ---- Librerías ----
library(terra)    
library(dplyr)    
library(ggplot2)  
library(factoextra)

# ---- 1. Cargar datos ----
# Raster stack de variables ambientales (recortadas en 02_descarga_capas.r)
files <- list.files("corte_m", pattern = "\\.tif$", full.names = TRUE)
capas <- rast(files)
names(capas)
# Puntos de presencia (csv con columnas species, lon, lat)
occs <- read.csv("cvenustus_clean_thin1.csv")

# Extraer valores ambientales en los puntos
valores <- terra::extract(capas, occs[, c("longitude", "latitude")])
valores <- valores[,-1]  # quitar columna ID que agrega `extract`

# Guardar tabla con presencias + valores ambientales
dir.create("outputs", showWarnings = FALSE)
write.csv(cbind(occs, valores), "outputs/valores_presencias.csv", row.names = FALSE)

# ---- 2. Estandarizar variables ----
# (No es necesario usar scale() si usaremos cor = TRUE, pero lo dejamos por claridad)
valores_std <- scale(valores)

# ---- 3. PCA con princomp (matriz de correlaciones) ----
pca <- prcomp(valores, center = TRUE, scale. = TRUE)
loadings_df <- pca$rotation
write.csv(loadings_df, "outputs/pca_loadings_prcomp.csv")


# Resumen de varianza explicada
summary(pca)

# Scree plot
# Extraemos la varianza explicada manualmente
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 60))


# ---- 4. Biplot ----
# Biplot base
biplot(pca, cex = 0.7)

fviz_pca_biplot(pca,
                repel = TRUE,
                col.var = "darkred",
                col.ind = "gray50") +
  ggtitle("PCA de variables ambientales (princomp)") +
  theme_minimal()


# ---- 5. Exportar correlaciones variables-ejes ----
# "loadings" son los coeficientes que relacionan variables con componentes
loadings_df <- unclass(loadings(pca))
write.csv(loadings_df, "outputs/pca_loadings_princomp.csv")
