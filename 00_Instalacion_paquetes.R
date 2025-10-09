##############################################
#  Verificación e instalación de paquetes
#  Curso: Introducción a los Modelos Predictivos de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################

# Lista COMPLETA de paquetes necesarios para todos los scripts
paquetes_necesarios <- c(
  # Para Script 01 (Limpieza de ocurrencias):
  "CoordinateCleaner", "spThin", "dplyr", "sf", "mapview",
  
  # Para Script 02 (Descarga de capas):
  "terra", "utils",
  
  # Para Script 03 (Análisis PCA):
  "ggplot2", "factoextra"
)

# Función para verificar e instalar paquetes
verificar_paquetes <- function(paquetes) {
  for (paquete in paquetes) {
    if (!requireNamespace(paquete, quietly = TRUE)) {
      cat("Instalando paquete:", paquete, "\n")
      install.packages(paquete, dependencies = TRUE)
    } else {
      cat("Paquete", paquete, "ya está instalado\n")
    }
  }
}

# Verificar e instalar paquetes faltantes
cat("=== VERIFICACIÓN DE PAQUETES ===\n")
verificar_paquetes(paquetes_necesarios)

# Cargar todos los paquetes
cat("\n=== CARGANDO PAQUETES ===\n")
for (paquete in paquetes_necesarios) {
  if (require(paquete, character.only = TRUE, quietly = TRUE)) {
    cat("✓", paquete, "cargado correctamente\n")
  } else {
    cat("✗ Error al cargar:", paquete, "\n")
  }
}

# Verificación final
cat("\n=== VERIFICACIÓN FINAL ===\n")
paquetes_cargados <- sapply(paquetes_necesarios, 
                            function(x) requireNamespace(x, quietly = TRUE))
if (all(paquetes_cargados)) {
  cat("✅ TODOS los paquetes están instalados y disponibles\n")
} else {
  cat("❌ Los siguientes paquetes tienen problemas:\n")
  print(paquetes_necesarios[!paquetes_cargados])
}