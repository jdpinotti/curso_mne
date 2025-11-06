##############################################
#  Instalación de Paquetes adicionales
#  Curso: Modelado de Nicho y Distribución
#  Autor: Juan Diego Pinotti
##############################################


########## Para instalar biomod2 #########

install.packages("biomod2", dependencies = TRUE)
library(biomod2)



######## Para instalar ntbox ############

if (!require('devtools')) install.packages('devtools')
devtools::install_github('luismurao/ntbox')

# Cuando pide instalar actualizaciones de paquetes apretar enter para que saltee ese paso
