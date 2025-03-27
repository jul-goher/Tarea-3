#Cargar paquetes
library(phyloseq)
library(microbiome)
library(ggplot2)
library(vegan)
library(dplyr)

#Cargar objeto phyloseq contenido en microbiome 
data("dietswap", package = "microbiome")
ps <- dietswap
ps

#¿Cuántas muestras y taxones contiene el objeto?

#¿Qué variables están disponibles en los metadatos de las muestras?