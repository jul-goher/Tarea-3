#Cargar paquetes
library(phyloseq)
library(microbiome)
library(ggplot2)
library(vegan)
library(dplyr)
            
            ###################################
#####       Características del objeto phyloseq         #####
            ###################################
#Cargar objeto phyloseq contenido en microbiome 
data("dietswap", package = "microbiome")
ps <- dietswap
ps

#¿Cuántas muestras y taxones contiene el objeto?
nsamples (ps) 
ntaxa(ps) 

#¿Qué variables están disponibles en los metadatos de las muestras?
sample_variables (ps)


                  #####################
#####             Curvas de rarefacción               #####
                  #####################
#Para hacer la curva, primero se debe de extraer la tabla de OTUs del objeto ps
otu_cr <- otu_table(ps)
otu_cr
#Cambiar a un data frame
otu_cr <- as.data.frame (t(otu_cr))
#Añadir nombres según los renglones del data-frame
sample_names <- rownames(otu_cr)
#Crear la curva de rarefacción
otu.rarecurve = rarecurve (otu_cr, step = 30)
#Dibujan una curva de rarefacción para cada renglón del input data. 
#Step es un intervalo en el que se calcula la riqueza conforme añade lecturas


### Selección del número correcto de step 
lect_totales <- rowSums(t(otu_cr)) #lecturas totales de cada muestra, que están en los renglones
#Para evitar seleccionar muestras donde hay 0 reads
#Tengo que calcular las lecturas totales:
summary(lect_totales)
#El resumen permite visualizar el mínimo de reads, que es 0
#El máximo de reads es 873314 según este resumen 

#Eliminar los 0 de la tabla de abundancias
# para saber cuantas de esas lecturas totales son igual a 0
sum(lect_totales == 0)
#Indica que 7 muestras tienen un valor de 0

ps_sin <- any (ps > 0) #esto está mal 


otu_cr <- otu_table(ps_filtered)
muest_totales <- rowSums(t(otu_cr))


#lect_totales
#min(lect_totales)




#Curva de rarefacción 2
pdf("Genomica Funcional/Tarea-3/curva_rarefaccion", width=13,height = 8)
dev.off()

#





#Apoyo: https://search.r-project.org/CRAN/refmans/vegan/html/rarefy.html

                  ################
#####             Diversidad alpha               #####
                  ################
#Calcula y grafica:
# Riqueza observada 

# ínidce de Shannon 

# índice de Simpson

                #########################
#####           Filtrado y Transformación         #####
                #########################
#géneros más abundantes 
#Ej., los que tienen más del 0.1% de abundancia relativa en al menos 10% de las muestras



#####             Diversidad beta               #####
#Ordención PCoA utilizando distancia Bray-Curtis
ordinate() 
plot_ordination()






