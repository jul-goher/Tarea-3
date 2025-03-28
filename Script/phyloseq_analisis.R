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
#Cargar paquetes
library(ggplot2)
library(data.table)

#Eliminar taxones menores a 1
ps_sin <- prune_taxa (taxa_sums(ps) > 1, ps) 

#Verificar los taxones con menos lecturas 
cuenta_lecturas <- data.table(as(sample_data(ps_sin), "data.frame"),
                        TotalReads = sample_sums(ps_sin), 
                        keep.rownames = TRUE) #Crea una tabla

setnames (cuenta_lecturas, "rn", "SampleID")

ggplot (cuenta_lecturas, aes(TotalReads)) + geom_histogram() + ggtitle("Cobertura Secuenciación")
#El histograma nos indica la distribución de las lecturas en las muestras
#Es posible ver que hay muestras que tienen muy pocas lecturas

#Tabla que ordene los conteos de lectura más bajos
head(cuenta_lecturas[order(cuenta_lecturas$TotalReads), c("SampleID", "TotalReads")])
#Con esta tabla, vemos que la muestra 56 tiene la menor cantidad de lecturas con 1,776 lecturas

## Curva de rarefacción 
otu_cr <- otu_table (ps_sin) #extrar la tabla de abundancias de los OTUs del objeto phyloseq
otu_cr <- as.data.frame (t(otu_cr)) #Cambiar a un data frame
sample_names <- rownames (otu_cr) #Añadir nombres según los renglones del data-frame

# Step es un intervalo en el que se calcula la riqueza conforme añade lecturas
# Step define cuantas lecturas se añaden en incrementos x al realizar la curva
# En este caso elegí 200 ya que computaba más rápido que a 100

otu.rarecurve = rarecurve (otu_cr, step = 200, label=FALSE)
#Eje-y -> riqueza de las especies
#Eje-x -> cobertura 

#Guardar la curva
pdf("Figuras/curva_rarefaccion.pdf", width = 13, height = 8)
otu.rarecurve = rarecurve (otu_cr, step = 200)
dev.off()

#Apoyo: https://search.r-project.org/CRAN/refmans/vegan/html/rarefy.html
#
                  ################
#####             Diversidad alpha               #####
                  ################
#Calcula y grafica:
# Riqueza observada + ínidce de Shannon & Simpson
plot_richness (ps_sin, x = "nationality", measures = c("Observed", "Shannon", "Simpson") ) + 
              geom_boxplot() + ggtitle ("Diversidad Alpha") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 7) ) 


#plot_richness(ps_sin, x = "bmi_group", color = "nationality" , measures = c("Observed", "Shannon", "Simpson"))
#plot_richness(ps_sin, x = "nationality", color = "sex" , measures = c("Observed", "Shannon", "Simpson"))




                #########################
#####           Filtrado y Transformación         #####
                #########################
#géneros más abundantes 
#Ej., los que tienen más del 0.1% de abundancia relativa en al menos 10% de las muestras




                #########################
#####                Diversidad beta              #####
                #########################
#Ordención PCoA utilizando distancia Bray-Curtis

data (esophagus)
distance(esophagus, "bray") 


ordinate() 
plot_ordination()


my.physeq <- import("Biom", BIOMfilename="myBiomFile.biom")
my.ord    <- ordinate(my.physeq)
plot_ordination(my.physeq, my.ord, color="myFavoriteVarible")


GloPa.pcoa = ordinate(GlobalPatterns, method="PCoA", distance=GPUF)

            #########################
#####         Gráfica Rank-abundance              #####
            #########################









