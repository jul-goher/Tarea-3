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

#Guardar la curva como pdf
pdf("Figuras/curva_rarefaccion.png", width = 13, height = 8)
otu.rarecurve = rarecurve (otu_cr, step = 200)
dev.off()
#Guardadr como pdf
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
div_alpha_plot <- plot_richness (ps_sin, x = "nationality", measures = c("Observed", "Shannon", "Simpson") ) + 
              geom_boxplot() + ggtitle ("Diversidad Alpha") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 7) ) 

print(div_alpha_plot)

plot_richness(ps_sin, x = "bmi_group", color = "nationality" , measures = c("Observed", "Shannon", "Simpson" )) +
  geom_boxplot()

#plot_richness(ps_sin, x = "nationality", color = "sex" , measures = c("Observed", "Shannon", "Simpson"))

# Guardar gráfica de diversidad alpha
pdf("Figuras/div_alpha.pdf", width = 13, height = 8)
div_alpha_plot
dev.off()

## Escribir la diversidad alpha en un csv 
# Estimar la riqueza con una función  
div_alpha <- estimate_richness (ps_sin, measures = c("Observed", "Shannon", "Simpson"))

# Del objeto con los índices:
div_alpha$nationality <- sample_data(ps_sin)$nationality
#se utiliza el sample_data que se utilizó para extraer los metadatos del objeto phyloseq cuando se creó la tabla de lecturas 
#Se elige $nationality porque sino devuelve todas las variables de las muestras 
  
# Escribir el csv
write.csv(div_alpha, file = "Data/div_alpha.csv", row.names = TRUE)
#Añadí row.names = TRUE porque sino no salen los nombres de las muestras 


                #########################
#####           Filtrado y Transformación         #####
                #########################
#Géneros más abundantes 
#Más del 0.1% de abundancia relativa en al menos 10% de las muestras

#Código tomado y adaptado de: https://joey711.github.io/phyloseq/preprocess.html
# y https://web.stanford.edu/class/bios221/labs/phyloseq/lab_phyloseq.html?utm

ps_relativ  = transform_sample_counts (ps_sin, function(x) x / sum(x) )
ps_filtr <- filter_taxa(ps_relativ, function(x) sum(x > 0.001) > (0.1 * length(x)), TRUE)
ps_filtr

                #########################
#####                Diversidad beta              #####
                #########################
#Ordención PCoA utilizando distancia Bray-Curtis

#data (ps_filtr) #no reconoce a ps_filter como un dataset por ser un objeto phyloseq 
#distance(ps_filtr, "bray") no me funciona para este caso

#Basado en: GloPa.pcoa = ordinate(GlobalPatterns, method="PCoA", distance=GPUF)

ps_pcoa = ordinate (ps_filtr, method="PCoA", distance = "bray")
ps_pcoa

#De (p12 <- plot_ordination(GlobalPatterns, GloPa.pcoa, "samples", color="SampleType") + 
#geom_point(size=5) + geom_path() + scale_colour_hue(guide = FALSE) )
#De: plot_ordination(my.physeq, my.ord, color="myFavoriteVarible")
p_ordin <- plot_ordination (ps_filtr, ps_pcoa, color="nationality")
p_ordin


pdf("Figuras/PCoA_dietswap.pdf", width = 13, height = 8)
p_ordin
dev.off()


            #########################
#####         Gráfica Rank-abundance              #####
            #########################

# Abundancia totalde cada uno de los taxones
abund_total <- taxa_sums(ps_sin) #uso de taxa_sums como cuando se eliminaron los taxones menores a 1
# Orden de mayor a menor
abundancias_ordenadas <- sort (abund_total, decreasing = TRUE)
abundancias_ordenadas
# Generar números para cada uno de los taxones,según la longitud de abund_total
secuencia_taxones <- seq_along (abund_total)
# Crear un data.frame del objecto phyloseq 
rank_df <- data.frame (
  Abundancia = abundancias_ordenadas, 
  Especies = secuencia_taxones
)


# Curva de rarefacción 
rank_abundance <- ggplot (rank_df, aes(x = Especies, y = Abundancia, colour = Especies)) +
  geom_point() + geom_line() +  #puntos y línea
  labs( title = "Rank-abundance", x = "Spp", y = "Abundance")

print(rank_abundance)


pdf("Figuras/rank_abundance.pdf", width = 13, height = 8)
rank_abundance
dev.off()


            ################################
#####       Gráfica apiladada de abundancias      #####
            ################################

#Gráficas apiladas de abundancia por taxón (género o phylum)
#Agrupa y grafica la composición de cada muestra como gráfica de barras apiladas.

phylum_apliada <- plot_bar(ps_sin, x = "nationality", y = "Abundance", fill = "Phylum") + 
                geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

pdf("Figuras/abundancias_apiladas.pdf", width = 13, height = 8)
phylum_apliada
dev.off()



      ############################################
#####               GlobalPatterns                ######
      ###########################################
data("GlobalPatterns")
gp <- GlobalPatterns

################
#     Preprocesamiento
#################

#Filtrar taxa con menos de 5 lecturas en al menos 20% de las muestras
gp_filtr <- filter_taxa(gp, function(x) sum(x > 5) > (0.2 * length(x)), TRUE)
gp_filtr

#Aglomerar a nivel de Familia

#Transformar a abundancias relativas (%)
gp_rel  = transform_sample_counts (gp_filtr, function(x) x / sum(x) )

#Subset para incluir solo muestras de: Soil, Feces, Skin

#Muestra las dimensiones del objeto resultante


################
#     Diversidad Alfa
#################

#Calcular 3 índices de diversidad alfa ( Shannon , Simpson , Observed )
#Crear boxplots comparativos de los índices entre tipos de muestra
#Realizar prueba estadística (Kruskal-Wallis) para diferencias entre grupos


################
#     Curvas de rango-abundancia
#################
#Crear gráficas de rango-abundancia para cada tipo de muestra Usar escala log10 en Y Comparar patrones
#entre ambientes


################
#     Perfil taxonómico
#################

#Crear gráfico apilado de abundancia a nivel de Phylum
#Mostrar solo los 5 phyla más abundantes
#Agrupar por tipo de muestra
#Usar facet_wrap para comparar ambientes


################
#     Diversidad Beta
#################
#Calcular distancia Bray-Curtis
#Realizar PCoA
#Visualizar con:
  #Colores por tipo de muestra
  #Elipses de confianza del 95%
  #Incluir stress plot
#Realizar PERMANOVA para diferencias entre grupos




