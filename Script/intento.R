# Código nuevo 

# FECES
gp_feces <- subset_samples (gp_subset, SampleType == "Feces") #utilizo subset generado previamente 
# Suma de las abundancias 
fe_abund_total <- taxa_sums(gp_feces)
# Acomodo de mayor a menor
fe_abundancias_ordenadas <- sort(fe_abund_total, decreasing = TRUE)
# Secuencias para cada taxón 
fe_secuencia_taxones <- seq_along(fe_abundancias_ordenadas)

##Generar el data frame a partir de los objetos generados 
fe_rank_df <- data.frame (
  Abundancia = fe_abundancias_ordenadas,
  Especies = fe_secuencia_taxones )

##Generar rank-abundance para heces 
fe_ra <- ggplot (fe_rank_df, aes (x = Especies, y = Abundancia, colour = Especies)) +
  geom_point() + geom_line() +  
  scale_y_log10() +
  labs (title = "Rank-Abundance - Heces ", x = "Spp", y = "Abundancia")

print (fe_ra) 


# SUELO 
gp_soil <- subset_samples (gp_subset, SampleType == "Soil") #utilizo subset generado previamente 
# Suma de las abundancias 
soil_abund_total <- taxa_sums(gp_soil)
# Acomodo de mayor a menor
soil_abundancias_ordenadas <- sort(soil_abund_total, decreasing = TRUE)
# Secuencias para cada taxón 
soil_secuencia_taxones <- seq_along(soil_abundancias_ordenadas)

##Generar el data frame a partir de los objetos generados 
soil_rank_df <- data.frame (
  Abundancia = soil_abundancias_ordenadas,
  Especies = soil_secuencia_taxones )

##Generar rank-abundance para suelo 
soil_ra <- ggplot (soil_rank_df, aes (x = Especies, y = Abundancia, colour = Especies)) +
  geom_point() + geom_line() +  
  scale_y_log10() +
  labs (title = "Rank-Abundance - Suelo ", x = "Spp", y = "Abundancia")

print (soil_ra) 



# SKIN
gp_skin <- subset_samples (gp_subset, SampleType == "Skin") #utilizo subset generado previamente 
# Suma de las abundancias 
skin_abund_total <- taxa_sums(gp_skin)
# Acomodo de mayor a menor
skin_abundancias_ordenadas <- sort(skin_abund_total, decreasing = TRUE)
# Secuencias para cada taxón 
skin_secuencia_taxones <- seq_along(skin_abundancias_ordenadas)

##Generar el data frame a partir de los objetos generados 
skin_rank_df <- data.frame (
  Abundancia = skin_abundancias_ordenadas,
  Especies = skin_secuencia_taxones )

##Generar rank-abundance para piel  
skin_ra <- ggplot (skin_rank_df, aes (x = Especies, y = Abundancia, colour = Especies)) +
  geom_point() + geom_line() +  
  scale_y_log10() +
  labs (title = "Rank-Abundance - Skin ", x = "Spp", y = "Abundancia")
  
print (skin_ra) 



