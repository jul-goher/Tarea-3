# Código nuevo 
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
fe_ra <- ggplot (fe_rank_df, aes (x = Especies, y = Abundancia)) +
  geom_point() + geom_line() +  
  labs (title = "Rank-Abundance - Heces ", x = "Spp", y = "Abundancia") +
  theme_minimal()

print (fe_ra)

plot (fe_ra)
