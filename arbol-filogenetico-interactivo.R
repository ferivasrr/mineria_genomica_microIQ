# Código para generar un html y hacer interactivo el árbol filogenético

install.packages("htmltools") 
install.packages("htmlwidgets")

library(htmltools)
library(htmlwidgets)
library(heattree)
library(ggtree)
library(dplyr)

# Árbol filogenético interactivo
YcaO_iqtree <- read.tree("resultados/YcaO/IQtree-arbol/clipkit_secuencias_unicas_YcaO.fasta.treefile")
secuencias_query <- readRDS(file = "raw-data/YcaO-selection/candidatos-seleccionados-YcaO.rsd")
metadata_ncbi_completa_YcaO <- read.csv(file = "resultados/YcaO/metadata_ncbi_completa_YcaO.csv", header = TRUE)

metadata_query <- metadata_ncbi_completa_YcaO |>
  mutate(query = if_else(metadata_ncbi_completa_YcaO$accession %in% c("MBN2909913.1", "WP_310192601.1", "WP_053602024.1", "WP_173617751.1", "WP_411789310.1"), 'Secuencia concenso (query)', 'Secuencia obtenida de BLASTP'))

arbol_html <- heat_tree(YcaO_iqtree, metadata = metadata_query)

# Guardar el html:
saveWidget(arbol_html, "anotaciones-bitacora/arbol-filogenetico-interactivo.html", selfcontained = TRUE)
