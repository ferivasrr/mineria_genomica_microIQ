# Este codigo es para el alineamiento multiple y generacion de arboles filogeneticos 
# de las secuencias fasta obtenidas de YcaO ya filtradas a partir de InterPro

#Paquetes que se necesitan:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install.packages("devtools")
install.packages("seqinr") #manejar secuencias 
BiocManager::install("msa") #Hacer el alineamiento
BiocManager::install("ggmsa") #Visualizacion del alineamiento
BiocManager::install("Biostrings") #Manejar secuencias
install.packages("ape") #Realizar el arbol filogenetico
BiocManager::install("ggtree") #Visualizar el arbol filogenetico
BiocManager::install("DECIPHER") #Hacer alineamiento y curacion 
BiocManager::install("microseq")
devtools::install_github('grunwaldlab/heattree')

#Cargar los paquetes instalados
library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)
library(microseq)
library(heattree)

#Leer las secuencias fasta

#############################################
############# SECUENCIA FASTA ###############
#############################################

secuencias_unicas_YcaO <- readAAStringSet("raw-data/YcaO-selection/blast-YcaO/secuencias-unicas-YcaO.fasta")
secuencias_unicas_YcaO #para visualizarlo 
summary(width(secuencias_unicas_YcaO)) #ver los datos generales de la longitud de secuencias 

# COMPLETARLO CON LAS DEMÁS 

#############################################
######## ALINEAMIENTO MÚLTIPLE (MSA) ########
#############################################

?msa() #ver los parámetros que puedes usar 
msa_secuencias_unicas_YcaO <- msa(secuencias_unicas_YcaO, method = 'Muscle') #alineamiento

# Visualizacion y descarga del MSA, se usa unmasked para hacer el archivo compatible
msa_secuencias_unicas_YcaO_compatible <- unmasked(msa_secuencias_unicas_YcaO)
BrowseSeqs(msa_secuencias_unicas_YcaO_compatible)

writeXStringSet(msa_secuencias_unicas_YcaO_compatible, "resultados/YcaO/msa_secuencias_unicas_YcaO.fasta")

#############################################
####### ALINEAMIENTOS CURADOS CLIPKIT #######
#############################################

# ANTES DE CARGAR EL ARCHIVO:
# Tienes que ir a la terminal y poner estos comandos
# conda activate MSA_conda
# clipkit ruta/al/archivo.fasta -o /ruta/nuevo/archivo/y/su/nuevo/nombre.fasta

# DESPUES DE LA CURACIÓN DE ALINEAMIENTO CON CLIPKIT
# Cargar el nuevo archivo curado 
clipkit_secuencias_curadas_YcaO <- readAAMultipleAlignment("resultados/YcaO/clipkit_secuencias_unicas_YcaO.fasta")
clipkit_secuencias_curadas_YcaO_compatible <- unmasked(clipkit_secuencias_curadas_YcaO)
BrowseSeqs(clipkit_secuencias_curadas_YcaO_compatible)

###########################################################
########### MATRICES DE DISTANCIA PARA ARBOL NJ ###########
###########################################################

# USA ESTE CÓGIDO SI TUS ÁRBOLES SON PEQUEÑOS
# Cambiar la clase del objeto con el que vamos a trabajar, de AAMultipleAlignment/Biostrings a Alignment
# tree_clipkit_sequnicas_YcaO <- msaConvert(clipkit_secuencias_curadas_YcaO, type = "seqinr::alignment")

#Matriz de distancia para NeightborJoining
# matriz_distancia_clipkit_sequnicas_YcaO <- dist.alignment(tree_clipkit_sequnicas_YcaO, matrix = "identity")

#################################
###########  ARBOL NJ ###########
#################################

# nj_tree_sequnicas_YcaO <- nj(matriz_distancia_clipkit_sequnicas_YcaO) 

# ggtree(nj_tree_sequnicas_YcaO, layout = "circular", size = 1) +
  #geom_tiplab(size = 2, aes(angle=angle)) +
  #geom_nodelab(geom = "label") +
  #hexpand(0.05)

####################################
########### ARBOL IQTREE ###########
####################################

YcaO_iqtree <- read.tree("resultados/YcaO/IQtree-arbol/clipkit_secuencias_unicas_YcaO.fasta.treefile")

heat_tree(tree = YcaO_iqtree,  layout = 'circular')

# Si tienes la metadata de tu árbol, corre este código
metadata_ncbi_completa_YcaO <- read.csv(file = "resultados/YcaO/metadata_ncbi_completa_YcaO.csv", header = TRUE)
heat_tree(tree = YcaO_iqtree, metadata = metadata_query, layout = 'circular')

# Señalar las secuencias query en el árbol 
secuencias_query <- readRDS(file = "raw-data/YcaO-selection/candidatos-seleccionados-YcaO.rsd")
secuencias_query

metadata_query <- metadata_ncbi_completa_YcaO |>
  mutate(query = if_else(metadata_ncbi_completa_YcaO$accession %in% c("MBN2909913.1", "WP_310192601.1", "WP_053602024.1", "WP_173617751.1", "WP_411789310.1"), 'query', 'no'))

####################################
############ CANDIDADOS ############
####################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Selección de candidatos del clado A
clado_a_candidatos <- c("Brevibacillus aydinogluensis", "Paenibacillus cymbidii", "Paenibacillus sp. J2R1T-4",
"Paenibacillus sp. P4_E11", "Paenibacillus sp. PL91", "Paenibacillus sp. Soil522", "Paenibacillus paridis", "Paenibacillus sp. P4_A7",
"Cohnella kolymensis", "Cohnella pontilimi", "Cohnella candida")

#Filtrar las filas de esos candidatos
genomas_filtrados <- metadata_ncbi_completa_YcaO %>%
  filter(scientific_name %in% clado_a_candidatos) %>%
  select(accession,protein_name,organism,strain,scientific_name,isolation_source,geo_loc_name,host)

write.csv(genomas_filtrados, file = "resultados/YcaO/candidatos/clado_A/candidatos_clado_A.csv", row.names = FALSE )
#Para extraer los genomas directo de la terminal 
ids_accession_clado_A <- genomas_filtrados$accession
writeLines(ids_accession_clado_A, "resultados/YcaO/candidatos/clado_A/ids_accession_clado_A.txt")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Selección de candidatos del clado B
clado_b_candidatos <- c("Paenibacillus sp. F024", "Paenibacillus thalictri", "Paenibacillus sp. MSJ-34", "Paenibacillus sp. GYB004",
"Paenibacillus sp. XY282", "Cohnella terricola")

genomas_filtrados_clado_B <- metadata_ncbi_completa_YcaO %>%
  filter(scientific_name %in% clado_b_candidatos) %>%
  select(accession,protein_name,organism,strain,scientific_name,isolation_source,geo_loc_name,host)

genomas_filtrados_clado_B <- genomas_filtrados_clado_B[c(1:7)[-6], ]

write.csv(genomas_filtrados_clado_B, file = "resultados/YcaO/candidatos/clado_B/candidatos_clado_B.csv", row.names = FALSE )
#Para extraer los genomas directo de la terminal 
ids_accession_clado_B <- genomas_filtrados_clado_B$accession
writeLines(ids_accession_clado_B, "resultados/YcaO/candidatos/clado_B/ids_accession_clado_B.txt")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Selección de candidatos del clado C
clado_c_candidatos <- c("Paenibacillus larvae subsp. larvae", "Paenibacillus periandrae", "Paenibacillus sp. P4_E11",
"Paenibacillus sp. F024", "Paenibacillus sp. Soil750", "Paenibacillus sp. J22R160", "Cohnella herbarum", "Cohnella soli",
"Cohnella sp. WQ 127256")

genomas_filtrados_clado_C <- metadata_ncbi_completa_YcaO %>%
  filter(scientific_name %in% clado_c_candidatos) %>%
  select(accession,protein_name,organism,strain,scientific_name,isolation_source,geo_loc_name,host)

genomas_filtrados_clado_C <- genomas_filtrados_clado_C[-c(1, 11), ]

write.csv(genomas_filtrados_clado_C, file = "resultados/YcaO/candidatos/clado_C/candidatos_clado_C.csv", row.names = FALSE )
#Para extraer los genomas directo de la terminal 
ids_accession_clado_C <- genomas_filtrados_clado_C$accession
writeLines(ids_accession_clado_C, "resultados/YcaO/candidatos/clado_C/ids_accession_clado_C.txt")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Selección de candidatos del clado D
clado_d_candidatos <- c("Bacillus sp. JJ1609", "Bacillus sp. JJ1122", "Rossellomorea oryzaecorticis", "Rossellomorea vietnamensis",
"Rossellomorea aquimaris")

genomas_filtrados_clado_D <- metadata_ncbi_completa_YcaO %>%
  filter(scientific_name %in% clado_d_candidatos) %>%
  select(accession,protein_name,organism,strain,scientific_name,isolation_source,geo_loc_name,host)

write.csv(genomas_filtrados_clado_D, file = "resultados/YcaO/candidatos/clado_D/candidatos_clado_D.csv", row.names = FALSE )
#Para extraer los genomas directo de la terminal 
ids_accession_clado_D <- genomas_filtrados_clado_D$accession
writeLines(ids_accession_clado_D, "resultados/YcaO/candidatos/clado_D/ids_accession_clado_D.txt")
