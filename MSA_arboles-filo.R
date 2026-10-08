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

