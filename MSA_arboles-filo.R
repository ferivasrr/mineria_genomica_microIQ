# Este codigo es para el alineamiento multiple y generacion de arboles filogeneticos 
# de las secuencias fasta obtenidas de YcaO y [4+2] cicloadición a partir de INTERPRO

#Paquetes que se necesitan:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install.packages("seqinr") #manejar secuencias 
BiocManager::install("msa") #Hacer el alineamiento
BiocManager::install("ggmsa") #Visualizacion del alineamiento
BiocManager::install("Biostrings") #Manejar secuencias
install.packages("ape") #Realizar el arbol filogenetico
BiocManager::install("ggtree") #Visualizar el arbol filogenetico
BiocManager::install("DECIPHER") #Hacer alineamiento y curacion 

#Cargar los paquetes instalados
library(Biostrings)
library(msa)
library(seqinr)
library(ggmsa)
library(ape)
library(ggtree)
library(DECIPHER)

#Leer las secuencias fasta

#####################################
############# GE2270A ###############
#####################################

GE2270A_tpdB <- readAAStringSet("raw-data/fasta/fasta-parametros_amplios/GE2270A-tpdB.txt")
GE2270A_tpdB #para visualizarlo 
summary(width(GE2270A_tpdB)) #ver los datos generales de la longitud de secuencias 

# COMPLETARLO CON LAS DEMÁS 

#####################################
###### ALINEAMIENTOS MÚLTIPLES ######
#####################################

?msa() #ver los parámetros que puedes usar 
msa_GE2270A_tpdB <- msa(GE2270A_tpdB, method = 'Muscle') #alineamiento

#visualizacion y descarga del msa
msa_GE2270A_tpdB_compatible <- unmasked(msa_GE2270A_tpdB)
BrowseSeqs(msa_GE2270A_tpdB_compatible)

writeXStringSet(msa_GE2270A_tpdB_compatible, "resultados/alineamientos-crudos/msa_GE2270A_tpdB_original.fasta")

#### Para clipkit en la terminal antes de cargarlo aquí para el arbol
# conda activate MSA_conda
# clipkit ruta/al/archivo -o /ruta/nuevo/archivo/y/su/nuevo/nombre

#############################################
####### ALINEAMIENTOS CURADOS CLIPKIT #######
#############################################

msa_GE2270A_tpdB_curada <- readAAMultipleAlignment("resultados/alineamientos-crudos/msa_GE2270A_tpdB_curadas_clipkit.fasta")
msa_GE2270A_tpdB_curada_compatible <- unmasked(msa_GE2270A_tpdB_curada)
BrowseSeqs(msa_GE2270A_tpdB_curada_compatible)

###########################################################
########### MATRICES DE DISTANCIA PARA ARBOL NJ ###########
###########################################################

#Cambiar la clase del objeto con el que vamos a trabajar, de AAMultipleAlignment/Biostrings a Alignment
tree_msa_GE2270A_tpdB_curada <- msaConvert(msa_GE2270A_tpdB_curada, type = "seqinr::alignment")

#Matriz de distancia para NeightborJoining
matriz_distancia_msa_GE2270A <- dist.alignment(tree_msa_GE2270A_tpdB_curada, matrix = "identity")

#################################
###########  ARBOL NJ ###########
#################################

library(ggtree)
library(ape)
library(ggplot2)

nj_tree_GE2270A <- nj(matriz_distancia_msa_GE2270A) 

p <- ggtree(nj_tree_GE2270A) +
  # Ajuste del tamaño de las letras
  # Reduce font size and adjust positioning of tip labels
  geom_tiplab(size = 1.5, color = "darkblue", offset = 0.005) + 
  theme_tree2() +
  # Add margin space on the right so tip text doesn't clip
  xlim(0, 0.6)

ggsave("phylogeny_tree.pdf", plot = p, width = 30, height = 30, limitsize = FALSE)

# Otra forma de hacer un árbol 
tree <- TreeLine(msa_GE2270A_tpdB_curada_compatible,
method = "ML",
model = MODELS,
showPlot = TRUE,
processors = NULL)
