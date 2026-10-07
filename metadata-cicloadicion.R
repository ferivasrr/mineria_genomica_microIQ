# Este código es para una parte de la metadata del analisis de la enzima 
# [4+2] cicloadicion con los archivos fasta, el obtener: WP..., enzima casificada
# y especie

#######################################################
# PRIMERA PARTE: JUNTAR LOS FASTA Y ELIMINAR REPETIDOS 

# Librerías requeridas
library(Biostrings)

# Definir la ruta de los fasta para evitar cargarlos uno por uno
archivos_fasta_cicloacion <- list.files(path = "resultados/cicloadicion/", pattern = "*txt", full.names = TRUE)

# Leer y combinar todas las recuencias de los archivos en un solo objeto
todas_secuencias_cicloadicion <- AAStringSet()

for (archivo in archivos_fasta_cicloacion) {
  seqs <- readAAStringSet(archivo)
  todas_secuencias_cicloadicion <- c(todas_secuencias_cicloadicion, seqs)
}

cat("Total de secuencias leídas (incluyendo duplicados):", length(todas_secuencias_cicloadicion), "\n")

# Eliminar las secuencias duplicadas
# unique() nos permite comparar la cadena de a.a. de cada entrada
secuencias_unicas_cicloadicion <- unique(todas_secuencias_cicloadicion)

cat("Total de secuencias únicas conservadas:", length(secuencias_unicas_cicloadicion), "\n")

# Guardar el nuevo archivo en formato fasta para el siguiente análisis (MSA y arbol filo)
writeXStringSet(secuencias_unicas_cicloadicion, filepath = "resultados/cicloadicion/secuencias_cicloadicion_metadata.fasta")

#######################################################
# SEGUNDA PARTE: EXTRACCIÓN DE LA PRIMERA FILA FASTA
 
