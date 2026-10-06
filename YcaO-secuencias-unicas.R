# Este es un código para eliminas las secuencias repetidas del fasta a partir
# del blastp realizado a las especies que se descargaron de InterPro (se 
# encuentran en el archivo `raw-data/YcaO-selection/blast-YcaO`)

# Despues de obtener los datos limpios de hace el alineamiento con 
# MSA_arboles-filo.R

# Librerías requeridas
library(Biostrings)

##############################################
############# SECUENCIAS FASTA ###############
##############################################

# Definir la ruta de los fasta para evitar cargarlos uno por uno 
archivos_fasta <- list.files(path = "raw-data/YcaO-selection/blast-YcaO/", pattern = "*.txt", full.names = TRUE)

# Leer y combinar todas las recuencias de los archivos en un solo objeto
todas_secuencias <- AAStringSet()

for (archivo in archivos_fasta) {
  seqs <- readAAStringSet(archivo)
  todas_secuencias <- c(todas_secuencias, seqs)
}

cat("Total de secuencias leídas (incluyendo duplicados):", length(todas_secuencias), "\n")

# Eliminar las secuencias duplicadas
# unique() nos permite comparar la cadena de a.a. de cada entrada
secuencias_unicas <- unique(todas_secuencias)

cat("Total de secuencias únicas conservadas:", length(secuencias_unicas), "\n")

# Guardar el nuevo archivo en formato fasta para el siguiente análisis (MSA y arbol filo)
writeXStringSet(secuencias_unicas, filepath = "raw-data/YcaO-selection/blast-YcaO/secuencias-unicas-YcaO.fasta")
