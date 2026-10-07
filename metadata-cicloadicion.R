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

library(dplyr)
library(tidyverse)
library(stringr)

# 1. Leer el FASTA y obtener los encabezados completos (los rownames)
metadatos_ciclacion <- readAAStringSet("resultados/cicloadicion/secuencias_cicloadicion_metadata.fasta")
cabeceras <- names(metadatos_ciclacion) # Usar names() en lugar de rownames() es más directo para Biostrings

# 2. Extraer WP, Nombre de la proteína y Especie
# Estructura típica de un FASTA NCBI: WP_XXXXX.1 protein name [Organism name]
df_metadatos <- data.frame(cabecera = cabeceras) %>%
  extract(
    col = cabecera,
    into = c("WP_ID", "Proteina", "Especie"),
    regex = "^(\\S+)\\s+(.*?)\\s*\\[(.*)\\]$",
    remove = FALSE
  )

# 3. Separar de la especie la sepa para aquellas que sí aplique
df_metadatos <- df_metadatos %>%
  mutate(
    # Capturar los códigos de la cepa
    Cepa = str_extract(Especie, "(?<=sp\\.\\s|strain\\s)[A-Za-z0-9_-]+.*$"),
    # Limpia la especie quitando el código de cepa
    Especie = ifelse(!is.na(Cepa), str_remove(Especie, paste0("\\s*", str_escape(Cepa))), Especie)
  ) %>%
  select(WP_ID, Proteina, Especie, Cepa)
