# Este código es para una parte de la metadata del analisis de la enzima 
# YcaO con los archivos fasta, el obtener: WP..., enzima casificada
# y especie

# Ya no se necesita hacer el filtrado de los repetidos, ya se realizó
# para poder hacer el alineamiento y la curación del alineamiento con el 
# código YcaO-secuencias-unicas.R

#######################################################
# Librerías requeridas
library(Biostrings)
library(dplyr)
library(tidyverse)
library(stringr)

#######################################################
# EXTRACCIÓN DE LA PRIMERA FILA FASTA

# 1. Leer el FASTA y obtener los encabezados completos (los rownames)
metadatos_YcaO <- readAAStringSet("resultados/YcaO/clipkit_secuencias_unicas_YcaO.fasta")
cabeceras <- names(metadatos_YcaO) # Usar names() en lugar de rownames() es más directo para Biostrings

# 2. Extraer WP, Nombre de la proteína y Especie
# Estructura típica de un FASTA NCBI: WP_XXXXX.1 protein name [Organism name]
df_metadatos_YcaO <- data.frame(cabecera = cabeceras) %>%
  extract(
    col = cabecera,
    into = c("WP_ID", "Proteina", "Especie"),
    regex = "^(\\S+)\\s+(.*?)\\s*\\[(.*)\\]$",
    remove = FALSE
  )

# 3. Separar de la especie la sepa para aquellas que sí aplique
df_metadatos_YcaO <- df_metadatos_YcaO %>%
  mutate(
    # Capturar los códigos de la cepa
    Cepa = str_extract(Especie, "(?<=sp\\.\\s|strain\\s)[A-Za-z0-9_-]+.*$"),
    # Limpia la especie quitando el código de cepa
    Especie = ifelse(!is.na(Cepa), str_remove(Especie, paste0("\\s*", str_escape(Cepa))), Especie)
  ) %>%
  select(WP_ID, Proteina, Especie, Cepa)

# Descarga 
write.csv(df_metadatos_YcaO, file = "resultados/YcaO/metadata_YcaO_parte1.csv")
