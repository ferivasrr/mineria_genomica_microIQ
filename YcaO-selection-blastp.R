# Codigo que sirve para hacer el filtrado de los mejores candidatos (secuencias)
# de organismos que tienen la YcaO y con ellos poder realizar mas analisis

# YcaO obtenida de INTERPRO, con ID: IPR003776 y número Pfam: PF02624
# ID obtenido a traves de la revision del articulo --> DOI: 10.1039/D0NP00027B

# =============================== PAQUETES ===================================
install.packages("tidyverse")
BiocManager::install("Biostrings")

library(tidyverse)
library(Biostrings)

# =========================== SUBIR DOCUMENTOS ===============================
tabla_interpro_YcaO <- read.csv(file = "raw-data/YcaO-selection/YcaO-all-PF02624.tsv",
sep = "\t")

# ========================== FILTRADO SECUENCIAS =============================
# Identificar la distribución de las longitudes
summary(tabla_interpro_YcaO$Length)

# Filtrado de secuencias 
candidatos_filtrados <- tabla_interpro_YcaO %>%
  filter(Length >= 512 & Length <= 647) %>% # filtrado de secuencias por tamaño
  group_by(Tax.Name) %>% # agrupar por especies para evitar cepas redundantes
  arrange(desc(Length)) %>% # dentro de cada especie elegir la más representativa
  slice(1) %>%
  ungroup()

tabla_candidatos_filtrados <- as.data.frame(candidatos_filtrados)
# Seleccion de 3-5 candidatos representativos 
set.seed(42)
candidatos_finales <- candidatos_filtrados %>%
  slice_sample(n = min(3, nrow(candidatos_filtrados)))

# Mostrar los candidatos seleccionados
print(candidatos_finales %>% select(Accession, Tax.Name, Length))

# ========================== FILTRADO FASTA =============================
# Extraer del documento .fasta los candidatos que se usaran para el blastp
# Cargar el archivo
tabla_fasta_YcaO <- readAAStringSet("raw-data/YcaO-selection/YcaO-all-PF02624.fasta")

# Obtener el listado de IDs
ids_candidatos <- candidatos_finales$Accession

# Filtrar los IDs que coincidan sólo con las secuencias
# se hace este filtrado ya que en el formato FASTA esta todo junto 
busqueda <- paste(ids_candidatos, collapse = "|")
indices_coinciden <- grep(busqueda, names(tabla_fasta_YcaO))

secuencias_filtradas_fasta <- tabla_fasta_YcaO[indices_coinciden]

# Ver las secuencias que se encontraron 
cat("Se encontraron", length(secuencias_filtradas_fasta), "secuencias de tus candidatos. \n")
print(secuencias_filtradas_fasta)

# Descargar el archivo con los candidatos seleccionados
writeXStringSet(secuencias_filtradas_fasta, "raw-data/YcaO-selection/YcaO-candidates-PF02624.fasta")
