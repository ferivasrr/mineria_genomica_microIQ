# Definir la ruta de tu carpeta BAGEL4
BAGEL4_PATH="/Users/fernandarivasromero/Library/CloudStorage/OneDrive-Personal/Documentos/Maestria/UNAM/Documentos/avances/mineria_genomica_microIQ/BAGEL4"

# 1. Crear la carpeta principal de resultados
mkdir -p clado_A/fastas_clado_A/bagel4

# 2. Iterar sobre cada archivo .fna en tu carpeta de origen
for file in clado_A/fastas_clado_A/*.fna; do
    # Extraer solo el nombre del archivo sin la ruta ni la extensión
    filename=$(basename "$file" .fna)
    
    echo "============================================="
    echo "Procesando genoma (con anotación PFAM): $filename"
    echo "============================================="
    
    # Definir la carpeta de sesión para este genoma
    SESSION_DIR="clado_A/fastas_clado_A/bagel4/$filename"
    
    # Crear la carpeta de sesión y la subcarpeta 'queryfolder' obligatoria que exige BAGEL4
    mkdir -p "$SESSION_DIR/queryfolder"
    
    # Copiar el archivo FASTA directamente dentro del queryfolder de la sesión
    cp "$file" "$SESSION_DIR/queryfolder/"
    
    # Ejecutar BAGEL4 apuntando al session_dir y al queryfolder interno
    perl -I"$BAGEL4_PATH/lib" \
         -I"/opt/homebrew/Caskroom/miniconda/base/envs/bagel4/lib/perl5/site_perl/5.22.0" \
         "$BAGEL4_PATH/bagel4_wrapper.pl" \
         -s "$SESSION_DIR" \
         -query "$SESSION_DIR/queryfolder" \
         -full
         
    echo "---------------------------------------------"
done

echo "¡Proceso completado con éxito para todos los genomas!"
