#!/bin/bash

#To run in the terminal:
#chmod +x run_prokka_clados.sh
#conda activate prokka
#./clado_x/fastas_clado_x/run_prokka_clados.sh

#!/bin/bash

#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define your target directory containing the .fna files
TARGET_DIR="clado_B/fastas_clado_B"

# Check if the target directory actually exists
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Directory '$TARGET_DIR' not found!"
    exit 1
fi

# Create the output 'prokka' directory inside that target folder
mkdir -p "$TARGET_DIR/prokka"

echo "Looking for .fna files in: $TARGET_DIR"

# Loop through all .fna files in the target directory
for file in "$TARGET_DIR"/*.fna; do
    # Ensure the file exists (prevents literal '*.fna' loop if none found)
    [ -e "$file" ] || continue
    
    # Extract just the filename without path and without extension (e.g., "sample1")
    filename=$(basename "$file" .fna)
    
    echo "=========================================="
    echo "Processing genome: $filename"
    echo "=========================================="
    
    # Run Prokka, outputting to clado_A/fastas_clado_A/prokka/filename
    prokka --outdir "$TARGET_DIR/prokka/$filename" --prefix "$filename" "$file"
    
    echo "Finished $filename. Output saved to $TARGET_DIR/prokka/$filename/"
done

echo "All Prokka annotations completed successfully!"

# After the analysis, run this in the terminal to have a directory just with the 
# .gbk files to use it for rodeo
# you're in the currect directory fastas_clado_x with the prokka directory
# mkdir prokka_gbk
#find prokka -type f -name "*.gbk" -exec cp {} prokka_gbk/ \;