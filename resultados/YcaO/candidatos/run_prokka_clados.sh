#!/bin/bash

#To run in the terminal:
#chmod +x run_prokka_clados.sh
#conda activate prokka
#./clado_x/fastas_clado_x/run_prokka_clados.sh

#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Get the directory where this script is actually located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Move into that directory so all relative paths work correctly
cd "$SCRIPT_DIR"

# Create the main output directory called 'prokka' inside fastas_clado_A
mkdir -p prokka

# Loop through all .fna files in this specific directory
for file in *.fna; do
    # Ensure the file exists (prevents errors if no .fna files are found)
    [ -e "$file" ] || continue
    
    # Extract the base name without the .fna extension
    basename=$(basename "$file" .fna)
    
    echo "=========================================="
    echo "Processing genome: $basename"
    echo "=========================================="
    
    # Run Prokka, putting each result in its own subfolder inside prokka/
    prokka --outdir "prokka/$basename" --prefix "$basename" "$file"
    
    echo "Finished $basename. Output saved to prokka/$basename/"
done

echo "All Prokka annotations completed successfully!"