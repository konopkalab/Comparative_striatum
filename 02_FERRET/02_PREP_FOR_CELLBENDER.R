#!/bin/bash

# Run CellBender
module load python/3.8.x-anaconda
module load gcc/12.2.0
source activate cellbender
source ~/load_modules.sh

# input dir
INPUT_DIR="/home2/gkonop/project/02_MATRIX_FOR_CELLBENDER/KRIENEN_FERRET"

### cellbender
# Directory where cellbender results will be stored
DIR="/home2/gkonop/project/03_MATRIX_FROM_CELLBENDER/FERRET_KRIENEN"

# Array to store folder names
ar1=()

# Collect folder names in INPUT_DIR
for folder in "${INPUT_DIR}"/*; do
  if [ -d "${folder}" ]; then  # Check if it's a directory
    ar1+=("$(basename "${folder}")")
  fi
done


# Iterate over the selected samples and run cellbender
for SAMPLE in "${ar1[@]}"
do
  # Define the path to the input files
  FILE_DIR="${INPUT_DIR}/${SAMPLE}/snRNA_count_data_top80000/"
 
  # Print the input file being processed
  echo "Input file directory: $FILE_DIR"

  # Run cellbender command
  echo "Running cellbender for sample: $SAMPLE"
  cellbender remove-background --cuda --total-droplets-included 80000 --input "${FILE_DIR}" --output "${DIR}/${SAMPLE}_cellbender_output_top80000.h5"

  echo "Cellbender completed for sample: $SAMPLE"
done

echo "Script execution complete."
