#!/bin/bash

# Run CellBender
module load python/3.8.x-anaconda
module load gcc/12.2.0
source activate cellbender

cd /home2/gkonop/project/02_MATRIX_FOR_CELLBENDER/KRIENEN_MOUSE

for dir in `ls -d $PWD/SRR* | grep -v SRR11921005`
do

cellbender remove-background --input $dir --output $dir"/CellBender_out.h5" --cuda

done
