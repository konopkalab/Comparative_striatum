#!/bin/bash

# Run CellBender
module load python/3.8.x-anaconda
module load gcc/12.2.0
source activate cellbender

for dir in `ls -d $PWD/*`
do

cellbender remove-background --input $dir --output $dir"/CellBender_out.h5" --cuda

done
