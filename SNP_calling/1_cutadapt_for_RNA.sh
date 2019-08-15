#!/bin/bash

#PBS -N trim_all
#PBS -l walltime=4:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load cutadapt/1.11

# create the directory where the output files are to be written 
OUTPUT=trimmed_data                                                                                                                                     
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

# Create a list of the files to be called                                                                                                            
FILES=*.fastq

for file in $FILES
do
	cutadapt -u 10 -o ${OUTPUT}/trimmed_${file} ${file}
done
