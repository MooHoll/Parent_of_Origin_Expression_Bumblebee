#!/bin/bash

#PBS -N trimmed_fastqc
#PBS -l walltime=04:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

cd $PBS_O_WORKDIR 

module load fastqc/0.11.5

OUTPUT=trimmed_fastqc                                                                                          

# create the directory where the output files are to be written                                                                                                                                     
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

# Create a list of the files to be called                                                                                                            
FILES=*.fq.gz

for file in $FILES
do
	fastqc -o $OUTPUT $file
done
