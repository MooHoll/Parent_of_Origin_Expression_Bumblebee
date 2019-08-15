#!/bin/bash

#PBS -N SNP_calling
#PBS -l walltime=00:05:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8
#PBS -q devel

# Run script in the working directory it was submitted in (6hrs)
cd $PBS_O_WORKDIR 

# Load software needed
module load freebayes/1.1.0

# Define file paths
REF_FILE=/scratch/monoallelic/hm257/genome/GCF_000214255.1_Bter_1.0_genomic.fa                                                                                         

# create the directory where the output files are to be written 
OUTPUT=vcf_files                                                                                                                                      
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

# Create a list of the bam files to be called                                                                                                            
queens=m*.bam
males=f*.bam

# Run freebayes for queen samples, min count 2 of alternative alleles, min 5 reads per SNP, ignore complex events, indels and mnps
for file in $queens
do
    freebayes \
        -f ${REF_FILE} \
        -C 2 \
        -! 5 \
        -u \
        -i \
        -X \
        -b ${file}\
        > ${OUTPUT}/${file}.freebayes.vcf.gz
done

# Run freebayes for male samples with ploidy of 1
for file in $males
do
    freebayes \
        -f ${REF_FILE} \
        -p 1 \
        -! 5 \
        -C 1 \
        -u \
        -i \
        -X \
        -b ${file} \
        > ${OUTPUT}/${file}.freebayes.vcf.gz
done

