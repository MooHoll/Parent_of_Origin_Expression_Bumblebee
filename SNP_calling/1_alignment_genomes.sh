#!/bin/bash

#PBS -N BWA_alignment_to_ref
#PBS -l walltime=16:45:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bwa/0.7.15
module load samtools/1.3.2

# Define file paths
REF_FA=/scratch/monoallelic/hm257/genome/GCF_000214255.1_Bter_1.0_genomic.fa

# Remember to index the ref genome 1st: bwa index <genome_file>


# As downloaded from SRA with -I flag (urgh) need to remove the names given to reads
# in order for BWA to process as paired reads
for file in $(ls *.fastq)
do
    base=$(basename $file ".fastq")
    sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" ${file} > ${base}_edited.fastq
done


# Align all fasta files to the reference, this needs 12:45hrs if n=1:ppn=8
for file in $(ls *1.fastq)
do
	base=$(basename $file "_1.fastq")
	bwa mem ${REF_FA} ${base}_1.fastq ${base}_2.fastq > ${base}.sam
done
                   
# samtools flatstat <file.sam> for mapping percentage                   
                                                                                                                                                                                                                                                                                                                                                                      
# Create a list of the sam files to be sorted                                                                                                            
sams=*.sam

# Sort the same files and convert to bam
for file in $sams
do
    samtools sort ${file} -O bam -o $file.sorted.bam
done


# Create a list of the bam files to have duplicates removed
bams=*sorted.bam

# Remove duplicates
for file in $bams
do
    samtools rmdup ${file} ${file}.final.bam
done



     