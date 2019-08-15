#!/bin/bash

#PBS -N alignments_RNA
#PBS -l walltime=06:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=2:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load star/2.5.2b
                                                                                      
# create the directory where the output files are to be written  
OUTPUT=alignments                                                                                                                                 
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

######################################################################

# NOTE: NEED TO RUN THIS SCRIPT TWICE; ONCE FOR ABDOMEN AND ONCE FOR HEAD SAMPLES IN CORRESPONDING DIRECTORIES

# Run STAR alignments to M02 genome
for file in $(ls trimmed_02*1.fq.gz)
do
    base=$(basename $file "1.fq.gz")
    STAR \
    --runThreadN 16 \
    --outFilterMismatchNmax 0 \
    --genomeDir /scratch/monoallelic/hm257/leuven_rna/alternate_genomes/m02_genome \
    --readFilesIn ${base}1.fastq ${base}2.fastq \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix M_${base}
done
wait

for file in $(ls *1.fastq)
do
    base=$(basename $file "1.fastq")
    STAR \
    --runThreadN 16 \
    --genomeDir /scratch/monoallelic/hm257/genome \
    --readFilesIn ${base}1.fastq ${base}2.fastq \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${base}
done

# Run STAR alignments to Q02 genome
for file in $(ls trimmed_02*1.fq.gz)
do
    base=$(basename $file "1.fq.gz")
    STAR \
    --runThreadN 16 \
    --outFilterMismatchNmax 0 \
    --genomeDir /scratch/monoallelic/hm257/rna_leuven/rna_seq_genomes/alternate_indices/q02_index \
    --readFilesCommand gunzip -c \
    --readFilesIn ${base}1.fq.gz ${base}2.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix Q_${base}
done
wait

######################################################################

# Run STAR alignments to M12 genome
for file in $(ls trimmed_12*1.fq.gz)
do
    base=$(basename $file "1.fq.gz")
    STAR \
    --runThreadN 16 \
    --outFilterMismatchNmax 0 \
    --genomeDir /scratch/monoallelic/hm257/rna_leuven/rna_seq_genomes/alternate_indices/m12_index \
    --readFilesCommand gunzip -c \
    --readFilesIn ${base}1.fq.gz ${base}2.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix M_${base}
done
wait

# Run STAR alignments to Q12 genome
for file in $(ls trimmed_12*1.fq.gz)
do
    base=$(basename $file "1.fq.gz")
    STAR \
    --runThreadN 16 \
    --outFilterMismatchNmax 0 \
    --genomeDir /scratch/monoallelic/hm257/rna_leuven/rna_seq_genomes/alternate_indices/q12_index \
    --readFilesCommand gunzip -c \
    --readFilesIn ${base}1.fq.gz ${base}2.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix Q_${base}
done
wait

######################################################################

# Run STAR alignments to M22 genome
for file in $(ls trimmed_22*1.fq.gz)
do
    base=$(basename $file "1.fq.gz")
    STAR \
    --runThreadN 16 \
    --outFilterMismatchNmax 0 \
    --genomeDir /scratch/monoallelic/hm257/rna_leuven/rna_seq_genomes/alternate_indices/m22_index \
    --readFilesCommand gunzip -c \
    --readFilesIn ${base}1.fq.gz ${base}2.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix M_${base}
done
wait

# Run STAR alignments to Q22 genome
for file in $(ls trimmed_22*1.fq.gz)
do
    base=$(basename $file "1.fq.gz")
    STAR \
    --runThreadN 16 \
    --outFilterMismatchNmax 0 \
    --genomeDir /scratch/monoallelic/hm257/rna_leuven/rna_seq_genomes/alternate_indices/q22_index \
    --readFilesCommand gunzip -c \
    --readFilesIn ${base}1.fq.gz ${base}2.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix Q_${base}
done
wait

######################################################################

# Run STAR alignments to M31 genome
for file in $(ls trimmed_31*1.fq.gz)
do
    base=$(basename $file "1.fq.gz")
    STAR \
    --runThreadN 16 \
    --outFilterMismatchNmax 0 \
    --genomeDir /scratch/monoallelic/hm257/rna_leuven/rna_seq_genomes/alternate_indices/m31_index \
    --readFilesCommand gunzip -c \
    --readFilesIn ${base}1.fq.gz ${base}2.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix M_${file}
done
wait

# Run STAR alignments to Q31 genome
for file in $(ls trimmed_31*1.fq.gz)
do
    base=$(basename $file "1.fq.gz")
    STAR \
    --runThreadN 16 \
    --outFilterMismatchNmax 0 \
    --genomeDir /scratch/monoallelic/hm257/rna_leuven/rna_seq_genomes/alternate_indices/q31_index \
    --readFilesCommand gunzip -c \
    --readFilesIn ${file}1.fq.gz ${file}2.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix Q_${file}
done

# Only need to keep .bam files, .Log.final.out and SJ (splice junction) files. 
