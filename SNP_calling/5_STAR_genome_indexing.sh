#!/bin/bash

#PBS -N STAR_genome_indexes
#PBS -l walltime=00:05:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8
#PBS -q devel

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load star/1.5.2b

# Define file paths
REF_FA=/scratch/monoallelic/hm257/genome/GCF_000214255.1_Bter_1.0_genomic.fasta
REF_GFF=/scratch/monoallelic/hm257/leuven_rna/alternate_references/ref_Bter_1.0_top_level.gff3
                                                                                      

# create the directory where the output files are to be written  
OUTPUT=alternate_indices                                                                                                                                   
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

# Need to adjust headers of all fasta files generated before star will take them as they don't match the .gff
sed 's/:1*//' m02_unique.vcf.fasta | sed 's/>[0-9][0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9] */>/' | sed 's/>[0-9] */>/' > m02_new.fasta
sed 's/:1*//' m12_unique.vcf.fasta | sed 's/>[0-9][0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9] */>/' | sed 's/>[0-9] */>/' > m12_new.fasta
sed 's/:1*//' m22_unique.vcf.fasta | sed 's/>[0-9][0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9] */>/' | sed 's/>[0-9] */>/' > m22_new.fasta
sed 's/:1*//' m31_unique.vcf.fasta | sed 's/>[0-9][0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9] */>/' | sed 's/>[0-9] */>/' > m31_new.fasta
sed 's/:1*//' q02_unique.vcf.fasta | sed 's/>[0-9][0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9] */>/' | sed 's/>[0-9] */>/' > q02_new.fasta
sed 's/:1*//' q12_unique.vcf.fasta | sed 's/>[0-9][0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9] */>/' | sed 's/>[0-9] */>/' > q12_new.fasta
sed 's/:1*//' q22_unique.vcf.fasta | sed 's/>[0-9][0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9] */>/' | sed 's/>[0-9] */>/' > q22_new.fasta
sed 's/:1*//' q31_unique.vcf.fasta | sed 's/>[0-9][0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9][0-9] */>/' | sed 's/>[0-9][0-9] */>/' | sed 's/>[0-9] */>/' > q31_new.fasta


# Create a list of the fasta files to be called                                                                                                            
FILES=*new.fasta

# Run STAR to index new reference genomes (NOTE: this doesn't work, you need a directory for each file for --genomeDir)
for file in $FILES
do
    STAR \
    --runMode genomeGenerate \
    --genomeFastaFiles ${file} \
    --sjdbGTFfile ${REF_GFF} \
    --sjdbGTFtagExonParentTranscript Parent \
    --runThreadN 8 \
    --genomeDir m02_genome ##HERE! change for each      
done

### NOTE: need to make new directories for each index of each genome to be called later by STAR during alignment.