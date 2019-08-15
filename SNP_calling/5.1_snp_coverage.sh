#!/bin/bash

#PBS -N snp_coverage
#PBS -l walltime=10:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=2:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load perl/5.24.0 
module load samtools/1.3.2
module load bam-readcount/20161006

REF_FA=/scratch/monoallelic/hm257/rna_leuven/genome/GCF_000214255.1_Bter_1.0_genomic.fa

#################################

# This script first indexes all genome .bam files 

for file in $(ls *.bam)
do
	samtools index ${file}
done


# Then creates a text file needed by bam-readcount of all specific SNP positions (perl script from bam-readcount)
for file in $(ls *.vcf)
do
    base=$(basename $file "_unique.vcf")
    perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${file} > ${base}_snp_sites_list.txt
done


# Then counts up coverage depth (amongst other things) for all the SNPs for each .bam alignment
# Can then get the average coverage for the informative SNPs
# -w suppresses warnings, -q1 is min quality of 1 

bam-readcount -w 1 -q1 -l m02_snp_sites_list.txt -f ${REF_FA} m02.bam > m02_snp.readcount
bam-readcount -w 1 -q1 -l m12_snp_sites_list.txt -f ${REF_FA} m12.bam > m12_snp.readcount
bam-readcount -w 1 -q1 -l m22_snp_sites_list.txt -f ${REF_FA} m22.bam > m22_snp.readcount
bam-readcount -w 1 -q1 -l m31_snp_sites_list.txt -f ${REF_FA} m31.bam > m32_snp.readcount
bam-readcount -w 1 -q1 -l q02_snp_sites_list.txt -f ${REF_FA} q02.bam > q02_snp.readcount
bam-readcount -w 1 -q1 -l q12_snp_sites_list.txt -f ${REF_FA} q12.bam > q12_snp.readcount
bam-readcount -w 1 -q1 -l q22_snp_sites_list.txt -f ${REF_FA} q22.bam > q22_snp.readcount
bam-readcount -w 1 -q1 -l q31_snp_sites_list.txt -f ${REF_FA} q31.bam > q31_snp.readcount

# Then cut -f1,2,3,4 and read into R to get means of the 4th coverage row