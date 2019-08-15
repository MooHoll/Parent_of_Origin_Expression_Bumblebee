#!/bin/bash

#PBS -N intersect_bed
#PBS -l walltime=03:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=2:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bedtools/2.25.0
                                                                                      
# create the directory where the output files are to be written  
OUTPUT=vcf_files_with_counts                                                                                                                                
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

##############################

#NOTE: Need to run this script in the directory for the head and the abdomen samples 

# Take all the informative snp positions within the vcf file and count the number of reads that map to that position
# Done here on a colony basis, can always filter later on SNPs that occur in all 4 males/queens from the crosses
for file in $(ls M_trimmed_02*)
do
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a m02_unique.vcf -b ${file} -c -header -sorted > ${base}male.vcf
done


for file in $(ls Q_trimmed_02*)
do 
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a q02_unique.vcf -b ${file} -c -header -sorted > ${base}queen.vcf
done


for file in $(ls M_trimmed_12*)
do
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a m12_unique.vcf -b ${file} -c -header -sorted > ${base}male.vcf
done



for file in $(ls Q_trimmed_12*)
do 
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a q12_unique.vcf -b ${file} -c -header -sorted > ${base}queen.vcf
done


for file in $(ls M_trimmed_22*)
do
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a m22_unique.vcf -b ${file} -c -header -sorted > ${base}male.vcf
done



for file in $(ls Q_trimmed_22*)
do 
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a q22_unique.vcf -b ${file} -c -header -sorted > ${base}queen.vcf
done


for file in $(ls M_trimmed_31*)
do
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a m31_unique.vcf -b ${file} -c -header -sorted > ${base}male.vcf
done



for file in $(ls Q_trimmed_31*)
do 
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a q31_unique.vcf -b ${file} -c -header -sorted > ${base}queen.vcf
done


### Now do it for the opposite queen/male snp file

for file in $(ls Q_trimmed_02*.bam)
do
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a m02_unique.vcf -b ${file} -c -header -sorted > ${base}male.vcf
done



for file in $(ls M_trimmed_02*.bam)
do 
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a q02_unique.vcf -b ${file} -c -header -sorted > ${base}queen.vcf
done


for file in $(ls Q_trimmed_12*.bam)
do
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a m12_unique.vcf -b ${file} -c -header -sorted > ${base}male.vcf
done



for file in $(ls M_trimmed_12*.bam)
do 
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a q12_unique.vcf -b ${file} -c -header -sorted > ${base}queen.vcf
done


for file in $(ls Q_trimmed_22*.bam)
do
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a m22_unique.vcf -b ${file} -c -header -sorted > ${base}male.vcf
done



for file in $(ls M_trimmed_22*.bam)
do 
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a q22_unique.vcf -b ${file} -c -header -sorted > ${base}queen.vcf
done


for file in $(ls Q_trimmed_31*.bam)
do
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a m31_unique.vcf -b ${file} -c -header -sorted > ${base}male.vcf
done



for file in $(ls M_trimmed_31*.bam)
do 
	base=$(basename $file "Aligned.sortedByCoord.out.bam")
	bedtools intersect -wa -a q31_unique.vcf -b ${file} -c -header -sorted > ${base}queen.vcf
done
wait

#### Make files that contain just the scaffold, SNP and count
for file in $(ls *.vcf)
do
  	base=$(basename ${file} ".vcf")
        grep -v "#" ${file} | cut -f1,2,11 > ${base}_counts_only.txt
done

