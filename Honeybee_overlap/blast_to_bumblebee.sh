#!/bin/bash

#PBS -q devel
#PBS -N blast
#PBS -l walltime=00:05:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load blast+/2.5.0

for file in $(ls *.fasta)
do
	base=$(basename $file ".fasta")
    blastx \
    -query ${file} \
    -db "nr" \
    -entrez_query "Bombus terrestris [organism]" \
    -remote \
    -max_target_seqs 1 \
    -outfmt 5 \
    -out ${base}.xml \
    -evalue 1e-3
done