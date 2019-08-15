#### -------------------------------------------------------------------------------------
#### Reciprocally Blasting the Honeybee and Bumblebee Genomes for Homologous Genes
#### -------------------------------------------------------------------------------------

# Already had a fasta file with every gene from Eamonn, 'final_GFF_bter_1.0.fa'

# Made an identical honeybee file, downloaded the CDS fasta file from
# http://hymenopteragenome.org/beebase/?q=download_sequences

# Used the following 'sed' commands to edit the headings to make them compatible with
# an .R script that pulls out the fasta data for a given gene list 

sed 's/.*GB/>GB/g' amel_OGSv3.2_cds.fa > fasta.fa
sed 's/-RA//g' fasta.fa > a.mellifera_genes_fasta.fa



# ON ALICE:

# Needed to remove duplicated sequences from Eamonn's bumblebee fasta file before the
# blast database command below would work, needed to request 20GB memory for this as an 
# interactive job, then it ran in seconds

module load bbmap/36.47
module load java/1.8 

bash dedupe.sh in=final_GFF_bter_1.0.fasta out=final_GFF_bter_1.0_dedup.fasta


# Turned the two above fasta files into blast databases using below code. This makes a ton
# of files so put these in their own folders
module load blast+/2.5.0
makeblastdb –in <input_fasta> –dbtype nucl –parse_seqids

# Run the reciprocal blasts to output an xml and a table, runs in a couple of mins when
# 20GB memory requested with interactive job
blastn \
-query final_GFF_bter_1.0.fasta \
-db /scratch/monoallelic/hm257/leuven_rna/honeybee_databse/a.mellifera_genes_fasta.fa \
-max_target_seqs 1 \
-outfmt 5 \
-out bumblebee_genome_to_hb.xml

blastn \
-query final_GFF_bter_1.0.fasta \
-db /scratch/monoallelic/hm257/leuven_rna/honeybee_databse/a.mellifera_genes_fasta.fa \
-max_target_seqs 1 \
-outfmt 6 \
-out bumblebee_genome_to_hb.txt

blastn \
-query a.mellifera_genes_fasta.fa \
-db /scratch/monoallelic/hm257/leuven_rna/bumblebee_database/final_GFF_bter_1.0_dedup.fasta \
-max_target_seqs 1 \
-outfmt 6 \
-out honeybee_genome_to_bb.txt

blastn \
-query a.mellifera_genes_fasta.fa \
-db /scratch/monoallelic/hm257/leuven_rna/bumblebee_database/final_GFF_bter_1.0_dedup.fasta \
-max_target_seqs 1 \
-outfmt 5 \
-out honeybee_genome_to_bb.xml


# Column headers for the .txt output:

# 1.	 qseqid	 query (e.g., gene) sequence id
# 2.	 sseqid	 subject (e.g., reference genome) sequence id
# 3.	 pident	 percentage of identical matches
# 4.	 length	 alignment length
# 5.	 mismatch	 number of mismatches
# 6.	 gapopen	 number of gap openings
# 7.	 qstart	 start of alignment in query
# 8.	 qend	 end of alignment in query
# 9.	 sstart	 start of alignment in subject
# 10.	 send	 end of alignment in subject
# 11.	 evalue	 expect value
# 12.	 bitscore	 bit score

















