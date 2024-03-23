#!/bin/bash
#
# Your job name
#$ -N index
#
# Use current working directory
#$ -cwd
#
# Join stdout and stderr
#$ -j y
#
# Run job through bash shell 
#$ -S /bin/bash
#
# Send an email after the job has finished
#$ -m e
#$ -M valery2602v@gmail.com
#
# If modules are needed, source modules environment (Do not delete the next line):
. /etc/profile.d/modules.sh
#
# Add any modules you might require:
module load star/2.7.9a 
#
# The amount of memory requested per core (Default for many applications is 2GB).
# -l virtual_free=4G
#
# Write your commands in the next line 

STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir /mnt/Guanina/bioinfo24/Equipo4/proyecto/SART_index \
--genomeFastaFiles /mnt/Archives/genome/human/GRCh38/ensembl76/chromosomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /mnt/Archives/genome/human/GRCh38/ensembl76/GTF-file/Homo_sapiens.GRCh38.76.gtf \
--sjdbOverhang 149

