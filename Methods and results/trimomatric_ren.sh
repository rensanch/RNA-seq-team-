#!/bin/bash
#
# Your job name
#$ -N trimo
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
module load trimmomatic/0.33
#
# The amount of memory requested per core (Default for many applications is 2GB).
# -l virtual_free=4G
#
# Write your commands in the next line 

cd /mnt/Guanina/bioinfo24/Equipo4/proyecto/Raw_data 

trimmomatic PE -threads 8 -phred33 SRR10333584_1.fastq SRR10333584_2.fastq \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333584_1_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333584_1_unpaired.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333584_2_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333584_2_unpaired.fq.gz \
ILLUMINACLIP:/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/adaptadores/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:80 

