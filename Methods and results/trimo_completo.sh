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

cd /mnt/Guanina/bioinfo24/Equipo4/proyect/Raw_data 

trimmomatic PE -threads 8 -phred33 SRR10333578_1.fastq SRR10333578_2.fastq \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333578_1_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333578_1_unpaired.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333578_2_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333578_2_unpaired.fq.gz \
ILLUMINACLIP:/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/adaptadores/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:80 

trimmomatic PE -threads 8 -phred33 SRR10333579_1.fastq SRR10333579_2.fastq \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333579_1_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333579_1_unpaired.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333579_2_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333579_2_unpaired.fq.gz \
ILLUMINACLIP:/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed//adaptadores/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:80

trimmomatic PE -threads 8 -phred33 SRR10333580_1.fastq SRR10333580_2.fastq \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333580_1_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333580_1_unpaired.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333580_2_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333580_2_unpaired.fq.gz \
ILLUMINACLIP:/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/adaptadores/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:80

trimmomatic PE -threads 8 -phred33 SRR10333581_1.fastq SRR10333581_2.fastq \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333581_1_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333581_1_unpaired.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333581_2_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333581_2_unpaired.fq.gz \
ILLUMINACLIP:/mnt/Guanina/bioinfo24/Equipo4/proyeto/Data_trimmed/adaptadores/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:80

trimmomatic PE -threads 8 -phred33 SRR10333582_1.fastq SRR10333582_2.fastq \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333582_1_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333582_1_unpaired.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333582_2_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333582_2_unpaired.fq.gz \
ILLUMINACLIP:/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/adaptadores/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:80

trimmomatic PE -threads 8 -phred33 SRR10333583_1.fastq SRR10333583_2.fastq \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333583_1_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333583_1_unpaired.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333583_2_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed//SRR10333583_2_unpaired.fq.gz \
ILLUMINACLIP:/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/adaptadores/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:80
