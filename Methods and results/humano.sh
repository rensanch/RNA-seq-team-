#!/bin/bash
#
# Your job name
#$ -N humano
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

#
# The amount of memory requested per core (Default for many applications is 2GB).
# -l virtual_free=4G
#
# Write your commands in the next line 
cd /mnt/Guanina/bioinfo24/Equipo4/proyecto/Raw_data/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/079/SRR10333579/SRR10333579_1.fastq.gz ###Usando wget estamos indicando que se 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/079/SRR10333579/SRR10333579_2.fastq.gz ###descargaran archivos de internet
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/082/SRR10333582/SRR10333582_1.fastq.gz ###y despues de wget se coloca la 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/082/SRR10333582/SRR10333582_2.fastq.gz ###direccion de donde se encuentra dicho
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/084/SRR10333584/SRR10333584_1.fastq.gz ###archivo
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/084/SRR10333584/SRR10333584_2.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/078/SRR10333578/SRR10333578_1.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/078/SRR10333578/SRR10333578_2.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/080/SRR10333580/SRR10333580_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/080/SRR10333580/SRR10333580_2.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/081/SRR10333581/SRR10333581_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/081/SRR10333581/SRR10333581_2.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/083/SRR10333583/SRR10333583_1.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/083/SRR10333583/SRR10333583_2.fastq.gz
