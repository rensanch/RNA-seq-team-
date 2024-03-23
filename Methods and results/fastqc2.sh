#!/bin/bash
#
# Your job name
#$ -N fastqc
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
module load fastqc/0.12.1
#
# The amount of memory requested per core (Default for many applications is 2GB).
# -l virtual_free=4G
#
# Write your commands in the next line 
cd /mnt/Guanina/bioinfo24/Equipo4/proyecto/TRIM_results

fastqc -o /mnt/Guanina/bioinfo24/Equipo4/proyecto/quality2 *.fq.gz 
