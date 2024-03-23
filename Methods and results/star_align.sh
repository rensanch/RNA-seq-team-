#!/bin/bash
#
# Your job name
#$ -N alignment
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

index=/mnt/Guanina/bioinfo24/Equipo4/proyecto/STAR_index
FILES=/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/*_1_trimmed.fq.gz
for f in $FILES
do
    base=$(basename $f _1_trimmed.fq.gz)
    echo $base
    STAR --runThreadN 12 --genomeDir $index --readFilesIn $f /mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/$base"_2_trimmed.fq.gz" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --readFilesCommand zcat \
    --outFileNamePrefix /mnt/Guanina/bioinfo24/Equipo4/proyecto/STAR_output/$base
done
