# RNA-seq 
## Made by Valeria M. Cuevas and Renata S. CHiw
Bioinformatic proyect with RNA-seq data from the papper "Identification of differentially expressed lncRNAs and mRNAs in luminal-B breast cancer by RNA-sequencing".
|Bioproject |Species            |Type of libraries|
|-----------|-------------------|-----------------|
|PRJNA579053|Homo sapiens(human)| paired-end      |

|Selection method     |Number of transcriptomes| Number of biological replicates|
|---------------------|------------------------|--------------------------------|
|PAXgene blood RNA kit|~30 milions             |4

|Sequencing instrument        |Distribution of samples  |Sequencing depth of each transcriptome|
|-----------------------------|-------------------------|--------------------------------------|
|Illumina Hiseq X-ten platform| 4 samples and 4 controls| 10G depth 

## Data download

> ### How to download the data?

The data for this article are available in the Gene Expression Omnibus (GEO) (GSE139274, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139274). The data were downloaded from the European Nucleotide Archive (ENA) in fastq.gz format with **the next command**:

```
cd /mnt/Guanina/bioinfo24/Equipo4/proyecto/Raw_data/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/079/SRR10333579/SRR10333579_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/079/SRR10333579/SRR10333579_2.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/082/SRR10333582/SRR10333582_1.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/082/SRR10333582/SRR10333582_2.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/084/SRR10333584/SRR10333584_1.fastq.gz 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/084/SRR10333584/SRR10333584_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/078/SRR10333578/SRR10333578_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/078/SRR10333578/SRR10333578_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/080/SRR10333580/SRR10333580_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/080/SRR10333580/SRR10333580_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/081/SRR10333581/SRR10333581_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/081/SRR10333581/SRR10333581_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/083/SRR10333583/SRR10333583_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/083/SRR10333583/SRR10333583_2.fastq.gz

```
Once the files are downloaded we unzip the files with the following command: 

```
gzip -d *.fastq.gz
```

The commands are in the [human.sh](https://github.com/rensanch/RNA-seq-team-/blob/main/Methods%20and%20results/humano.sh) script that is in the Methods and results directory.

The rute that directs you to the .fastq files in the Raw_data folder is : 

```
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Raw_data
```

## Data quality

In the fastqc reports from the raw data, we observed that all samples presented good quality in general. There are errors in the modules: Per base sequence content, Sequence Duplication Levels, Per sequence GC content and Overrepresented sequences since we are working with RNA-seq and these modules tend to fail with this type of experiment. However, errors in the Per sequence GC content and Overrepresented sequence modules may be correlated. 

When we visualized the fastqc and searched for some of these overrepresented sequences in blast https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome we found that they are not adapters but are already highly conserved sequences that occur in other species besides humans. These conserved sequences that are overrepresented are indicators of contamination which explains the peaks in the per sequence GC content. Some of these sequences are lnRNA, which are the elements of main interest in the article.

We decided to preserve these sequences because in a significant number of files they represent 5% to 10%. 

In the multiqc, in the same way as in the fastqc, good quality is observed and what was previously reported but in a more general way. A point to highlight is that we observed a considerable percentage of duplicates in all files.

> ### How to create fastqc files and multiqc file?

We use the next command to generate fastqc files:

```
module load fastqc/0.12.1
cd /mnt/Guanina/bioinfo24/Equipo4/proyecto/Raw_data
fastqc -o /mnt/Guanina/bioinfo24/Equipo4/proyecto/FastQC_rawData *.fastq
```
The commands for the fastqc are in the [fastqc.sh](https://github.com/rensanch/RNA-seq-team-/blob/main/Methods%20and%20results/fasqc.sh) script that is in the Methods and results directory.

We use the next command to generate multiqc file:

```
module load multiqc/1.5
cd /mnt/Guanina/bioinfo24/Equipo4/proyecto/FastQC_rawData
multiqc .
```

The rute that directs you to the .fastqc files and multiqc in the FastQC_rawData folder is : 

```
/mnt/Guanina/bioinfo24/Equipo4/proyecto/FastQC_rawData
```

## Trimming



