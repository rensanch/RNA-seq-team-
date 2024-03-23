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

> ### How to download the data

The data for this article are available in the Gene Expression Omnibus (GEO) (GSE139274, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139274). The data were downloaded from the European Nucleotide Archive (ENA) in fastq.gz format with the next command:

`
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

`



