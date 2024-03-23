# RNA-seq 
> ### Made by Valeria M. Cuevas and Renata S. CHiw
Bioinformatic project with RNA-seq data from the papper "Identification of differentially expressed lncRNAs and mRNAs in luminal-B breast cancer by RNA-sequencing".

## Data description

|Bioproject |Species            |Type of libraries|
|-----------|-------------------|-----------------|
|PRJNA579053|Homo sapiens(human)| paired-end      |

|Selection method     |Number of transcriptomes| Number of biological replicates|
|---------------------|------------------------|--------------------------------|
|PAXgene blood RNA kit|~30 milions             |4

|Sequencing instrument        |Distribution of samples  |Sequencing depth of each transcriptome|
|-----------------------------|-------------------------|--------------------------------------|
|Illumina Hiseq X-ten platform| 4 samples and 4 controls| 10G depth 

## Abstract

We are comparing genes in people without illuminal-b breast cancer and genes in people with illuminal-b cancer, with the purpose of observing which genes are overexpressed in people with illuminal-b breast cancer. This comparison is through pair-end RNA-seq data from which lnRNA and mRNA are specifically analyzed. For data trimming, alignment, enrichment, differential expression and data visualization, the trimmomatic, STAR and Rstudio programs were used.

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

> ### How to do the trimming process?

For the trimming process the trimmomatic program is used, which includes several steps to remove adapters and filtering the reads by quality. The next command is an example of what we do:

```
cd /mnt/Guanina/bioinfo24/Equipo4/proyecto/Raw_data

trimmomatic PE -threads 8 -phred33 SRR10333578_1.fastq SRR10333578_2.fastq \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333578_1_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333578_1_unpaired.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333578_2_trimmed.fq.gz \
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed/SRR10333578_2_unpaired.fq.gz \
ILLUMINACLIP:/mnt/Guanina/bioinfo24/Equipo4/proyecto/TData_trimmed/adaptadores/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:2$
```
The commands for the trimming are in the [Trimming files](https://github.com/rensanch/RNA-seq-team-/blob/main/Methods%20and%20results/trimo_completo.sh) folder that is in the Methods and results directory.

The rute that directs you to the trimmed files  in the FastQC_rawData folder is : 

```
/mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed
```


## Data quality after trimming

We observed a slight improvement in the quality of the data, since, as mentioned above, the samples initially had a very good quality. Trimmomatic also removed some of the overrepresented sequences and adapters. In general we still have a high degree of duplicate sequences, irregular behavior in the GC content.

> ### How to create fastqc files and multiqc file?

We use the next command to generate fastqc files:

```
module load fastqc/0.12.1
cd /mnt/Guanina/bioinfo24/Equipo4/proyecto/Data_trimmed
fastqc -o /mnt/Guanina/bioinfo24/Equipo4/proyecto/FastQC_trimmed *.fq.gz
```
The commands for the fastqc are in the [fastqc.sh](https://github.com/rensanch/RNA-seq-team-/blob/main/Methods%20and%20results/fastqc2.sh) script that is in the Methods and results directory.

We use the next command to generate multiqc file:

```
module load multiqc/1.5
cd /mnt/Guanina/bioinfo24/Equipo4/proyecto/FastQC_trimmed
multiqc .
```

The rute that directs you to the .fastqc files and multiqc in the FastQC_trimmed folder is : 

```
/mnt/Guanina/bioinfo24/Equipo4/proyecto/FastQC_trimmed
```


## Alignment

> How to align with STAR?

First we make the index of the genome, to create the index, the reference genome and the annotation file are needed. These files are available in the cluster in the following paths:

- Reference genome:

```/mnt/Archives/genome/human/GRCh38/ensembl76/chromosomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa``` 

-Annotations file: 

```/mnt/Archives/genome/human/GRCh38/ensembl76/GTF-file/Homo_sapiens.GRCh38.76.gtf```

-The index is available in: 

```/mnt/Guanina/bioinfo24/Equipo4/proyecto/STAR_index```

The next commands were used to do the index:

```
module load star/2.7.9a
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
```

The rute that directs you to the alignment files in the  folder is : 

```
/mnt/Guanina/bioinfo24/Equipo4/proyecto/STAR_index
```

The commands for the index are in the [Star_index](https://github.com/rensanch/RNA-seq-team-/blob/main/Methods%20and%20results/star_index2.sh) script that is in the Methods and results directory.

The STAR program was used for sequence alignment with the next commands: 

```
STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir /mnt/Guanina/bioinfo24/Equipo4/proyecto/STAR_index \
--genomeFastaFiles /mnt/Archives/genome/human/GRCh38/ensembl76/chromosomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /mnt/Archives/genome/human/GRCh38/ensembl76/GTF-file/Homo_sapiens.GRCh38.76.gtf \
--sjdbOverhang 149
```

The commands for the alignment are in the [Star_align](https://github.com/rensanch/RNA-seq-team-/blob/main/Methods%20and%20results/star_align.sh) script that is in the Methods and results directory.

The rute that directs you to the alignment files in the  folder is : 

```
/mnt/Guanina/bioinfo24/Equipo4/proyecto/STAR_output
```


## DEG 
To perform the differential expression analysis, the control group was taken against the group of samples with luminal B breast cancer. The following scripts were used: 

- [load_data_inR_proyecto.R](https://github.com/rensanch/RNA-seq-team-/blob/main/Methods%20and%20results/load_data_inR_proyecto.R)   
- [DEG_analysis_proyecto.R](https://github.com/rensanch/RNA-seq-team-/blob/main/Methods%20and%20results/DEG_analysis_proyecto.R) 
- [VisualizacionDatos_proyecto.R](https://github.com/rensanch/RNA-seq-team-/blob/main/Methods%20and%20results/VisualizacionDatos_proyecto.R)

## Results 
 
Differential expression analysis of the read count data using DESeq2 identified a total of 148 significant DEGs (padj<0.05 and LFC≥2) between smaples with luminal B breast cancer and controls, of them 41 genes were found to be downregulated, while 107 genes were found to be upregulated ([up_genes.csv](https://github.com/rensanch/RNA-seq-team-/blob/main/Visualization/up_genes.csv) y [down_genes.csv](https://github.com/rensanch/RNA-seq-team-/blob/main/Visualization/down_genes.csv)). 

## PCA 

It was normalized for the second time with the methods: normalization of the accounts by logarithm (rlog) and normalization of the accounts with respect to the size of the library (vsd). What we can observe is that there are clusters between the samples that are controls and the samples with luminal B breast cancer. We do not observe Batch effect.

Principal component analysis ([PCA](https://github.com/rensanch/RNA-seq-team-/blob/main/Visualization/PCA_rlog.png)) of luminal B breast cancer and controls used in this study. The PC1 captured 43% of the total variance and PC2 captured 27% of the total variance. Each dot represents one sample

Principal component analysis ([PCA](https://github.com/rensanch/RNA-seq-team-/blob/main/Visualization/PCA_vsd.png)) of luminal B breast cancer and surrounding normal samples used in this study. The PC1 captured 50% of the total variance and PC2 captured 25% of the total variance. Each dot represents one sample

## Volcano plot 

The [Volcano plot](https://github.com/rensanch/RNA-seq-team-/blob/main/Visualization/VolcanoPlot_BreastCancer_vs_control.png) representing differential expression pattern of Luminal-B breast cancer vs controsls. Red colour represents up regulated genes. Blue colour represents down regulated genes. 

## Heatmaps 

We made a heatmap based on the results of normalization by vsd with the name of the 20 genes with the most significant p value.  

In the [heatmap](https://github.com/rensanch/RNA-seq-team-/blob/main/Visualization/Heatmap_vsd_topgenes.png) we can see the 4 patients (on the x axis the samples with ending 77, 78, 79 and 80) and the 4 controls (on the controls we observed that most of the genes have a low p-value, so we could say that they are not being expressed in the same way as in the patients, of which the genes seem to be overexpressed except for the gene ENSG00000109107 which seems to be overexpressed in controls and not in patients, so when suffering from luminal-b breast cancer the expression of that gene is less significant.

## Functional analysis 

Functional enrichment analysis for 148 DEGs obtained in the above step was carried out using gprofiler2. Using a padj cut off of 0.05, genes were enriched for KEGG pathways and [GO terms](https://github.com/rensanch/RNA-seq-team-/blob/main/Methods%20and%20results/GOterms_analysis_proyecto.R).   

In the [manhattan](https://github.com/rensanch/RNA-seq-team-/blob/main/Visualization/ManhattanGO_BreastCancer_vs_control.png) we can see the down regulated genes and we do not observe anything significant. On the other hand, in the second panel we have the Up regulated genes, we can see that we have genes with significant p-valeus in the category of biological process. 

In this [graph](https://github.com/rensanch/RNA-seq-team-/blob/main/Visualization/barplotUP_GO_BreastCancer_vs_control.png) we observed the p-value of different biological processes, also the biological processes with a higher p-value are related to the functions of certain genes that are overexpressed in patients with luminal-b breast cancer, for example; ENSG00000181847 which is a T-cell immunoreceptor with Ig and ITIM domains that is responsible for suppresses T-cell activation by promoting the generation of mature immunoregulatory dendritic cells. 

In this [graph](https://github.com/rensanch/RNA-seq-team-/blob/main/Visualization/barplotDOWN_GO_BreastCancer_vs_control.png) we observed the p value of the down regulated genes, these are related to various biological processes, but we can see a greater association with processes related to glycolysis and lipid processing. Examples of genes belonging to these categories found by differential expression were PPP1R1A and ALDO, both coding for proteins that participate in glycolysis and gluconeogenesis. 

## Conclusions 

The results obtained in the article and the results achieved with our analysis are very similar. In the article we are provided with the ten main up regulated and down regulated genes, our analysis obtained practically the same genes in these sections. Regarding the functional analysis, the biological processes related to these differentially expressed genes also coincide. In both cases we have a strong expression of genes related to the immune system.  

We compare our results with the following articles: Expression profiling of luminal B breast tumor in Indian women (Ulaganathan et al., 2023), The immune system and inflammation in breast cancer (Jiang & Shapiro, 2014) and Immune System Effects on Breast Cancer ( Amens et al., 2021).  In the first article the authors identified lncRNAs. In our analysis we also identified lnRANs although our data do not coincide, it seems to be possible to see lnRNA characteristics of this subtype of breast cancer. 

Regarding the other two articles, there is evidence that breast cancer generates a greater immune system response compared to other types of cancer, "breast tumor progression is either prevented by the action of antitumor immunity or exacerbated by proinflammatory cytokines released mainly by the immune cells" (Amens et al., 2021).  In conclusion, our data are consistent with the selected article and with the bibliography consulted on the topic.





