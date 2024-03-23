######
# Script : Analisis de expresion diferencial
# Author: Sofia Salazar, Diego Ramirez y Evelia Coss
# Date: 27/02/2024
# Description: El siguiente script nos permite realiza el Analisis de expresion Diferencial
# a partir de los datos provenientes del alineamiento de STAR a R,
# Primero correr el script "load_data_inR.R"
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: Cargar la variable raw_counts.RData que contiene la matriz de cuentas y la metadata
#   - Output: DEG
#######

# qlogin 
# module load r/4.0.2
# R

# --- Load packages ----------
library(DESeq2)

# --- Load data -----
# Cargar archivos
outdir <- "/mnt/Guanina/bioinfo24/Equipo4/proyecto/DEG_results/"
figdir <- '/mnt/Guanina/bioinfo24/Equipo4/proyecto/DEG_results/'

#Cargar variable "counts", proveniente del script "load_data_inR.R"
load("/mnt/Guanina/bioinfo24/Equipo4/proyecto/DEG_results/counts/raw_counts.RData") 
samples <- metadata$sample_id # Extraer los nombres de los Transcriptomas
metadata$type <- as.factor(metadata$type) # convertir a factor

# --- DEG ---- 
#Prefiltrado
counts <- counts[which(rowSums(counts) > 10),] #Seleccionamos genes con mas de 10 cuentas

# Convertir al formato dds
dds <- DESeqDataSetFromMatrix(countData =  counts,
            colData = metadata, design = ~type) #Se hace un DESeqDataSet para realizar un analisis

dim(dds) # checar las dimensiones
#[1] 33207     8

##  -- Asignar la referencia y generar contrastes -----
# Las comparaciones se realizan por pares
#Si no se indica de manera explicita que se va a comparara, lo va a tomar de manera alfabetica, 
# en este caso se indica que control es la referencia, 
dds$type <- relevel(dds$type, ref = "control") 

## --- Obtener archivo dds ----

dds <- DESeq(dds)

# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

# Obtener la lista de coeficientes o contrastes
resultsNames(dds)

# [1] "Intercept"                     "type_breast_cancer_vs_control"

# Guardar la salida del diseno
save(metadata, dds, file = paste0(outdir, 'dds_BreastCancer_vs_control.RData'))

## --- Normalizacion de los datos ---------
# Opcion 1. log2(n + 1)
# ntd <- normTransform(dds)

# Opcion 2. regularized logarithm or rlog
# Normalizacion de las cuentas por logaritmo y podrias hacer el analisis usando este objeto en lugar del dds
ddslog <- rlog(dds, blind = F) 

# Opcion 3. vsd
# Estima la tendencia de dispersion de los datos y calcula la varianza, hace una normalizacion de las 
# cuentas con respecto al tamaño de la libreria
vsdata <- vst(dds, blind = F) 

## --- Deteccion de batch effect ----

# Almacenar la grafica
png(file = paste0(figdir, "PCA_rlog.png"))
plt <- plotPCA(ddslog, intgroup = "type")
print(plt)
dev.off()

# plotPCA
# Arguments:
# Input: a ‘DESeqTransform’ object, with data in ‘assay(x)’, produced
# for example by either ‘rlog’ or ‘varianceStabilizingTransformation’.


# Almacenar la grafica
png(file = paste0(figdir, "PCA_vsd.png"))
plt <- plotPCA(vsdata, intgroup = "type")
print(plt)
dev.off()

# Guardar la salida del diseno (vsdata)
save(metadata, vsdata, file = paste0(outdir, 'vst_BreastCancer_vs_control.RData'))


## ---- Obtener informacion del contraste 1 ----
# results(dds, contrast=c("condition","treated","untreated"))
breastCancer <- results(dds, name = "type_breast_cancer_vs_control")
breastCancer

summary(breastCancer)

# out of 33470 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1657, 5%
# LFC < 0 (down)     : 811, 2.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 18169, 54%
# (mean count < 12)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Guardar los resultados
write.csv(breastCancer, file=paste0(outdir, 'DE_BreastCancer_vs_control.csv'))

# Guardar informacion de ejecucion
sessionInfo() 
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRblas.so
# LAPACK: /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     
# 
# other attached packages:
#   [1] DESeq2_1.30.1               SummarizedExperiment_1.20.0
# [3] Biobase_2.50.0              MatrixGenerics_1.2.1       
# [5] matrixStats_1.0.0           GenomicRanges_1.42.0       
# [7] GenomeInfoDb_1.26.7         IRanges_2.24.1             
# [9] S4Vectors_0.28.1            BiocGenerics_0.36.1        
# 
# loaded via a namespace (and not attached):
#   [1] genefilter_1.72.1      locfit_1.5-9.4         tidyselect_1.2.0      
# [4] splines_4.0.2          lattice_0.20-41        generics_0.1.0        
# [7] colorspace_2.0-2       vctrs_0.5.1            utf8_1.2.2            
# [10] blob_1.2.1             XML_3.99-0.6           survival_3.5-5        
# [13] rlang_1.0.6            pillar_1.6.2           glue_1.4.2            
# [16] DBI_1.1.1              BiocParallel_1.24.1    bit64_4.0.5           
# [19] RColorBrewer_1.1-2     GenomeInfoDbData_1.2.4 lifecycle_1.0.3       
# [22] zlibbioc_1.36.0        munsell_0.5.0          gtable_0.3.0          
# [25] memoise_2.0.0          labeling_0.4.2         geneplotter_1.68.0    
# [28] fastmap_1.1.0          AnnotationDbi_1.52.0   fansi_0.5.0           
# [31] Rcpp_1.0.7             xtable_1.8-4           scales_1.1.1          
# [34] cachem_1.0.5           DelayedArray_0.16.3    annotate_1.68.0       
# [37] XVector_0.30.0         farver_2.1.0           bit_4.0.4             
# [40] digest_0.6.27          ggplot2_3.3.5          dplyr_1.0.10          
# [43] grid_4.0.2             cli_3.6.0              tools_4.0.2           
# [46] bitops_1.0-7           magrittr_2.0.1         RCurl_1.98-1.3        
# [49] RSQLite_2.2.7          tibble_3.1.3           pkgconfig_2.0.3       
# [52] crayon_1.4.1           ellipsis_0.3.2         Matrix_1.3-4          
# [55] assertthat_0.2.1       httr_1.4.2             R6_2.5.0              
# [58] compiler_4.0.2        
