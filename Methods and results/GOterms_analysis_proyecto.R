#######
# Script : Analisis de terminos GO
# Author: Sofia Salazar, Diego Ramirez y Evelia Coss
# Date: 01/03/2024
# Description: El siguiente script nos permite realiza la Determinacion funcional de los genes diferencialmente expresados
# a partir de los datos provenientes del alineamiento de STAR a R,
# Primero correr el script "load_data_inR.R"
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: metadata.csv, cuentas de STAR (Terminacion ReadsPerGene.out.tab)
#   - Output: Matriz de cuentas (CSV y RData)
#######

# qlogin 
# module load r/4.0.2
# R

# --- Load packages ----------
library(gprofiler2)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(dplyr)

# --- Load data -----
# Cargar archivos
indir <- "/mnt/Guanina/bioinfo24/Equipo4/proyecto/DEG_results/"
outdir <- "/mnt/Guanina/bioinfo24/Equipo4/proyecto/DEG_results/Go_terms/"
figdir <- '/mnt/Guanina/bioinfo24/Equipo4/proyecto/DEG_results/'

# ---- Analisis de terminos Go ----
# Seleccionar bases de datos
sources_db <- c("GO:BP", "KEGG", "REAC", "TF", "MIRNA", "CORUM", "HP", "HPA", "WP")

# Seleccionar solo archivos CSV
files <- dir(indir, pattern = "^DE_(.+)\\.csv$")

# ---- Ejemplo de UN SOLO ARCHIVO --------
# Extraer el nombre del primer archivo
plot_name <- gsub("^DE_(.+)\\.csv$", "\\1",  files[1]) #name

# Cargar archivo
df <- read.csv(file = paste0(indir, files[1]), row.names = 'X')
head(df)

# baseMean log2FoldChange     lfcSE       stat    pvalue
# ENSG00000223972  4.668899    -1.89245500 1.6709911 -1.1325345 0.2574098
# ENSG00000227232 62.673986    -0.06534508 0.5778994 -0.1130735 0.9099723
# ENSG00000278267  2.561158    -1.65959750 1.7891499 -0.9275900 0.3536203
# ENSG00000238009  1.948419     2.08774590 2.6846657  0.7776558 0.4367719
# ENSG00000233750  5.225040    -0.39508275 2.5358410 -0.1557995 0.8761911
# ENSG00000268903  8.645805     0.46902662 1.5152980  0.3095276 0.7569202
# padj
# ENSG00000223972 0.8292930
# ENSG00000227232 0.9936289
# ENSG00000278267        NA
# ENSG00000238009        NA
# ENSG00000233750 0.9891925
# ENSG00000268903 0.9778736

# Agregar informacion sobre la expresion
abslogFC <- 2 # Corte de 2 log2FoldChange
df <- df %>% 
  dplyr::mutate(Expression = case_when(log2FoldChange >= abslogFC & padj < 0.05 ~ "Up-regulated",
                                       log2FoldChange <= -(abslogFC) & padj < 0.05 ~ "Down-regulated",
                                       TRUE ~ "Unchanged")) 

# Obtener los nombres de los genes
# > UP 
up_genes <- df %>% filter(Expression == 'Up-regulated') %>% 
  arrange(padj, desc(abs(log2FoldChange)))
# Extraer solo el nombre de los genes
up_genes <- rownames(up_genes) 

head(up_genes)
# [1] "ENSG00000188404" "ENSG00000181847" "ENSG00000099953" "ENSG00000224950"
# [5] "ENSG00000138755" "ENSG00000166509"
write.csv(up_genes, file=paste0(outdir, 'up_genes.csv'))

# > Down
down_genes <- df %>% filter(Expression == 'Down-regulated') %>% 
  arrange(padj, desc(abs(log2FoldChange)))
# Extraer solo el nombre de los genes
down_genes <- rownames(down_genes) 

head(down_genes)
# [1] "ENSG00000250961" "ENSG00000109107" "ENSG00000135447" "ENSG00000223375"
# [5] "ENSG00000159387" "ENSG00000109846"

write.csv(down_genes, file=paste0(outdir, 'down_genes.csv')) # Genes menos expresados 

# 
multi_gp <- gost(list("Upregulated" = up_genes, 
                      "Downregulated" = down_genes), 
                 correction_method = "fdr", user_threshold = 0.05,
                 multi_query = F, ordered_query = T, 
                 sources = sources_db, 
                 evcodes = TRUE,  # intersection = intersection - a comma separated list of genes 
                 # from the query that are annotated to the corresponding term
                 organism = 'hsapiens') # para humano es hsapiens


## ---- colors ---
# paleta de colores
Category_colors <- data.frame(
  category = c("GO:BP", "GO:CC", "GO:MF", "KEGG",
               'REAC', 'TF', 'MIRNA', 'HPA', 'CORUM', 'HP', 'WP'), 
  label = c('Biological Process', 'Cellular Component', 'Molecular Function',  "KEGG",
            'REAC', 'TF', 'MIRNA', 'HPA', 'CORUM', 'HP', 'WP'),
  colors =  c('#FF9900', '#109618','#DC3912', '#DD4477',
              '#3366CC','#5574A6', '#22AA99', '#6633CC', '#66AA00', '#990099', '#0099C6'))

## ----manhattan plot--------
gostp1 <- gostplot(multi_gp, interactive = FALSE)

# Guardar grafica
ggsave(paste0(figdir, "ManhattanGO_", plot_name, ".png"),
       plot = gostp1, dpi = 300)

## ----Dataframe de todos los datos --------
# Convertir a dataframe
gost_query <- as.data.frame(multi_gp$result)

# Extarer informacion en modo matriz de todos los resultados
bar_data <- data.frame("term" = as.factor(gost_query$term_name), "condition" = gost_query$query, 
                       "count" = gost_query$term_size, "p.adjust" = gost_query$p_value, 
                       'category' = as.factor(gost_query$source), "go_id" = as.factor(gost_query$term_id),
                       'geneNames' = gost_query$intersection
)


## ---- DOWN genes ----
bar_data_down <- subset(bar_data, condition == 'Downregulated')

# Ordenar datos y seleccion por pvalue
bar_data_down <-head(bar_data_down[order(bar_data_down$p.adjust),],40) # order by pvalue
bar_data_down_ordered <- bar_data_down[order(bar_data_down$p.adjust),] # order by pvalue
bar_data_down_ordered<- bar_data_down_ordered[order(bar_data_down_ordered$category),] # order by category
bar_data_down_ordered$p.val <- round(-log10(bar_data_down_ordered$p.adjust), 2)
bar_data_down_ordered$num <- seq(1:nrow(bar_data_down_ordered)) # num category for plot

# Guardar dataset
save(bar_data_down_ordered, file = paste0(outdir, "DOWN_GO_", plot_name, ".RData"))

# agregar colores para la grafica
bar_data_down_ordered_mod <- left_join(bar_data_down_ordered, Category_colors, by= "category")

### ---- DOWN genes (barplot) ----
# Generar la grafica
g.down <- ggplot(bar_data_down_ordered_mod, aes(p.val, reorder(term, -num), fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = p.val),
    color = "black",
    hjust = 0,
    size = 2.2,
    position = position_dodge(0)
  ) +
  labs(x = "-log10(p-value)" , y = NULL) +
  scale_fill_manual(name='Category', 
                    labels = unique(bar_data_down_ordered_mod$label), 
                    values = unique(bar_data_down_ordered_mod$colors)) +
  theme(
    legend.position = "right",
    axis.title.y = element_blank(),
    strip.text.x = element_text(size = 11, face = "bold"),
    strip.background = element_blank()
  )+ theme_classic()


# Guardar la figura
ggsave(paste0(figdir,"barplotDOWN_GO_", plot_name, ".png"),
       plot = g.down + theme_classic(), dpi = 600, width = 10, height = 5)



## ---- UP genes ----
bar_data_up <- subset(bar_data, condition == 'Upregulated')

# Ordenar datos y seleccion por pvalue
bar_data_up <-head(bar_data_up[order(bar_data_up$p.adjust),],40) # order by pvalue
bar_data_up_ordered <- bar_data_up[order(bar_data_up$p.adjust),] # order by pvalue
bar_data_up_ordered<- bar_data_up_ordered[order(bar_data_up_ordered$category),] # order by category
bar_data_up_ordered$p.val <- round(-log10(bar_data_up_ordered$p.adjust), 2)
bar_data_up_ordered$num <- seq(1:nrow(bar_data_up_ordered)) # num category for plot

# Guardar dataset
save(bar_data_up_ordered, file = paste0(outdir, "UP_GO_", plot_name, ".RData"))

# agregar colores para la grafica
bar_data_up_ordered_mod <- left_join(bar_data_up_ordered, Category_colors, by= "category")

### ---- UP genes (barplot) ----
# Generar la grafica
g.up <- ggplot(bar_data_up_ordered_mod, aes(p.val, reorder(term, -num), fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = p.val),
    color = "black",
    hjust = 0,
    size = 2.2,
    position = position_dodge(0)
  ) +
  labs(x = "-log10(p-value)" , y = NULL) +
  scale_fill_manual(name='Category', 
                    labels = unique(bar_data_up_ordered_mod$label), 
                    values = unique(bar_data_up_ordered_mod$colors)) +
  theme(
    legend.position = "right",
    # panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x = element_text(size = 11, face = "bold"),
    strip.background = element_blank() 
  ) + theme_classic()

# Guardar la figura
ggsave(paste0(figdir, "barplotUP_GO_", plot_name, ".png"),
       plot = g.up + theme_classic(), dpi = 600, width = 10, height = 5)

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
#   [1] forcats_0.5.1               stringr_1.4.0              
# [3] purrr_0.3.4                 readr_1.4.0                
# [5] tidyr_1.2.1                 tibble_3.1.3               
# [7] tidyverse_1.3.0             clusterProfiler_3.18.1     
# [9] DOSE_3.16.0                 enrichplot_1.10.2          
# [11] gprofiler2_0.2.0            ggplot2_3.3.5              
# [13] pheatmap_1.0.12             dplyr_1.0.10               
# [15] DESeq2_1.30.1               SummarizedExperiment_1.20.0
# [17] Biobase_2.50.0              MatrixGenerics_1.2.1       
# [19] matrixStats_1.0.0           GenomicRanges_1.42.0       
# [21] GenomeInfoDb_1.26.7         IRanges_2.24.1             
# [23] S4Vectors_0.28.1            BiocGenerics_0.36.1        
# 
# loaded via a namespace (and not attached):
#   [1] fgsea_1.16.0           colorspace_2.0-2       ellipsis_0.3.2        
# [4] qvalue_2.22.0          XVector_0.30.0         fs_1.5.0              
# [7] rstudioapi_0.13        farver_2.1.0           graphlayouts_0.7.1    
# [10] ggrepel_0.9.1          bit64_4.0.5            AnnotationDbi_1.52.0  
# [13] fansi_0.5.0            scatterpie_0.1.6       lubridate_1.7.9.2     
# [16] xml2_1.3.2             splines_4.0.2          cachem_1.0.5          
# [19] GOSemSim_2.16.1        geneplotter_1.68.0     polyclip_1.10-0       
# [22] jsonlite_1.7.2         broom_0.7.9            annotate_1.68.0       
# [25] GO.db_3.12.1           dbplyr_2.2.1           ggforce_0.3.3         
# [28] BiocManager_1.30.21    compiler_4.0.2         httr_1.4.2            
# [31] rvcheck_0.1.8          backports_1.2.1        assertthat_0.2.1      
# [34] Matrix_1.3-4           fastmap_1.1.0          lazyeval_0.2.2        
# [37] cli_3.6.0              tweenr_1.0.2           htmltools_0.5.1.1     
# [40] tools_4.0.2            igraph_1.2.8           gtable_0.3.0          
# [43] glue_1.4.2             GenomeInfoDbData_1.2.4 reshape2_1.4.4        
# [46] DO.db_2.9              fastmatch_1.1-0        Rcpp_1.0.7            
# [49] cellranger_1.1.0       vctrs_0.5.1            ggraph_2.0.5          
# [52] rvest_0.3.6            lifecycle_1.0.3        XML_3.99-0.6          
# [55] zlibbioc_1.36.0        MASS_7.3-53            scales_1.1.1          
# [58] tidygraph_1.2.0        hms_1.1.0              RColorBrewer_1.1-2    
# [61] memoise_2.0.0          gridExtra_2.3          downloader_0.4        
# [64] stringi_1.6.2          RSQLite_2.2.7          genefilter_1.72.1     
# [67] BiocParallel_1.24.1    rlang_1.0.6            pkgconfig_2.0.3       
# [70] bitops_1.0-7           lattice_0.20-41        htmlwidgets_1.5.3     
# [73] labeling_0.4.2         cowplot_1.1.1          shadowtext_0.0.8      
# [76] bit_4.0.4              tidyselect_1.2.0       plyr_1.8.6            
# [79] magrittr_2.0.1         R6_2.5.0               generics_0.1.0        
# [82] DelayedArray_0.16.3    DBI_1.1.1              pillar_1.6.2          
# [85] haven_2.3.1            withr_2.4.2            survival_3.5-5        
# [88] RCurl_1.98-1.3         modelr_0.1.8           crayon_1.4.1          
# [91] utf8_1.2.2             plotly_4.9.3           viridis_0.5.1         
# [94] readxl_1.3.1           locfit_1.5-9.4         grid_4.0.2            
# [97] data.table_1.14.0      blob_1.2.1             reprex_1.0.0          
# [100] digest_0.6.27          xtable_1.8-4           munsell_0.5.0         
# [103] viridisLite_0.4.0     

  


