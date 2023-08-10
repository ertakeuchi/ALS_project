#' ---
#' title: "01 - Primary annotation"
#' author: "ertakeuchi"
#' date: "10 August 2023"
#' output: 
#'   html_document:
#'     # code_folding: "hide"
#'     # fig_caption: "true"
#' theme: 
#'   sketchy: "true"
#'   highlight: "tango"
#' warnings: "false"
#' params: 
#'   tissue_name: 
#'     value: brain
#'     choices:
#'       - spinalcord
#'       - brain
#' ----
#' 
## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' ## 01. Set environment
#' 
#' ### 1. Load libraries
## ---- include=FALSE-----------------------------------------------------------
library("Seurat")
library("SeuratDisk")
library("tidyverse")
library("reticulate")
use_python("/usr/bin/python3")
library("patchwork")
library("gridExtra")
library("leiden")
library("Signac")
set.seed(1234)

#' 
#' ### 2. Set directories
## -----------------------------------------------------------------------------
homeDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq")
NGSDir <- file.path(homeDir, "../../NGS_original")
metadataDir <- file.path(homeDir, "metadata")

# # Get command-line arguments
# args <- commandArgs(trailingOnly = TRUE)
# # Check if the correct number of arguments is provided
# if (length(args) != 5) {
#   stop("Usage: Rscript test.R tissue output_directory", call. = FALSE)
# }
# # Extract the tissue argument
# tissue <- args[1]
# directory <- args[2]
# feature <- args[3]
# dim <- args[4]
# resolution <- args[5]

tissue <- "brain"
directory <- "230810_experiment2/out_2"
outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)


#' 
#' 
#' ################################################################################
#' # Annotation 01: Primary annotation 
#' ################################################################################
#' 
#' ## 02. Primary Annotation
## -----------------------------------------------------------------------------
# 01: Set environment
# |
# |-- 02: Primary annotation 
# |   |
# |   |-- 01: Load Seurat object 
# |   |
# |   |-- 02: Add Primary annotation
# |   |
# |   |-- 03: Dimplot with the primary annotation
# |   |
# |   |-- 04: Save Rdata, h5seurat, h5ad

#' 
#' 
#' ### 01. Load Seurat object
## -----------------------------------------------------------------------------
# set decieded parameters
feature <- 5000
dim <- 50
resolution <- 0.6

# Read the harmonized Seurat object
so <- readRDS(file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, "_SCT.rds")))
# # levels(so)
# sorted_levels <- sort(as.numeric(levels(so))) %>% as.character()
# so$seurat_clusters <- factor(x = so$seurat_clusters, levels = sorted_levels)

f <- factor(
  so@meta.data$seurat_clusters,
  levels = seq_along(unique(so@meta.data$seurat_clusters))
  )
head(f)
identical(as.vector(f), as.vector(so@meta.data$SCT_snn_res.0.6))

so@meta.data$seurat_clusters <- f

p <- DimPlot(so, reduction = "umap_HM_01", label = TRUE, repel = TRUE, group.by = "seurat_clusters", raster = TRUE)
p


marker_genes <- list(
  'Neuron' = c("RBFOX3"),
  'GABAergic' = c("GAD1", "GAD2"), 
  'Glutamatergic' = c("SLC17A7", "SATB2"), 
  'OPC' = c("PDGFRA", "VCAN"), 
  'Astro' = c("AQP4", "GFAP", "ATP1A2", "SLC1A2"),
  'Oligo' = c("PLP1", "MOBP", "OPALIN", "MBP"), 
  'Perivasc_Macro' = c("MRC1"), 
  'Tcells' = c("PTPRC"), 
  'Vasc_smooth' = c("PDGFRB"), 
  'Vasc_endo' = c("FLT1"), 
  'Vasc_fibro' = c("DCN"), 
  'Micro' = c("APBB1IP", "CX3CR1", "P2RY12", "TREM2", "CSF1R", "CD74", "C3")
)
marker_gene_lst <- c(unlist(marker_genes, use.names = F))

# Dotplot
p1 <- DotPlot(
  so,
  features = marker_gene_lst,
  group.by = "seurat_clusters",
  dot.scale = 10
  ) &
  labs("Canonical markers") &
  theme(
    axis.text.y=element_text(hjust=0, size = 15),
    axis.title.y=element_blank(), 
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 15)
  )
p1

ggsave(
  file.path(outDir, paste0("MarkerGenesDotPlot.pdf")),
  plot = p1,
  height = 10, width = 10)


#' 
#' ### 02. Add Primary annotation
## ---- fig.height=4, fig.width=8-----------------------------------------------

cluster_annotations <- list(
  '1' = 'Oligo',
  '2' = 'Oligo',
  '3' = 'EX',
  '4' = 'Oligo',
  '5' = 'Astro',
  '6' = 'EX',
  '7' = 'Micro',
  '8' = 'INH',
  '9' = 'Oligo',
  '10' = 'Opc',
  '11' = 'INH',
  '12' = 'INH',
  '13' = 'Oligo',
  '14' = 'EX',
  '15' = 'EX',
  '16' = 'EX',
  '17' = 'INH',
  '18' = 'EX',
  '19' = 'Endo',
  '20' = 'EX',
  '21' = 'INH',
  '22' = 'EX',
  '23' = 'EX'

)

# Add the annotation to seurat metadata
so$PrimaryAnnotation <- unlist(cluster_annotations[so$seurat_clusters])

celltypes <- unique(so$PrimaryAnnotation)
m <- so@meta.data

for (celltype in celltypes) {
  if(celltype == "Oligo" | celltype == "Opc"){
    m[which(m$PrimaryAnnotation == celltype), "subPrimaryAnnotation"] <- "Oligos"
  } else if(celltype == "EX" | celltype == "INH"){
    m[which(m$PrimaryAnnotation == celltype), "subPrimaryAnnotation"] <- "Neurons"
  } else if(celltype == "Astro" ){
    m[which(m$PrimaryAnnotation == celltype), "subPrimaryAnnotation"] <- "Astrocytes"
  } else if(celltype == "Micro" ){
    m[which(m$PrimaryAnnotation == celltype), "subPrimaryAnnotation"] <- "Micros"
  } else {
    m[which(m$PrimaryAnnotation == celltype), "subPrimaryAnnotation"] <- "Others"
  } 
}
  
head(m)
# unique(m$subPrimaryAnnotation)

so@meta.data <- m
head(so@meta.data)

#' 
#' ### 03. Dimplot with the primary annotation
## ---- fig.height=10, fig.width=8----------------------------------------------
clusters <- c("PrimaryAnnotation", "subPrimaryAnnotation")

for (i in seq_along(clusters)){
  # Dimplot
  group <- clusters[i]
  print(group)

  p1 <- DimPlot(
    so,
    reduction = "umap_HM_01",
    repel = TRUE,
    label = TRUE,
    group.by = group,
    raster = TRUE
  ) & 
    ggtitle('UMAP colored by primary annotations') 
  # p1

  ggsave(
    file.path(outDir, paste0(group, "_DimPlot.pdf")),
    plot = p1,
    height = 10, width = 10
  )
}


#' 
#' 
#' ### 04. Save Rdata, h5seurat, h5ad
#' 
## -----------------------------------------------------------------------------
saveRDS(
  so,
  file = file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, "_PrimaryAnnotation.rds"))
  )

SaveH5Seurat(
  so,
  filename = file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, "_PrimaryAnnotation.h5seurat")),
  overwrite = TRUE
  )

Convert(
  file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, "_PrimaryAnnotation.h5seurat")),
  dest = "h5ad",
  overwrite = TRUE
  )

#' 
#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../01_AnnotationPrimary.Rmd"), documentation = 2)

