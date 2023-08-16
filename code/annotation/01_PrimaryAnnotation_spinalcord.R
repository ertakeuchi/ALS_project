#' ---
#' title: "01 - Primary annotation"
#' author: "ertakeuchi"
#' date: "15 August 2023"
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
#'     value: spinalcord
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
library("Signac")
set.seed(1234)

#' 
#' ### 2. Set directories
## -----------------------------------------------------------------------------
homeDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq")
NGSDir <- file.path(homeDir, "../../NGS_original")
metadataDir <- file.path(homeDir, "metadata")

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if the correct number of arguments is provided
if (length(args) != 5) {
  stop("Usage: Rscript test.R tissue output_directory", call. = FALSE)
}
# Extract the tissue argument
tissue <- args[1]
directory <- args[2]
feature <- args[3]
dim <- args[4]
resolution <- args[5]

# # set parameters
# tissue <- "spinalcord"
# directory <- "230815_experiment1/out_1"
# feature <- 3000
# dim <- 50
# resolution <- 0.6

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
#' ### 01. Load Seurat object
## -----------------------------------------------------------------------------
# Read the harmonized Seurat object
so <- readRDS(file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, "_withoutRN01_SCT.rds")))

#' 
#' ### 02. Add Primary annotation
## ---- fig.height=4, fig.width=8-----------------------------------------------

cluster_annotations <- list(
  '1' = 'Oligo',
  '2' = 'Oligo',
  '3' = 'Astro',
  '4' = 'Micro',
  '5' = 'Astro',
  '6' = 'Opc',
  '7' = 'Oligo',
  '8' = 'Oligo',
  '9' = 'Oligo',
  '10' = 'Neuron',
  '11' = 'Neuron',
  '12' = 'Oligo',
  '13' = 'Micro',
  '14' = 'Neuron',
  '15' = 'Endo',
  '16' = 'Neuron',
  '17' = 'Endo',
  '18' = 'Micro',
  '19' = 'Lymphocyte',
  '20' = 'MNs',
  '21' = 'Ependymal',
  '22' = 'Oligo'
)

# Add the annotation to seurat metadata
so$PrimaryAnnotation <- unlist(cluster_annotations[so$seurat_clusters])

celltypes <- unique(so$PrimaryAnnotation)
m <- so@meta.data

for (celltype in celltypes) {
  if(celltype == "Oligo" | celltype == "Opc"){
    m[which(m$PrimaryAnnotation == celltype), "subPrimaryAnnotation"] <- "Oligos"
  } else if(celltype == "Neuron" | celltype == "MNs"){
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
#' 
#' ### 03. Dimplot with the primary annotation
## ---- fig.height=10, fig.width=8----------------------------------------------
clusters <- c("PrimaryAnnotation", "subPrimaryAnnotation")

for (i in seq_along(clusters)){
  # Dimplot
  # i = 1
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
  file = file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, "_withoutRN01_PrimaryAnnotation.rds"))
  )

SaveH5Seurat(
  so,
  filename = file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, "_withoutRN01_PrimaryAnnotation.h5seurat")),
  overwrite = TRUE
  )

Convert(
  file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, "_withoutRN01_PrimaryAnnotation.h5seurat")),
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
knitr::purl(file.path(outDir, "../01_PrimaryAnnotation.Rmd"), documentation = 2)

