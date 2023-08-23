#' ---
#' title: "03 - Secondary annotation"
#' author: "ertakeuchi"
#' date: "23 August 2023"
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
#' # 1. Set environment
#' 
#' ## 0-0. Load libraries
## ---- include=FALSE-----------------------------------------------------------
library("Seurat")
library("SeuratDisk")
library("Signac")
library("tidyverse")
library("reticulate")
use_python("/usr/bin/python3")
library("patchwork")
library("progressr")
library("gridExtra")
library("future")
set.seed(1234)


#' 
#' 
#' ## 0-1. Set a common directory (FOR CONTAINER)
## -----------------------------------------------------------------------------
homeDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq")
# homeDir <- file.path(".")
# NGSDir <- file.path(homeDir, "../../NGS_original")
metadataDir <- file.path(homeDir, "metadata")

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("Usage: Rscript test.R tissue output_directory", call. = FALSE)
}
# Extract the tissue argument
# tissue <- args[1]
# directory <- args[2]
# feature <- args[3]
# dim <- args[4]
# resolution <- args[5]

tissue <- "spinalcord"
directory <- "230823_experiment2/out_2"
ref_annotation <- "subtype_annotation"
outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)
resolutions <- c(0.2, 0.4, 0.6, 0.8)

#' 
#' ################################################################################
#' # Annotation 03: Secondary annotation 
#' ################################################################################
#' 
#' # 02. Secondary annotation
#' 
## -----------------------------------------------------------------------------
# 01: Set environment
# |
# |-- 02: Secondary annotation 
# |   |
# |   |-- 01: Load Seurat object using metadata
# |   |
# |   |-- 02: Add Secondary annotation
# |   |  | 
# |   |  |-- 01: Astrocytes
# |   |  |-- 02: Micros
# |   |  |-- 03: Oligodendrocytes
# |   |  |-- 04: Neurons
# |   |  |-- 05: Others


#' 
#' ## 01: Load Seurat object with Primary Annotation
#' 
#' ### 01-1: Decide PCA dimentions and the number of variable genes
## -----------------------------------------------------------------------------
subsets <- read.table(file = file.path(rdataDir, paste0(tissue, "_subsets.txt")), header = FALSE, stringsAsFactors = FALSE)$V1


#' 
#' ## 01-2: Set meta data
## -----------------------------------------------------------------------------
featDimList <- read.table(file = file.path(metadataDir, "featureDimList.tsv"), header = TRUE, stringsAsFactors = FALSE)
featDimList <- featDimList[featDimList$tissue == tissue & featDimList$levels == 2, ]
head(featDimList)

markers_df <- read.table(file = file.path(metadataDir, "canonical_markers.tsv"), header = TRUE, stringsAsFactors = FALSE)
mark1 <- markers_df[markers_df$Level == "Top", ]
mark2 <- markers_df[markers_df$Level == paste0("Sub_", tissue), ]
mark3 <- markers_df[markers_df$Level == paste0("SubPrimary_", tissue), ]

# file.copy(file.path(outDir, "../../230816_experiment2/out_2", "refDataRefMapping.csv"), file.path(rdataDir, "refDataRefMapping.csv"))
metadata <- read.csv(
        file = file.path(rdataDir, "refDataRefMapping.csv"), 
        sep = ",",
        header = TRUE,
        row.names = 1
        )
metadata$Project <- NULL
metadata$barcode <- rownames(metadata)

head(metadata)

#' 
#' 
#' ## 02: Add Secondary annotation
#' 
#' ### 00: Set Secondary Annotation function
## -----------------------------------------------------------------------------
runAnalysis <- function(subPrimaryAnnotation) {
  
  filename <- featDimList$filename[featDimList$subPrimaryAnnotation == subPrimaryAnnotation]
  resolution <- featDimList$resolution[featDimList$subPrimaryAnnotation == subPrimaryAnnotation]
  filtered_resolutions <- resolutions[resolutions != resolution]

  so <- readRDS(file.path(rdataDir, filename))
  so@meta.data$seurat_clusters <- factor(
    so@meta.data$seurat_clusters,
    levels = seq_along(unique(so@meta.data$seurat_clusters))
  )


  m_filtered <- so@meta.data %>%
    select(-contains(c(paste0("SCT_snn_res.", filtered_resolutions))))
  m_filtered$seurat_clusters <- as.factor(m_filtered[[paste0("SCT_snn_res.", resolution)]])
  m_filtered$Secondary_SeuratClusters <- m_filtered$seurat_clusters
  m_filtered$barcode <- rownames(m_filtered)

  m <- m_filtered %>% left_join(metadata, by = "barcode")
  rownames(m) <- m$barcode
  m$barcode <- NULL
  head(m)
  so@meta.data <- m
  Idents(so) <- "Secondary_SeuratClusters"

  p01 <- DimPlot(so, reduction = "umap_HM_02", label = TRUE, raster = TRUE)
  p02 <- DimPlot(so, reduction = "umap_HM_02", label = TRUE, group.by = "Condition", raster = TRUE)
  p03 <- DimPlot(so, reduction = "umap_HM_02", label = TRUE, group.by = "Project", raster = TRUE)
  p04 <- DimPlot(so, reduction = "umap_HM_02", label = TRUE, group.by = ref_annotation, raster = TRUE)

  p1 <- FeaturePlot(so, features = c("nCount_RNA", "nFeature_RNA"), reduction = "umap_HM_02", raster = TRUE)
  p2 <- FeaturePlot(so, features = mark1$gene, cols = c("grey", "red"), reduction = "umap_HM_02", raster = TRUE)
  p3 <- DotPlot(so, features = mark1$gene) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  p0 <- p01 + p04 / p02 + p03
  p4 <- p1 / p3

  markers <-  mark3[mark3$Category == subPrimaryAnnotation, ]$Genes
  # markers
  p5 <- DotPlot(so, features = markers) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  p6 <- FeaturePlot(so, features = markers, cols = c("grey", "red"), reduction = "umap_HM_02", raster = TRUE)

  return(list(so, p0, p2, p4, p5, p6))

}

plotAnalysis <- function(subPrimaryAnnotation){
  p1 <- resultsList[[2]]
  ggsave(
    plot = p1,
    file.path(outDir, paste0(subPrimaryAnnotation, "_withoutRN01_Dimplot.pdf")),
    width = 40, height = 20
    )
  p2 <- resultsList[[3]]
  ggsave(
    plot = p2,
    file.path(outDir, paste0(subPrimaryAnnotation, "_withoutRN01_Featureplot.pdf")),
    width = 20, height = 20
    )
  p3 <- resultsList[[4]]
  ggsave(
    plot = p3,
    file.path(outDir, paste0(subPrimaryAnnotation, "_withoutRN01_Dotplot.pdf")),
    width = 25, height = 20
    )

  p4 <- resultsList[[5]]
  ggsave(
    plot = p4,
    file.path(outDir, paste0(subPrimaryAnnotation, "_withoutRN01_Dotplot_canonical.pdf")),
    width = 25, height = 20
    )
  p5 <- resultsList[[6]]
  ggsave(
    plot = p5,
    file.path(outDir, paste0(subPrimaryAnnotation, "_withoutRN01_Featureplot_canonical.pdf")),
    width = 20, height = 20
    )
  
}


saveAnalysis <- function(subPrimaryAnnotation){
 
  so$SecondaryAnnotation <- unlist(annotations_2nd[so$Secondary_SeuratClusters])

  p6_1 <- DimPlot(
    so, reduction = "umap_HM_02",
    group.by = "SecondaryAnnotation",
    label = TRUE, raster = TRUE
    )
  p6_2 <- DimPlot(
    so, reduction = "umap_HM_02",
    group.by = "Project",
    label = TRUE, raster = TRUE
    )
  p6 <- p6_1 + p6_2
  ggsave(
    file.path(outDir, paste0(subPrimaryAnnotation, "_withoutRN01_SecondaryAnnotation.pdf")),
    plot = p6,
    width = 30, height = 20
    )

  # Save RData
  saveRDS(
    so,
    file = file.path(rdataDir, paste0(subPrimaryAnnotation, "_so_withoutRN01_SecondaryAnnotation.rds")))

  # Save Singlets
  so.singlet <- subset(so, subset = SecondaryAnnotation == "Doublet", invert = TRUE)

  saveRDS(
    so.singlet,
    file = file.path(rdataDir, paste0(subPrimaryAnnotation, "_so_withoutRN01_SecondaryAnnotation_singlet.rds")))
  
}


#' 
#' ### 01: Visualization of the subPrimaryAnnotation clusters
## -----------------------------------------------------------------------------
subsets <- read.table(file = file.path(rdataDir, paste0(tissue, "_subsets.txt")), header = FALSE, stringsAsFactors = FALSE)$V1
resultsLists <- list()

for (i in seq_along(subsets)){
  subPrimaryAnnotation <- subsets[[i]]

  resultsList <- runAnalysis(subPrimaryAnnotation)
  so <- resultsList[[1]]
  plotAnalysis(subPrimaryAnnotation)

  resultsLists[[i]] <- resultsList
  names(resultsLists)[[i]] <- subPrimaryAnnotation

}

#' 
#' ### 01: Astrocytes
## -----------------------------------------------------------------------------
subPrimaryAnnotation <- "Astrocytes"

resultsList <- resultsLists[[subPrimaryAnnotation]]
so <- resultsList[[1]]

annotations_2nd <- list(
  "1" = "Astrocytes",
  "2" = "Doublet", # Oligos
  "3" = "Astrocytes",
  "4" = "Astrocytes",
  "5" = "Astrocytes",
  "6" = "Doublet", # Micros
  "7" = "Doublet" # Opcs
)

# Save RData with Secondary annotation
saveAnalysis(subPrimaryAnnotation)

#' 
#' ### 02: Micros
## -----------------------------------------------------------------------------
subPrimaryAnnotation <- "Micros"

resultsList <- resultsLists[[subPrimaryAnnotation]]
so <- resultsList[[1]]

annotations_2nd <- list(
  "1" = "Micros",
  "2" = "Micros", 
  "3" = "Micros",
  "4" = "Doublet", # Oligos
  "5" = "Doublet", # Astrocytes
  "6" = "Doublet", # Oligos
  "7" = "Micros", 
  "8" = "Macrophages",
  "9" = "Micros", 
  "10" = "Doublet" # Opcs
)

# Save RData with Secondary annotation
saveAnalysis(subPrimaryAnnotation)


#' 
#' 
#' ### 03: Others
## -----------------------------------------------------------------------------
subPrimaryAnnotation <- "Others"

resultsList <- resultsLists[[subPrimaryAnnotation]]
so <- resultsList[[1]]

annotations_2nd <- list(
  "1" = "Endothelial",
  "2" = "Lymphocytes",
  "3" = "Meninges",
  "4" = "Doublet", # Oligos
  "5" = "Pericytes",
  "6" = "Meninges", 
  "7" = "Ependymal", 
  "8" = "Schwann"
)

# Save RData with Secondary annotation
saveAnalysis(subPrimaryAnnotation)

#' 
#' ### 04: Oligodendrocytes
## -----------------------------------------------------------------------------
subPrimaryAnnotation <- "Oligos"

resultsList <- resultsLists[[subPrimaryAnnotation]]
so <- resultsList[[1]]

annotations_2nd <- list(
  "1" = "Oligos",
  "2" = "Oligos",
  "3" = "Oligos",
  "4" = "Opc", 
  "5" = "Oligos",
  "6" = "Doublet", # Astrocytes
  "7" = "Doublet", # Micros
  "8" = "Oligos",  
  "9" = "Doublet", # Astrocytes
  "10" = "OligoProjenitor",
  "11" = "Doublet", # Micros/alpha
  "12" = "Doublet" # Micros
)

saveAnalysis(subPrimaryAnnotation)

#' 
#' ### 05: Neurons
## -----------------------------------------------------------------------------
subPrimaryAnnotation <- "Neurons"

resultsList <- resultsLists[[subPrimaryAnnotation]]
so <- resultsList[[1]]

annotations_2nd <- list(
  "1" = "Doublet", # Oligos
  "2" = "EX",
  "3" = "EX",
  "4" = "EX",
  "5" = "EX",
  "6" = "INH",
  "7" = "EX",
  "8" = "EX",
  "9" = "INH",
  "10" = "Doublet", # Astrocytes
  "11" = "EX",
  "12" = "INH",
  "13" = "EX",
  "14" = "EX",
  "15" = "INH",
  "16" = "INH",
  "17" = "EX",
  "18" = "MotorNeurons",
  "19" = "Doublet",
  "20" = "EX",
  "21" = "INH",
  "22" = "Doublet" # Micros
)

saveAnalysis(subPrimaryAnnotation)

#' 
#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../03_SecondaryAnnotation.Rmd"), documentation = 2)

