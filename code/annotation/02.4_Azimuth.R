#' ---
#' title: "02 - Preparing Secondary Annotation, Azimuth embedding"
#' author: "ertakeuchi"
#' date: "17 August 2023"
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
library("gridExtra")
library("leiden")
library("Azimuth")
set.seed(1234)


#' 
#' 
#' ## 0-1. Set a common directory
## -----------------------------------------------------------------------------
homeDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq")
metadataDir <- file.path(homeDir, "metadata")

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("Usage: Rscript test.R tissue output_directory", call. = FALSE)
}
# Extract the tissue argument
tissue <- args[1]
directory <- args[2]

tissue <- "brain"
directory <- "230817_experiment3/out_3"
# feature <- 3000
# dim <- 50
outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)


#' 
#' 
#' ################################################################################
#' # Annotation 02.1: Emberdding Azimuth 
#' ################################################################################
#' 
## -----------------------------------------------------------------------------
# 01: Set environment
# |
# |-- 02: AzimuthEmbedding
# |   |
# |   |-- 01: Load Secondary Annotation objects, and split by Project
# |   |
# |   |-- 02: SCTransform
# |   |
# |   |-- 03: Save RData
# |   



#' 
#' 
## -----------------------------------------------------------------------------
# Download from https://zenodo.org/record/4546932#.Y83g8OzP1f0
reference <- LoadReference(path = file.path(homeDir, "../../reference/Azimuth/Human_M1_cortex"))
# following URL was expired!!
# reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_motorcortex")

#' 
#' ## load the query SeuratObject
## -----------------------------------------------------------------------------
feature <- 3000
dim <- 50

query_obj <- readRDS(file.path(outDir, "query_obj.rds" ))

#' 
#' ## Find anchors 
## -----------------------------------------------------------------------------
# find anchors
anchors <- FindTransferAnchors(
  reference = reference$map,
  query = query_obj,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT", 
#   reduction = "pcaproject",
  features = intersect(
        rownames(x = reference$map),
        VariableFeatures(object = query_obj)
        ),
  dims = 1:50,
  n.trees = 50, # default = 50
  mapping.score.k = 100
  )

#' 
#' ## Transfer cell type labels and impute gene expression
## -----------------------------------------------------------------------------

clusterL1 <- "subclass"
clusterL2 <- "cluster"

refdata <- lapply(X = clusterL1, function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- clusterL1
if (FALSE) {
  refdata[["impADT"]] <- GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )
}

refdata_cluster <- lapply(X = clusterL2, function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata_cluster) <- clusterL2
if (FALSE) {
  refdata_cluster[["impADT"]] <- GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )
}

query <- TransferData(
  reference = reference$map,
  query = query_obj,
  dims = 1:50,
  anchorset = anchors,
  refdata = c(refdata, refdata_cluster),
  n.trees = 20,
  store.weights = TRUE
)



#' 
#' ## Integrate embeddings and calculate mapping scores
## -----------------------------------------------------------------------------
# Calculate the embedding of the query on the reference
query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference$map,
  query = query, 
  new.reduction.name = "integrated_dr",
  reductions = "pcaproject",
  reuse.weights.matrix = TRUE
)

# Calculate the query neighbors 
query[["query_ref.nn"]] <- FindNeighbors(
  object = Embeddings(reference$map[["refDR"]]),
  query = Embeddings(query[["integrated_dr"]]),
  return.neighbor = TRUE,
  l2.norm = TRUE
)

query <- Azimuth:::NNTransform(
        object = query,
        meta.data = reference$map[[]]
        )

# Project the query to the reference UMAP.
query[["proj.umap"]] <- RunUMAP(
  object = query[["query_ref.nn"]],
  reduction.model = reference$map[["refUMAP"]],
  reduction.key = 'UMAP_'
)

# Calculate mapping score and add to metadata
query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = anchors),
  col.name = "mapping.score"
)

#' 
#' ## Visualization
## -----------------------------------------------------------------------------

clusterL1 <- "predicted.subclass"
clusterL2 <- "predicted.cluster"
reduction_name <- "proj.umap"

  
p1 <- DimPlot(
        query,
        reduction = reduction_name,
        group.by = clusterL1,
        shuffle = TRUE, 
        label = TRUE, 
        repel = TRUE,
        raster = TRUE
                ) +
        labs(title = "Mapped query (Predicted clusterL1)")
p2 <- DimPlot(
        query,
        reduction = reduction_name,
        group.by = clusterL2,
        shuffle = TRUE, 
        label = TRUE, 
        repel = TRUE,
        raster = TRUE
                ) +
        labs(title =  "Mapped Query (Predicted clusterL2)")

p <- p1 / p2

ggsave(
        file.path(outDir, "AzimuthMappingtoQuery.pdf"),
        plot =  p,
        width = 40,
        height = 20
        )

d1 <- query@meta.data[c(clusterL1, paste0(clusterL1, ".score"))]
d2 <- query@meta.data[c(clusterL2, paste0(clusterL2, ".score"))]
df <- cbind(d1, d2)
# head(df)
# df$barcode <- rownames(df)
# df$barcode <- strsplit(df$barcode, "_") %>% sapply("[", 1)
# df$barcode <- paste0(df$Project, "_", df$barcode)
# head(df)
# rownames(df) <- df$barcode
# df$barcode <- NULL

# output files
write.csv(
        df,
        file = file.path(outDir, "refDataRefMapping_Azimuth.csv"),
        row.names = TRUE,
        quote = FALSE
        )

#' 
#' 
#' ## Add reference annotation to query with PrimaryAnnotation
#' 
## -----------------------------------------------------------------------------
# Load reference mapping as metadata
metadata <- read.csv(
        file = file.path(outDir, "refDataRefMapping_Azimuth.csv"), 
        sep = ",",
        header = TRUE,
        row.names = 1
        )

feature <- 3000
dim <- 50

so <- readRDS(file.path(
        rdataDir,
        paste0(
                tissue, "_som_HVG", feature, "_pca", dim, "_withoutRN01_PrimaryAnnotation.rds"
                )
        )
        )

# Add metadata to SeuratObject with 1st annotation
so <- AddMetaData(so, metadata = metadata)
head(so@meta.data)

saveRDS(so, file = file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_PrimaryAnnotationAzimuth.rds")))

#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
## -----------------------------------------------------------------------------
# check

#' 
