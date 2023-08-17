#' ---
#' title: "03 - Embedding Public data using Symphony2"
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
#' # 1. Set environment
#' 
#' ## 0. Load libraries
## ---- include=FALSE-----------------------------------------------------------
library("Seurat")
library("SeuratDisk")
library("Signac")
library("tidyverse")
library("reticulate")
use_python("/usr/bin/python3")
library("patchwork")
library("gridExtra")
library("cowplot")
library("symphony")
library("magrittr")
library("data.table") 
library("leiden")
set.seed(1234)

#' 
#' 
#' ## 1. Set a common directory (FOR CONTAINER)
## -----------------------------------------------------------------------------
homeDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq")
metadataDir <- file.path(homeDir, "metadata")

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if the correct number of arguments is provided
if (length(args) != 4) {
  stop("Usage: Rscript test.R tissue output_directory", call. = FALSE)
}
# # Extract the tissue argument
# tissue <- args[1]
# directory <- args[2]
# gse <- args[3]
# feature <- args[4]

tissue <- "spinalcord"
directory <- "230816_experiment2/out_2"
gse <- "GSE190442"
outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)
NGSDir <- file.path(homeDir, "../../NGS_public", gse)
# feature <- 5000
# dim <- 50

# Call Symphony 
source(file.path(homeDir, 'R/Symphony/utils_seurat.R'))

#' 
#' ## 2. Contents
## -----------------------------------------------------------------------------
# 01: Set environment
# | 
# |-- 02: Build Reference
# |   |
# |   |-- 2-1: Load data for ref
# |   |
# |   |-- 2-2: Do seurat Analysis
# |   |       |-- 01: Create Seurat object from counts
# |   |       |-- 02: Run SCTransform
# |   |
# |   |-- 2-3: Make Symphony ref object
# |   
# |-- 03: MapQuery
# |   |
# |   |-- 03-1: Load query data and MapQuery
# |   |
# |   |-- 03-2: UMAP
# |   |
# |   |-- 03-3: Predict Clusters
# |
# |-- 04: MapQuery and predict clusters
# |
# |-- 05: Add reference annotation to query SeuratObject with Primary Annotation


#' 
#' ################################################################################
#' # Annotation 02: Symphony 
#' ################################################################################
#' 
#' # 02. Build Reference
#' 
#' ## 2-1. Load data for ref
#'  
## -----------------------------------------------------------------------------
# Load public count data
counts <- read.csv(
        file.path(NGSDir, paste0(gse, "_aggregated_counts_postqc.csv")),
        header = TRUE,
        row.names = 1
        )
metadata <- read.csv(
        file.path(NGSDir, paste0(gse, "_aggregated_metadata_postqc.csv")),
        row.names = 1
        )

# head(metadata)
# counts[1:3, 1:3]

# modify cell ids for creating seurat object
rownames(metadata) <- sub("-", ".", rownames(metadata))
head(metadata)


#' 
#' ## 2-2. Do seurat Analysis
## -----------------------------------------------------------------------------
.verbose <- TRUE

suppressWarnings({ # to suppress the SCTransform warning: "iteration limit reached"
        ref_obj <- CreateSeuratObject(
                counts = counts,
                meta.data = metadata) 
                # %>%
#                 SCTransform(
#                         assay = "RNA",
#                         vst.flavor = "v2",
#                         method = "glmGamPoi",
#                         verbose = .verbose
#                 )
# }) 

        ref_list <- SplitObject(ref_obj, split.by = "sample")

                ## 03. SCTransform
        sct_list <- lapply(ref_list, FUN = function(x){
                x <- SCTransform(
                        x,
                        assay = "RNA",
                        vst.flavor = "v2",
                        method = "glmGamPoi",
                        verbose = .verbose
                )
                return(x)
        })
})

so.features <- SelectIntegrationFeatures(
        object.list = sct_list
        ) # default HVG = 2000

sct_list <- PrepSCTIntegration(
        object.list = sct_list,
        anchor.features = so.features
        )

ref_obj <- merge(
        sct_list[[1]],
        y = sct_list[2:length(sct_list)],
        project = gse,
        merge.data = TRUE
        )

VariableFeatures(ref_obj) <- so.features 

# ref_obj
dim <- 30
feature <- 2000
.verbose <- FALSE
# 03: Run PCA
ref_obj <- ref_obj %>% 
        ScaleData(verbose = .verbose) %>%
        RunPCA(
                assay = "SCT",
                npcs = dim,
                features = so.features,
                verbose = .verbose
                ) 

ref_obj@meta.data$sample <- ref_obj@meta.data$sample %>% as.character()

ref_obj <- ref_obj %>% 
        RunHarmony.Seurat(
                group.by.vars = "sample",
                assay.use = "SCT",
                reduction = "pca",
                dim.use = 1:dim,
                project.dim = FALSE,
                plot_convergence = TRUE,
                verbose = .verbose
                ) %>%
        FindNeighbors(
                assay = "SCT", 
                reduction = "harmony",
                dims = 1:dim,
                verbose = .verbose
                ) %>%
        FindClusters(
                resolution = 0.6,
                algorithm = 4, # default = 3, 4 = leiden
                method = "igraph",
                verbose = .verbose
                )

# UMAP
ref_obj[['umap']] <- RunUMAP2(
        Embeddings(ref_obj, "harmony")[, 1:dim], 
        assay = 'RNA',
        verbose = FALSE,
        umap.method = 'uwot', 
        return.model = TRUE
        )

options(repr.plot.height = 4, repr.plot.width = 6)
p <- DimPlot(
        ref_obj,
        reduction = 'umap',
        group.by = 'seurat_clusters',
        shuffle = TRUE
        )

ggsave(
        file.path(outDir, paste0(gse, "_varsSample_SymphonyRefUmap.pdf")),
        plot = p,
        width = 6,
        height = 4
        )


#' 
#' 
#' ## 2-3. Make Symphony ref object
#' 
## -----------------------------------------------------------------------------
ref <- buildReferenceFromSeurat(
        ref_obj, 
        assay = "SCT",
        verbose = TRUE, 
        save_umap = TRUE, 
        save_uwot_path = "cache_symphony_sct.uwot"
        )

ref$normalization_method = "SCTransform"

saveRDS(ref_obj, file = file.path(outDir, paste0(gse, "_ref_obj.rds")))
# ref_obj <- readRDS(file = file.path(outDir, "../", paste0( gse, "_ref_obj.rds")))

#' 
#' 
#' # 03. MapQuery
#' 
#' ## 3-1. Load query data and MapQuery
## -----------------------------------------------------------------------------
# # Load query data
# soSCTList <- readRDS(file = file.path(rdataDir, "soSCTList_ForSymphony.rds"))
# # names(soSCTList)
# sct_list <- soSCTList[["withoutRN01"]]

# rm(soSCTList)

# so.features <- SelectIntegrationFeatures(
#         object.list = sct_list
#         ) # default HVG = 2000

# sct_list <- PrepSCTIntegration(
#         object.list = sct_list,
#         anchor.features = so.features
#         )

# query_obj <- merge(
#         sct_list[[1]],
#         y = sct_list[2:length(sct_list)],
#         project = gse,
#         merge.data = TRUE
#         )

# saveRDS(query_obj, file = file.path(rdataDir, "query_obj.rds"))

query_obj <- readRDS(file = file.path(rdataDir, "query_obj.rds"))

# Map query
query <- mapQuery(
        exp_query = query_obj@assays$SCT@scale.data, 
        metadata_query = query_obj@meta.data,
        ref_obj = ref, 
        vars = c("SubjectID"), 
        do_normalize = FALSE,
        return_type = "Seurat"
        )


#' 
#' ## 3-2. UMAP
## -----------------------------------------------------------------------------
options(repr.plot.height = 4, repr.plot.width = 10)
p1 <- DimPlot(
        ref_obj,
        reduction = 'umap',
        shuffle = TRUE
        ) + labs(title = 'Original Reference (Clusters)') 

p2 <- DimPlot(
        query,
        reduction = 'umap',
        group.by = "Project",
        shuffle = TRUE
        ) + labs(title = 'Mapped Query (Donors)')

p <- p1 + p2

ggsave(
        file.path(outDir, "refMappingUMAP_RefMapping.pdf"),
        plot = p,
        width = 10,
        height = 4)

#' 
#' ## 3-3. Predict Clusters
## -----------------------------------------------------------------------------
query <- knnPredict.Seurat(query, ref, 'seurat_clusters')

options(repr.plot.height = 4, repr.plot.width = 10)
p1 <- DimPlot(
        ref_obj,
        reduction = 'umap',
        group.by = 'seurat_clusters',
        shuffle = TRUE,
        raster = TRUE
        ) + labs(title = 'Original Reference (Clusters)') 
p2 <- DimPlot(
        query,
        reduction = 'umap',
        group.by = 'seurat_clusters',
        shuffle = TRUE,
        raster = TRUE
        ) + labs(title = 'Mapped Query (Predicted Clusters)')

p <- p1 + p2
# p
ggsave(
        file.path(outDir, "refMappingPredictCluster_RefMapping.pdf"),
        plot = p,
        width = 10,
        height = 4
)

#' 
#' # 04. MapQuery and predict clusters
## ---- fig.height=5, fig.width=30----------------------------------------------

## Predict clusters
clusterL1 <- "top_level_annotation"
clusterL2 <- "subtype_annotation"

queryL1 <- knnPredict.Seurat(
        query_obj = query,
        ref_obj = ref,
        label_transfer = clusterL1
        )

queryL2 <- knnPredict.Seurat(
        query_obj = query,
        ref_obj = ref,
        label_transfer = clusterL2
        )

# UMAP by predict clusters
options(repr.plot.height = 10, repr.plot.width = 20)

p1 <- DimPlot(
        queryL1,
        reduction = "umap",
        group.by = clusterL1,
        shuffle = TRUE, 
        label = TRUE, 
        repel = TRUE,
        raster = TRUE
        ) +
        labs(title = paste0("Mapped Query Predicted", clusterL1))

p2 <- DimPlot(
        queryL2,
        reduction = "umap",
        group.by = clusterL2,
        shuffle = TRUE, 
        label = TRUE, 
        repel = TRUE,
        raster = TRUE
        ) +
        labs(title = paste0("Mapped Query Predicted", clusterL2))

ggsave(
        file.path(outDir, paste0(clusterL1, "_refMappingPredictCluster_RefMapping.pdf")),
        plot = p1,
        width = 30,
        height = 20
        )
ggsave(
        file.path(outDir, paste0(clusterL2, "_refMappingPredictCluster_RefMapping.pdf")),
        plot = p2,
        width = 30,
        height = 20
        )

p3 <- DimPlot(
        ref_obj,
        reduction = 'umap',
        group.by = clusterL1,
        label = TRUE,
        repel = TRUE,
        shuffle = TRUE,
        raster = TRUE
        ) +
        labs(title = paste0("Original Reference",  clusterL1))

p4 <- DimPlot(
        ref_obj,
        reduction = 'umap',
        group.by = clusterL2,
        label = TRUE,
        repel = TRUE,
        shuffle = TRUE,
        raster = TRUE
        ) +
        labs(title = paste0("Original Reference",  clusterL2))

ggsave(
        file.path(outDir, paste0(clusterL1, "_RefMapping.pdf")),
        plot = p3,
        width = 30,
        height = 20
        )
ggsave(
        file.path(outDir, paste0(clusterL2, "_RefMapping.pdf")),
        plot = p4,
        width = 30,
        height = 20
        )

d1 <- queryL1@meta.data[c(clusterL1, paste0(clusterL1, "_prob"), "Project")]
d2 <- queryL2@meta.data[c(clusterL2, paste0(clusterL2, "_prob"))]
df <- cbind(d1, d2)

df$barcode <- rownames(df)
df$barcode <- strsplit(df$barcode, "_") %>% sapply("[", 1)
df$barcode <- paste0(df$Project, "_", df$barcode)
head(df)
rownames(df) <- df$barcode
df$barcode <- NULL

# output files
write.csv(
        df,
        file = file.path(outDir, "refDataRefMapping.csv"),
        row.names = TRUE,
        quote = FALSE
        )

# query

#' 
#' # 05. Add reference annotation to query SeuratObject with 1st annotation
## -----------------------------------------------------------------------------
# Load reference mapping as metadata
metadata <- read.csv(
        file = file.path(outDir, "refDataRefMapping.csv"), 
        sep = ",",
        header = TRUE,
        row.names = 1
        )
# head(metadata)

feature <- 3000
dim <- 50

# Load local Seurat object 
so <- readRDS(
        file = file.path(
                rdataDir,
                paste0(
                        tissue, "_som_HVG", feature, "_pca", dim, "_withoutRN01_PrimaryAnnotation.rds"
                        )
                )
        )

head(metadata)

# Add metadata to SeuratObject with 1st annotation
so <- AddMetaData(so, metadata = metadata)
head(so@meta.data)

saveRDS(so, file = file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_PrimaryAnnotationSymphony.rds")))

#' 
#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../03_Symphony2.Rmd"), documentation = 2)

