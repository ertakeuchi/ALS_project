#' ---
#' title: "02 - Preparing for Secondary annotation"
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
library("harmony")
library("leiden")
set.seed(1234)


#' 
#' 
#' ## 0-1. Set a common directory (FOR CONTAINER)
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

# tissue <- "brain"
# directory <- "230816_experiment4/out_4"
feature <- 3000
dim <- 50
outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)


#' 
#' ################################################################################
#' # Annotation 02: Secondary annotation 
#' ################################################################################
#' 
#' # 02. Devide to primary annotation subsets
#' 
## -----------------------------------------------------------------------------
# 01: Set environment
# |
# | 
# |-- 02: Devide to Primary Annotation subsets (cell type)
# |   |
# |   |-- 01: Load Seurat object with Primary Annotation
# |   |
# |   |-- 02: Devide to subsets, and save
# | 
# |-- 03: SCTransform to each subset
# |   |
# |   |-- 01: Load subset Seurat object, and split by Project
# |   |
# |   |-- 02: SCTransform
# |   |
# |   |-- 03: Save RData
# |   
# |-- 04: Integration using Harmony to each subset
# |   |
# |   |-- 01: Preprocessing 
# |   |   
# |   |-- 02: Merge Seurat objects
# |   |
# |   |-- 03: Run PCA
# |   |
# |   |-- 04: Run Harmony
# |   |
# |   |-- 05: Run UMAP
# |   |
# |   |-- 06: Find Neighbors
# |   |
# |   |-- 07: Find Clusters
# |   |
# |   |-- 08: Visualize results and save figures
# |   |
# |   |-- 09: Save RData


#' 
#' # 02: Devide to Primary Annotation subsets
#' 
#' ## 01: Load Seurat object with Primary Annotation
#' 
## -----------------------------------------------------------------------------
# so <- readRDS(
#   file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca",  dim, "_withoutRN01_PrimaryAnnotation.rds")))

# #' 
# #' ## 02: Devide to subsets, and save
# ## -----------------------------------------------------------------------------
# subsets <- unique(so$subPrimaryAnnotation)
# # subsets
# write.table(subsets, file = file.path(rdataDir, paste0(tissue, "_subsets.txt")), row.names = FALSE, col.names = FALSE, quote = FALSE)

# for (s in subsets){

#   so.sub <- subset(so, subset = subPrimaryAnnotation == s)

#   saveRDS(
#     so.sub, 
#     file = file.path(
#       rdataDir,
#       paste0(s, "_so_withoutRN01_subPrimaryAnnotation.rds")
#       )
#     )
# }

#' 
#' 
#' # 03: SCTransform to each subset
#' 
#' ## 01: Split Objects by Project
## -----------------------------------------------------------------------------
# Get a Seurat object file list
# file_list <- list.files(
#   rdataDir,
#   pattern = "_so_withoutRN01_subPrimaryAnnotation.rds",
#   full.name = FALSE
#   )
# # file_list

# # Get a project name list
# subName_list <- lapply(file_list, function(file_name) {
#   prefix <- gsub("_so_withoutRN01_subPrimaryAnnotation.rds", "", file_name)
#   return(prefix)
# })

# # Load seurat object, split by Project for each cell type
# so_list <- lapply(file_list, function(file_name) {
#   so <- readRDS(file.path(rdataDir, file_name))
#   return(so)
# })

# for (i in seq_along(so_list)){
#   # i = 1
#   obj <- so_list[[i]]

#   split_list <- SplitObject(
#     object = obj,
#     split.by = "Project"
#     )
#   subPrimaryAnnotation <- unique(obj@meta.data$subPrimaryAnnotation)

#  saveRDS(
#    split_list, 
#    file = file.path(
#      rdataDir,
#      paste0(subPrimaryAnnotation, "_withoutRN01_split_list.rds")
#      )
#    )

# } 



#' 
#' ## 02: SCTransform
## -----------------------------------------------------------------------------
# subsets <- read.table(file = file.path(rdataDir, paste0(tissue, "_subsets.txt")), header = FALSE, stringsAsFactors = FALSE)$V1
# print(subsets)

# for (i in seq_along(subsets)){
#   # i = 4
#   subPrimaryAnnotation <- subsets[i]
#   split_list <- readRDS(file.path(rdataDir, paste0(subPrimaryAnnotation, "_withoutRN01_split_list.rds")))

#   sct_list <- lapply(split_list, function(x){

#     x <- SCTransform(
#       x,
#       assay = "RNA",
#       vst.flavor = "v2",
#       method = "glmGamPoi",
#       vars.to.regress = c("percent.mt", "pct_chrY"),
#       verbose = FALSE)
#   })

#   saveRDS(
#     sct_list, 
#     file = file.path(
#       rdataDir,
#       paste0(subPrimaryAnnotation, "_withoutRN01_sct_list.rds")
#       )
#     )
# }
    

#' 
#' # 04: Integration using Harmony to each subset
#' 
#' 
## -----------------------------------------------------------------------------
subsets <- read.table(file = file.path(rdataDir, paste0(tissue, "_subsets.txt")), header = FALSE, stringsAsFactors = FALSE)$V1 

for (j in seq_along(subsets)){

  sct_list <- readRDS(file.path(rdataDir, paste0(subsets[j], "_withoutRN01_sct_list.rds")))
  subPrimaryAnnotation <- subsets[j]
  print(paste0(subPrimaryAnnotation, " is loaded."))

  nfeat_sct <- c(1000, 2000, 3000)
  pca_dim <- c(10, 15, 30)

  for (feature in nfeat_sct){
    # 01: Preprocessing 
    sosub_features <- SelectIntegrationFeatures(
      object.list = sct_list,
      nfeatures = feature, # default = 2000
      fvf.nfeatures = feature # default = 2000
      )

    sct_list <- PrepSCTIntegration(
      object.list = sct_list,
      anchor.features = sosub_features
      )

    # 02: Merge Seurat objects
    so.m <- merge(
      sct_list[[1]],
      y = sct_list[2:length(sct_list)],
      project = tissue,
      merge.data = TRUE
      )

    VariableFeatures(so.m) <- sosub_features 
    write.csv((sosub_features), file = file.path(outDir, paste0(subPrimaryAnnotation, "_HVG", feature, "_withoutRN01_soSubFeatures.csv")))

  for (dim in pca_dim){

      # 03: Run PCA
      so.m <- RunPCA(
        object = so.m,
        assay = "SCT",
        features = sosub_features,
        npcs = dim # default = 50
        )

      # # visualization of PCA results as DimHeatmap
      # for (i in seq(from = 1, to = dim, by = 10)){

      #   if(i + 9 <= dim){
      #     p1 <- DimHeatmap(
      #       so.m,
      #       dims = i:(i + 9),
      #       cells = 500,
      #       raster = TRUE,
      #       balanced = TRUE,
      #       ncol = 5
      #       )
      #     #  p1
      #     ggsave(file.path(
      #       outDir, paste0(subPrimaryAnnotation, "_HVG", feature, "_PCA_", i, "_", (i+9), "_withoutRN01_DimHeat.png")),
      #       plot = p1,
      #       height = 20, width = 20)
      #     } else{

      #     }
          
        # }
      # # visualization of PCA results as ElbowPlot
      # p2 <- ElbowPlot(
      #   so.m,
      #   ndims = dim # default = 20
      #   )
      # # p2
      # ggsave(file.path(
      #   outDir, paste0(subPrimaryAnnotation, "_HVG", feature, "_PCA_", dim, "_withoutRN01_Elbow.png")),
      #   plot = p2,
      #   height = 20, width = 20
      #   )
      print(paste0(subPrimaryAnnotation, " is runnning Harmony.")) 
      # 04: RUN Harmony
      so.m <- so.m %>%
        RunHarmony(
          assay.use = "SCT", 
          reduction = "pca", # default
          dims.use = 1:dim, # default = all
          group.by.vars = "SubjectID", 
          plot_convergence = TRUE # default = FALSE
        ) %>%

      # 05: Run UMAP
        RunUMAP(
          assay = "SCT",
          reduction = "harmony",
          dims = 1:dim,
          reduction.name = "umap_HM_02"
        ) %>%
        
      # 06: Find neighbors
        FindNeighbors(
          assay = "SCT", 
          reduction = "harmony",
          dims = 1:dim
          ) 

      # 07: Find clusters
        resolutions <- c(0.4, 0.6, 0.8)

        for (resolution in resolutions){

          print(paste0(subPrimaryAnnotation, " is runnning the leiden Clustering")) 
          so.m <- FindClusters(
            object = so.m,
            resolution = resolution,
            algorithm = 4, # default = 3, 4 = leiden
            method = "igraph"
            )

          # 08: Visualize results and save figures
          p1 <- DimPlot(
            so.m,
            reduction = "umap_HM_02",
            repel = TRUE,
            label = TRUE,
            raster = TRUE,
            shuffle = TRUE
            )

          ggsave(
            file.path(outDir, paste0(subPrimaryAnnotation, "_", resolution, "_DimPlot_HVG", feature, "_pca", dim, "_withoutRN01.png")),
            plot = p1,
            height = 20,
            width = 20
            )

          # 09: Save RData
          saveRDS(
            object = so.m,
            file = file.path(rdataDir, paste0(subPrimaryAnnotation, "_som_HVG", feature, "_pca", dim, "_withoutRN01_preSecondaryAnnotation.rds"))
          )
        }
    }
  }
}


#' 
#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../02_prepSecondaryAnnotation.Rmd"), documentation = 2)

