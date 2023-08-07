#' ---
#' title: "03 - Integration using Harmony considering batch effect"
#' author: "ertakeuchi"
#' date: "05 August 2023"
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
#' ## 0-0. Load libraries
## ---- include=FALSE-----------------------------------------------------------
library("Seurat")
library("tidyverse")
library("harmony")
library("reticulate")
use_python("/usr/bin/python3")
library("patchwork")
library("leiden")
set.seed(1234)

#' 
#' 
#' ## 0-1. Set a common directory
## -----------------------------------------------------------------------------
homeDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq")
NGSDir <- file.path(homeDir, "../../NGS_original")
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

# tissue <- "spinalcord"
# directory <- "230805_experiment2/out_2"
outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)

# read sample information
meta_datas <- read.table(file.path(metadataDir, "SampleInformation.tsv"), header = TRUE, sep = "\t")
meta_data <- meta_datas[meta_datas$Tissue == tissue,]
# meta_data <- meta_datas[,c("Project", "SubjectID", "SampleID", "Tissue", "Condition", "Sex", "Age", "Onset")]
print(meta_data)

#' 
#' 
#' ################################################################################
#' # Preprocessing 03: Integration using Harmony considering batch effect
#' ################################################################################
#' 
#' ## 1. Set environment
## -----------------------------------------------------------------------------
# Load a merged but unintegrated Seurat object  
so_list <- readRDS(file.path(rdataDir, paste0(tissue, "_singlet_list_ncells0.rds")))
# Set SCT assay
so_list <- lapply(so_list, function(x) {DefaultAssay(x) <- "SCT"; x})

# 01: Set environment
# |
# |-- 02: Integration test
# |   |
# |   |-- 01: Preprocessing
# |   |
# |   |-- 02: Add metadata
# |   |
# |   |-- 03: Run PCA
# |   |
# |   |-- 04: RUN Harmony
# |   |
# |   |-- 05: Run UMAP
# |   |
# |   |-- 06: Find neighbors
# |   |
# |   |-- 07: Find clusters
# |   |
# |   |-- 08: Visualize results and save figures
# |
# |-- 03: Integration specified parameters
#     |
#     |-- 01:Save RData

#' 
#' 
#' ## 2. Integration test
## -----------------------------------------------------------------------------
nfeat_sct <- c(3000, 5000)
pca_dim <- c(15, 30, 50)

for (feature in nfeat_sct){

  # 01: Preprocessing
  so.features <- SelectIntegrationFeatures(
    object.list = so_list,
    nfeatures = feature, # default = 2000
    fvf.nfeatures = feature # default = 2000
    )

  so_list <- PrepSCTIntegration(
    object.list = so_list,
    anchor.features = so.features
    )

  # Merge Seurat objects
  prefixes <- lapply(
    so_list, function(x) unique(x@meta.data$Project)
  ) %>% unlist()
  # prefixes

  so.m <- merge(
    so_list[[1]],
    y = so_list[2:length(so_list)],
    project = tissue,
    add.cell.ids = prefixes,
    merge.data = TRUE
    )

  VariableFeatures(so.m) <- so.features 
  write.csv((so.features), file = file.path(outDir, "so_features.csv"))

  # 02: Add metadata
  meta_data$SubjectID <- as.factor(meta_data$SubjectID)

  som_meta_data <- so.m@meta.data[,c("SubjectID", "Project")] %>%
    rownames_to_column("barcodes")
  som_meta_data$SampleID <- paste(som_meta_data$Project, som_meta_data$SubjectID, sep = "_")

  merged_meta_data <- left_join(x = som_meta_data, y = meta_data, by = "SampleID") %>% 
    column_to_rownames("barcodes") %>% 
    select(-SampleID)

  so.m <- AddMetaData(object = so.m, metadata = merged_meta_data)

    for (dim in pca_dim){

      # 03: Run PCA
      so.m <- RunPCA(
        object = so.m,
        assay = "SCT",
        features = so.features,
        npcs = dim # default = 50
        )

      # visualization of PCA results as DimHeatmap
      for (i in seq(from = 1, to = dim, by = 10)){
        # print(i)

        if(i + 9 <= dim){
          p1 <- DimHeatmap(
            so.m,
            dims = i:(i + 9),
            cells = 500,
            raster = TRUE,
            balanced = TRUE,
            ncol = 5
            )
          #  p1
          ggsave(file.path(
            outDir, paste0(tissue, "_HVG", feature, "_PCA_", i, "_", (i+9), "DimHeat.png")),
            plot = p1,
            height = 20, width = 20)
          } else{

          }
          
        }
      # visualization of PCA results as ElbowPlot
      p2 <- ElbowPlot(
        so.m,
        ndims = dim # default = 20
        )
      # p2
      ggsave(file.path(
        outDir, paste0(tissue, "_HVG", feature, "_PCA_", dim, "Elbow.png")),
        plot = p2,
        height = 20, width = 20
        )
      
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
          reduction.name = "umap_HM_01"
        ) %>%
        
      # 06: Find neighbors
        FindNeighbors(
          assay = "SCT", 
          reduction = "harmony",
          dims = 1:dim
          ) 

      # 07: Find clusters
        resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)

        for (resolution in resolutions){

          so.m <- FindClusters(
            object = so.m,
            resolution = resolution,
            algorithm = 4, # default = 3, 4 = leiden
            method = "igraph"
            )

        # 08: Visualize results and save figures
          col <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

          p1 <- VlnPlot(object = so.m, group.by = "Project", features = col, pt.size = 0.1, ncol = 3) &
                  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) 
          p2 <- FeatureScatter(so.m, col[2], col[1], group.by = "Project", pt.size = 0.5)
          p3 <- DimPlot(so.m, reduction = "umap_HM_01", repel = TRUE, label = TRUE, raster = TRUE, shuffle = TRUE)
          p4 <- DimPlot(so.m, reduction = "umap_HM_01", repel = TRUE, label = TRUE, split.by = "Project", raster = TRUE, shuffle = TRUE)
          p5 <- p1 / p2
          ggsave(file.path(outDir, paste0(tissue, "_", resolution, "_VlnScat_HVG", feature, "_pca", dim, ".png")), plot = p5, height = 20, width = 20)
          ggsave(file.path(outDir, paste0(tissue, "_", resolution, "_DimPlot_HVG", feature, "_pca", dim, ".png")), plot = p3, height = 20, width = 20)
          ggsave(file.path(outDir, paste0(tissue, "_", resolution, "_DimPlot_HVG", feature, "_pca", dim, "_splitbyProject.png")), plot = p4, height = 20, width = 20)


      # 09: Save RData
          saveRDS(
            object = so.m,
            file = file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, ".rds"))
          )
        }
    }

  }



#' 
#' ## 3. Integration specified parameters
#' 
## -----------------------------------------------------------------------------
# # define parameters
# feature <- 1000
# dim <- 15
# resolution <- 0.5

# # 01: Preprocessing
# so.features <- SelectIntegrationFeatures(
#   object.list = so_list,
#   nfeatures = feature, # default = 2000
#   fvf.nfeatures = feature # default = 2000
#   )

# so_list <- PrepSCTIntegration(
#   object.list = so_list,
#   anchor.features = so.features
#   )

# # Merge Seurat objects
# prefixes <- lapply(
#   so_list, function(x) unique(x@meta.data$Project)
# ) %>% unlist()
# # prefixes

# so.m <- merge(
#   so_list[[1]],
#   y = so_list[2:length(so_list)],
#   project = tissue,
#   add.cell.ids = prefixes,
#   merge.data = TRUE
#   )

# VariableFeatures(so.m) <- so.features 
# write.csv((so.features), file = file.path(outDir, "so_features.csv"))

# # 02: Add metadata
# meta_data$SubjectID <- as.factor(meta_data$SubjectID)

# som_meta_data <- so.m@meta.data[,c("SubjectID", "Project")] %>%
#   rownames_to_column("barcodes")
# som_meta_data$SampleID <- paste(som_meta_data$Project, som_meta_data$SubjectID, sep = "_")

# merged_meta_data <- left_join(x = som_meta_data, y = meta_data, by = "SampleID") %>% 
#   column_to_rownames("barcodes") %>% 
#   select(-SampleID)

# so.m <- AddMetaData(object = so.m, metadata = merged_meta_data)

# # 03: Run PCA
# so.m <- RunPCA(
#   object = so.m,
#   assay = "SCT",
#   features = so.features,
#   npcs = dim # default = 50
#   )
  
# # 04: RUN Harmony
# so.m <- so.m %>%
#   RunHarmony(
#     assay.use = "SCT", 
#     reduction = "pca", # default
#     dims.use = 1:dim, # default = all
#     group.by.vars = "SubjectID", 
#     plot_convergence = TRUE # default = FALSE
#   ) %>%

# # 06: Run UMAP
#   RunUMAP(
#     assay = "SCT",
#     reduction = "harmony",
#     dims = 1:dim,
#     reduction.name = "umap_HM_01"
#   ) %>%
  
# # 05: Find neighbors
#   FindNeighbors(
#     assay = "SCT", 
#     reduction = "harmony",
#     dims = 1:dim
#     ) %>%
  
# # 06: Find clusters
#   FindClusters(
#     resolution = resolution,
#     algorithm = 4, # default = 3, 4 = leiden
#     method = "igraph"
#     )

# # # 07: Save RData
# saveRDS(
#   object = so.m,
#   file = file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, ".rds"))
# )


#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../03_integrationHarmony.Rmd"), documentation = 2)

