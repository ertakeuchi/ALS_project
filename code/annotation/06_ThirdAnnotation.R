#' ---
#' title: "06 - Third annotation"
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
# knitr::opts_chunk$set(echo = TRUE)

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
tissue <- args[1]
directory <- args[2]

# tissue <- "spinalcord"
# directory <- "230823_experiment4/out_4"
groupby_annotation <- "ThirdAnnotation"
ref_annotation <- "subtype_annotation"

outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)
resolutions <- c(0.2, 0.4, 0.6, 0.8)

#' 
#' ################################################################################
#' # Annotation 05: Third annotation 
#' ################################################################################
#' 
#' # 02. Third annotation
#' 
## -----------------------------------------------------------------------------
# 01: Set environment
# |
# |-- 02: Third annotation 
# |   |
# |   |-- 01: Load Seurat object with Secondary Annotation using metadata
# |   |  |
# |   |  |-- 01: Load marker gene set for each annotation
# |   |
# |   |-- 02: Add Third annotation
# |   |  |
# |   |  |-- 00: Set Secondary Annotation function
# |   |  |-- 01: Visualization of the preThirdAnnotation clusters
# |   |  |-- 02: Astrocytes
# |   |  |-- 03: Micros
# |   |  |-- 04: Oligodendrocytes
# |   |  |-- 05: Neurons
# |   |  |-- 06: Others
# |   |
# |   |-- 03: Add Third annotation to all cell clusters, and save RData


#' 
#' ## 01: Load Seurat object with Secondary Annotation using metadata
#' 
## -----------------------------------------------------------------------------
featDimList <- read.table(file = file.path(metadataDir, "featureDimList.tsv"), header = TRUE, stringsAsFactors = FALSE)
featDimList <- featDimList[featDimList$tissue == tissue & featDimList$levels == 3, ] # levels-3 is for Third annotation
head(featDimList)

#' 
#' ### 01: Load marker gene set for each annotation
## -----------------------------------------------------------------------------
markers_df <- read.table(file = file.path(metadataDir, "canonical_markers.tsv"), header = TRUE, stringsAsFactors = FALSE)
mark1 <- markers_df[markers_df$Level == "Top", ]
mark2 <- markers_df[markers_df$Level == paste0("Sub_", tissue), ]
mark3 <- markers_df[markers_df$Level == paste0("SubPrimary_", tissue), ]

#' 
#' 
#' ## 02: Add Third annotation
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
  m_filtered$Third_SeuratClusters <- m_filtered$seurat_clusters

  so@meta.data <- m_filtered
  Idents(so) <- "Third_SeuratClusters"

  p01 <- DimPlot(so, reduction = "umap_HM_03", label = TRUE, raster = TRUE)
  p02 <- DimPlot(so, reduction = "umap_HM_03", label = TRUE, group.by = "Condition", raster = TRUE)
  p03 <- DimPlot(so, reduction = "umap_HM_03", label = TRUE, group.by = "Project", raster = TRUE)
  p04 <- DimPlot(so, reduction = "umap_HM_03", label = TRUE, group.by = ref_annotation, raster = TRUE)

  p1 <- FeaturePlot(so, features = c("nCount_RNA", "nFeature_RNA"), reduction = "umap_HM_03", raster = TRUE)
  p2 <- FeaturePlot(so, features = mark1$Genes, cols = c("grey", "red"), reduction = "umap_HM_03", raster = TRUE)
  p3 <- DotPlot(so, features = mark1$Genes) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  p0 <- p01 + p04 / p02 + p03
  p4 <- p1 / p3

  markers <-  mark3[mark3$Category == subPrimaryAnnotation, ]$Genes
   # markers
  p5 <- DotPlot(so, features = markers) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  p6 <- FeaturePlot(so, features = markers, cols = c("grey", "red"), reduction = "umap_HM_03", raster = TRUE)

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
 
  so$ThirdAnnotation <- unlist(annotations_3rd[so$Third_SeuratClusters])
  head(so@meta.data)

  p6_1 <- DimPlot(
    so, reduction = "umap_HM_03",
    group.by = groupby_annotation,
    label = TRUE, raster = TRUE
    )
  p6_2 <- DimPlot(
    so, reduction = "umap_HM_03",
    group.by = "Project",
    label = TRUE, raster = TRUE
    )
  p6 <- p6_1 + p6_2
  ggsave(
    file.path(outDir, paste0(subPrimaryAnnotation, "_withoutRN01_ThirdAnnotation.pdf")),
    plot = p6,
    width = 30, height = 20
    )

  # Save RData
  saveRDS(
    so,
    file = file.path(rdataDir, paste0(subPrimaryAnnotation, "_so_withoutRN01_ThirdAnnotation.rds")))

  m <- so@meta.data
  m_filtered <- m[, c("SecondaryAnnotation", "Third_SeuratClusters", "ThirdAnnotation")]

  # Save Third Annotation data with cell-barcodes as csv
  write.table(
    m_filtered,
    file = file.path(outDir, paste0(subPrimaryAnnotation, "_withoutRN01_ThirdAnnotation.csv")),
    row.names = TRUE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
    )

# }


# findDEMarkers <- function(subPrimaryAnnotation){

  feature <- featDimList$feature[featDimList$subPrimaryAnnotation == subPrimaryAnnotation]

  so <- PrepSCTFindMarkers(so)

  feature_path <- file.path(outDir, "../../230817_experiment1/out_1", paste0(subPrimaryAnnotation, "_HVG", feature, "_second_withoutRN01_soSubFeatures.csv"))

  so_features <- read_csv(file = feature_path)
  so_features <- so_features[-1]  # "X"の列を削除

  colnames(so_features) <- "gene"
  # head(so_features)
  so_features <- so_features$gene

  so.markers <- FindAllMarkers(
    so,
    assay = "SCT",
    only.pos = TRUE,
    test.use= "wilcox", # default
    min.pct = 0.1, # default 
    logfc.threshold = 0.25, # default = 0.25,
    features = so_features
    )

  write.table(
    so.markers,
    file = file.path(outDir, paste0(subPrimaryAnnotation, "_ThirdAnnotation_AllMarkers.tsv")),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
    )

  Idents(so) <- "ThirdAnnotation"
  clusters <- unique(so@meta.data$ThirdAnnotation)

  for(i in seq_along(clusters)){

    ThirdAnnotation <- clusters[[i]]
    celltype_markers <- FindMarkers(
    so,
    ident.1 = ThirdAnnotation,
    min.pct = 0.1, # default 
    logfc.threshold = 0.25, # default = 0.25,
    test.use= "wilcox", # default
    only.pos = TRUE,
    features = so_features

    )

    write.table(
      celltype_markers,
      file = file.path(outDir, paste0(subPrimaryAnnotation, "_", ThirdAnnotation, "_Markers.tsv")),
      row.names = TRUE,
      col.names = TRUE,
      sep = "\t",
      quote = FALSE
    )

  }
  

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
#' ### 02: Astrocytes
## -----------------------------------------------------------------------------
subPrimaryAnnotation <- "Astrocytes"

resultsList <- resultsLists[[subPrimaryAnnotation]]
so <- resultsList[[1]]

annotations_3rd <- list(
  "1" = "Astro_1",
  "2" = "Astro_2", # Oligos
  "3" = "Astro_3",
  "4" = "Astro_1",
  "5" = "Astro_2",
  "6" = "Astro_1"
)

# Save RData with Secondary annotation
saveAnalysis(subPrimaryAnnotation)
# findDEMarkers(subPrimaryAnnotation)

#' 
#' ### 03: Micros
## -----------------------------------------------------------------------------
subPrimaryAnnotation <- "Micros"

resultsList <- resultsLists[[subPrimaryAnnotation]]
so <- resultsList[[1]]

annotations_3rd <- list(
  "1" = "Micro_1",
  "2" = "Micro_2", 
  "3" = "Micro_3",
  "4" = "Micro_4", 
  "5" = "Micro_5", 
  "6" = "Macro_1",
  "7" = "Micro_6"
)

# Save RData with Secondary annotation
saveAnalysis(subPrimaryAnnotation)
# findDEMarkers(subPrimaryAnnotation)

#' 
#' ### 04: Others
## -----------------------------------------------------------------------------
subPrimaryAnnotation <- "Others"

resultsList <- resultsLists[[subPrimaryAnnotation]]
so <- resultsList[[1]]

annotations_3rd <- list(
  "1" = "Endothelial_1",
  "2" = "Meninges_2",
  "3" = "Lymphocyte_1",
  "4" = "Pericyte_1", 
  "5" = "Ependymal_1",
  "6" = "Meninges_1"
)

# Save RData with Secondary annotation
saveAnalysis(subPrimaryAnnotation)
# findDEMarkers(subPrimaryAnnotation)

#' 
#' ### 05: Oligodendrocytes
## -----------------------------------------------------------------------------
subPrimaryAnnotation <- "Oligos"

resultsList <- resultsLists[[subPrimaryAnnotation]]
so <- resultsList[[1]]

annotations_3rd <- list(
  "1" = "Oligo_3",
  "2" = "Oligo_2",
  "3" = "Opc_1",
  "4" = "Oligo_1", 
  "5" = "Oligo_4",
  "6" = "Oligo_1",
  "7" = "OligoProg_1"
)

saveAnalysis(subPrimaryAnnotation)
# findDEMarkers(subPrimaryAnnotation)

#' 
#' ### 06: Neurons
## -----------------------------------------------------------------------------
subPrimaryAnnotation <- "Neurons"

resultsList <- resultsLists[[subPrimaryAnnotation]]
so <- resultsList[[1]]

annotations_3rd <- list(
  "1" = "EX_9", 
  "2" = "INH_7",
  "3" = "INH_6",
  "4" = "INH_4",
  "5" = "INH_1",
  "6" = "EX_4",
  "7" = "EX_8",
  "8" = "EX_12",
  "9" = "INH_5",
  "10" = "INH_3", 
  "11" = "EX_5",
  "12" = "MN_1",
  "13" = "EX_2",
  "14" = "EX_10",
  "15" = "INH_9",
  "16" = "EX_11",
  "17" = "EX_3",
  "18" = "EX_6",
  "19" = "EX_1",
  "20" = "EX_7",
  "21" = "INH_2"
)

saveAnalysis(subPrimaryAnnotation)
# findDEMarkers(subPrimaryAnnotation)


#' 
#' ## 03: Add Third annotation to all cell clusters
#' 
## -----------------------------------------------------------------------------
subsets <- read.table(file = file.path(rdataDir, paste0(tissue, "_subsets.txt")), header = FALSE, stringsAsFactors = FALSE)$V1

mList <- list()
for (i in seq_along(subsets)){
  subPrimaryAnnotation <- subsets[[i]]

  m <- read.table(
    file = file.path(outDir, paste0(subPrimaryAnnotation, "_withoutRN01_ThirdAnnotation.csv")),
    header = TRUE,
    stringsAsFactors = FALSE
    )
  mList[[i]] <- m

}

metadata <- do.call(rbind, mList)

so.m <- readRDS(file.path(rdataDir, paste0(tissue, "_som_HVG3000_PrimaryAnnotationSymphony.rds")))
so.m <- AddMetaData(so.m, metadata)
unique(so.m$ThirdAnnotation)

so.m <- subset(so.m, subset = ThirdAnnotation != "NA")

saveRDS(so.m, file = file.path(rdataDir, paste0(tissue, "_som_HVG3000_PrimaryAnnotationSymphony_ThirdAnnotation.rds")))


#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../06_ThirdAnnotation.Rmd"), documentation = 2)

