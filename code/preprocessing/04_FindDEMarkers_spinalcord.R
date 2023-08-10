#' ---
#' title: "04 - Find DE marker genes"
#' author: "ertakeuchi"
#' date: "08 August 2023"
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
library("SeuratDisk")
library("tidyverse")
library("reticulate")
use_python("/usr/bin/python3")
library("patchwork")
library("Signac")
library("progressr")
library("gridExtra")  
# library("leiden")
set.seed(1234)

#' 
#' 
#' ## 0-1. Set a common directory (FOR CONTAINER)
## -----------------------------------------------------------------------------
# homeDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq")
homeDir <- file.path(".")
# NGSDir <- file.path(homeDir, "../../NGS_original")
metadataDir <- file.path(homeDir, "metadata")

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("Usage: Rscript test.R tissue output_directory", call. = FALSE)
}
# # Extract the tissue argument
tissue <- args[1]
directory <- args[2]
# feature <- args[3]
# dim <- args[4]
# resolution <- args[5]

# tissue <- "spinalcord"
# directory <- "230808_experiment2/out_2"
outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)

resolutions <- c(0.4, 0.6, 0.8, 1.0)

#' 
#' ################################################################################
#' # Preprocessing 04: Find DE marker genes 
#' ################################################################################
#' 
#' ## 02. Find DE marker genes in each cluster
## -----------------------------------------------------------------------------
# 01: Set environment
# |
# | NOTE: AT THE FOLLOWING SECTION <02> USE SEURAT DEVELOP VERSION (CONTAINER)
# | 
# |-- 02: FindDEMarkers in each cluster 
# |   |
# |   |-- 01: Load Seurat object and Filtering meta.data
# |   |
# |   |-- 02: Preparing SCT FindMarkers and save Seurat object
# |   |
# |   |-- 03: FindAllMarkers assay = "SCT"
# |   |
# |   |-- 04: Visualize top 10 genes
# | 
# |-- 03: Visualize cannonical markers
# |   |
# |   |-- 01: Set canonical markers
# |   |
# |   |-- 02: Dotplot
# |   |
# |   |-- 03: Featureplot
# |   |
# |   |-- 04: Dimplot

#' 
#' ### 01: Load Seurat object and Filtering meta.data 
## -----------------------------------------------------------------------------
# set decieded parameters
feature <- 5000
dim <- 50
resolution <- 0.6

filtered_resolutions <- resolutions[resolutions != resolution]

# Read the harmonized Seurat object
so <- readRDS(file.path(rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, ".rds")))
m <- so@meta.data

# metadataを整理する
# "DF_class_"を含む列を削除
m_filtered <- m %>%
  select(
    -contains(c("DF.classifications", "pANN_", "MULTI_", "SubjectID.x", 
                "SubjectID.y", "Project.x", "Project.y", "_HTO",
                paste0("SCT_snn_res.", filtered_resolutions)))
  )
m_filtered$seurat_clusters <- as.factor(m_filtered[[paste0("SCT_snn_res.", resolution)]])
# head(m_filtered)

so@meta.data <- m_filtered
Idents(so) <- "seurat_clusters"


#' 
#' ### 02. Preparing SCT FindMarkers
## ---- results='hide', include=FALSE-------------------------------------------
# handlers(global = TRUE) 

so <- PrepSCTFindMarkers(so)

so_features <- read_csv(file.path(outDir, "so_features.csv"))
so_features <- so_features[-1]  # "X"の列を削除

colnames(so_features) <- "gene"
# head(so_features)
so_features <- so_features$gene

saveRDS(
  so,
  file.path(
    rdataDir, paste0(tissue, "_som_HVG", feature, "_pca", dim, "_SCT.rds")
    )
  )

#' 
#' ### 03. FindAllMarkers assay = "SCT"
## -----------------------------------------------------------------------------
so.markers <- FindAllMarkers(
  so,
  assay = "SCT",
  only.pos = TRUE,
  test.use= "wilcox", # default
  min.pct = 0.1, # default 
  logfc.threshold = 1, # default = 0.25,
  features = so_features
  )

write.table(
  so.markers,
  file = file.path(outDir, "so_AllMarkersSCT.tsv"),
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t",
  quote = FALSE
  )


#' 
#' ### 04. Visualize top 10 genes
## -----------------------------------------------------------------------------
top10 <- so.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
head(top10)

write.table(
  top10,
  file = file.path(outDir, "so_top10_AllMarkersSCT.tsv"),
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t",
  quote = FALSE
  )

so.subsampled <- so[, sample(colnames(so), size =2999, replace=F)]

p <- DoHeatmap(
  so.subsampled,
  features = top10$gene,
  raster = TRUE
  ) + 
  NoLegend()

ggsave(
  file.path(outDir, "Heatmap_top10.pdf"), 
  plot = p,
  height = 20, width = 40
  )


#' 
#' 
#' ## 03: Visualize cannonical markers
#' 
#' ### 01: Set canonical markers
## ---- fig.height=10, fig.width=8----------------------------------------------
# http://doi.org/10.1038/s41586-021-03465-8

canonical_markers <- list(
  'Neuron' = c("RBFOX3", "SNAP25"), 
  'Oligo' = c("MBP", "MAG", "MOG", "MOBP", "PLP1"),
  'Oligo_prog' = c("TNR", "UST"),
  'Astro' = c("GFAP", "AQP4", "SLC7A10", "SLC1A2", "CCDC85A", "CPAMD8", "ITPRID1", "SLC6A11", "GABBR2", "HAP1"),
  'Micro' = c("MYO1F", "CSF1R", "LY86", "APBB1IP", "PRDM11", "CX3CR1", "PPARG", "NFATC1", "TMEM163"),
  'Endo' = c("FLT1", "PDGFRB"), # removed "LY6C"
  'Meninges' = c("COL1A2", "DCN", "COL6A3", "ALDH1A2", "LAMC3"),
  'Opc' = c("PDGFRA", "C1QL1"),
  'Pericyte' = c("ABCC9", "NOTCH3"),
  'Ependymal' = c("CFAP43", "DNAH12"), 
  'Schwann' = c("MPZ", "PMP22"),
  'Lymphoctye' = c("PTPRC"),
  'MNs' = c('CHAT', 'SLC5A7', 'ACLY', 'PRPH', 'NEFH', 'NEFM', 'STMN2'),
  'Proliferating Microglia' = c('POLQ', 'TOP2A', 'MKI67', 'MRC1', 'LYVE1', 'MARCO', 'F13A1')
  
)

marker_genes <- list(
  'Neuron' = c("RBFOX3"),
  'GABAergic' = c("GAD1", "GAD2"), 
  'glutamatergic' = c("SLC17A7", "SATB2"), 
  'OPC' = c("PDGFRA"), 
  'astro' = c("AQP4"),
  'oligo' = c("PLP1", "MOBP"), 
  'perivasc_macro' = c("MRC1"), 
  'Tcells' = c("PTPRC"), 
  'vasc_smooth' = c("PDGFRB"), 
  'vasc_endo' = c("FLT1"), 
  'vasc_fibro' = c("DCN"), 
  'micro' = c("APBB1IP")
)

# https://doi.org/10.1038/s41467-021-25125-1 p6 right 
MN_marker_genes <- list(
  'MNa' = c("SPP1", "POLN"),
  'MNg' = c("ESRRG", "HTR1F"), 
  'PGC' = c("GFRA3", "NOS1", "FBN2")
  
)

marker_gene_lst <- c(unlist(canonical_markers, use.names = F))
# marker_gene_lst

#' 
#' ### 02: Dotplot
## ---- fig.height=10, fig.width=8----------------------------------------------
# Dotplot
p <- DotPlot(
  so,
  features = marker_gene_lst,
  group.by = "seurat_clusters",
  dot.scale = 8
  ) +
  theme(
    axis.text.y=element_text(hjust=0, size = 15),
    axis.title.y=element_blank(), 
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 15)
  )

ggsave(
  file.path(outDir, paste0("MarkerGenesDotPlot.pdf")),
  height = 10, width = 12)

#' 
#' ### 03: Featureplot
## ---- fig.height=10, fig.width=8----------------------------------------------
# Featureplot
p_sumList <- list()

for(i in seq_along(canonical_markers)){
  # i = 1
  genes <- canonical_markers[[i]]
  celltype <- names(canonical_markers)[[i]]
  
  p_list <- list()
  for (j in seq_along(genes)){
    # j = 1
    gene <- genes[[j]]
    p <- FeaturePlot(
      so,
      features = gene,
      reduction = "umap_HM_01",
      raster = TRUE
    ) +
      ggtitle(paste0(celltype, "_", gene))

    p_list[[j]] <- p
  }

  p_sumList[[i]] <- do.call(rbind, p_list) 


  # グラフのサイズを指定
  graphWidth <- 4
  graphHeight <- 4

  # ページごとのグラフ数を指定
  graphsPerPage <- 16

  # ページ数を計算
  n_pages <- ceiling(length(p_list) / graphsPerPage)

  # 1ページずつPDFに保存
  for (page in 1:n_pages) {
    # グリッドに配置するグラフを選択
    startIdx <- (page - 1) * graphsPerPage + 1
    endIdx <- min(startIdx + graphsPerPage - 1, length(p_list))
    graphList <- p_list[startIdx:endIdx]

    # リスト内のグラフを4×4のグリッドに配置
    n_graphs <- length(graphList)
    n_cols <- ifelse(n_graphs < 4, n_graphs, 4)
    n_rows <- ceiling(n_graphs / n_cols)
    graphGrid <- arrangeGrob(
      grobs = graphList,
      nrow = n_rows,
      ncol = n_cols,
      widths = rep(graphWidth, n_cols),
      heights = rep(graphHeight, n_rows))

    # PDFに保存
      ggsave(
      file.path(outDir, paste0(celltype, page, "_FeaturePlot.pdf")),
      plot = graphGrid,
      width = graphWidth * n_cols,
      height = graphHeight * n_rows)
  }

}

#' 
#' ### 04: Dimplot
## ---- fig.height=10, fig.width=8----------------------------------------------
# visualize as Dimplot
p <- DimPlot(
  so,
  reduction = "umap_HM_01",
  repel = TRUE,
  label = TRUE,
  group.by = paste0("SCT_snn_res.", resolution),
  raster = TRUE
 ) + 
  ggtitle('UMAP colored by primary annotations') 

ggsave(
  file.path(outDir, "PrimaryDimPlot.pdf"),
  plot = p,
  height = 10, width = 10)

#' 
#' 
#' # Memo: Count cell numbers in a specific cluster
## ---- fig.height=4, fig.width=8-----------------------------------------------
s <- subset(so, subset = seurat_clusters == 18)
# table(s@meta.data$Condition)
t <- table(s@meta.data$SubjectID) %>% as.data.frame()
colnames(t) <- c("SubjectID", "count")
t$SubjectID <- as.character(t$SubjectID)
t$Condition <- lapply(strsplit(t$SubjectID, "_"), "[", 1)
# head(t)
# Condition列をファクターに変換
t$Condition <- factor(t$Condition, levels = c("ALS", "HC"))

# プロットの描画
p <- ggplot(t, aes(x = SubjectID, y = count, fill = Condition)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = count),
    position = position_stack(vjust = 0.5),
    color = "white", size = 6) +
  coord_flip()
# p
ggsave(
  file.path(outDir, "MNsCellNumberSubjectID.pdf"),
  plot = p,
  height = 5,
  width = 10)

#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
## -----------------------------------------------------------------------------
knitr::purl(file.path(outDir, "../04_FindDEMarkers.Rmd"), documentation = 2)

