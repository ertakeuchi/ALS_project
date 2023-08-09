#' ---
#' title: "04 - Find DE marker genes"
#' author: "ertakeuchi"
#' date: "09 August 2023"
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
# Extract the tissue argument
tissue <- args[1]
directory <- args[2]
feature <- args[3]
dim <- args[4]
resolution <- args[5]

# tissue <- "brain"
# directory <- "230809_experiment2/out_2"
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
# # set decieded parameters
# feature <- 5000
# dim <- 50
# resolution <- 0.6

filtered_resolutions <- resolutions[resolutions != resolution]

t <- read.table(outDir, "../log.txt")

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

so_features <- read_csv(file.path(outDir, "../so_features.csv"))
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
  logfc.threshold = 0.5, # default = 0.25,
  features = so_features
  )

write.table(
  so.markers,
  file = file.path(outDir, "so_AllMarkersSCT_LFC05.tsv"),
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
  'Oligo' = c("MBP", "MOBP", "OPALIN", "PDGFRA"),
  'Astro' = c("GFAP", "AQP4", "ATP1A2", "SLC1A2"), 
  'Micro' = c("CX3CR1", "TREM2", "P2RY12", "FYB", "CSF1R", "CD74", "C3"),
  'Endo' = c("CLDN5", "SLC2A1", "FLT1"), 
  'Opc' = c("ADARB2", "PDGFRA", "VCAN", "CSPG4"), 
  'Pericyte' = c("ACTA2", "KCNJ8", "PDGFRB", "CD248", "ANPEP", "DES", "DLK1", "ZIC1", "ABCC9", "RGS5", "CSPG4"), 
  'VLMC' = c("COLEC12", "ITIH5", "COL1A2", "TBX18", "EBF1", "C7", "COL6A2", "SRPX2", "FLVCR2", "FMO2"), 
  'Endo' = c("EBF1", "ABCG2", "CLDN5", "FLI1", "LEF1", "EMCN", "IFI27", "HLA.E", "ADGRL4", "CLEC3B"), 
  'Lamp5' = c("FGF13", "PTPRT", "PRELID2", "GRIA4", "RELN", "PTCHD4", "EYA4", "MYO16", "FBXL7", "LAMP5"), 
  'Pvalb' = c("ADAMTS17", "ERBB4", "DPP10", "ZNF804A", "MYO16", "BTBD11", "GRIA4", "SLIT2", "SLIT2", "SDK1", "PVALB"), 
  'Sncg' = c("CNR1", "SLC8A1", "ASIC2", "CXCL14", "MAML3", "ADARB2", "NPAS3", "CNTN5", "FSTL5", "SNCG"), 
  'Sst' = c("GRIK1", "PALYL", "SST", "TRHDE", "GRID2", "NXPH1", "COL25A1", "SLC8A1", "SOX6", "ST6GALNAC5"), 
  'Sst_Chodl' = c("NPY", "FAM46A", "STAC", "OTOF", "NPY2R", "CRHBP", "ANKRD34B", "NOS1", "SST", "CHODL"), 
  'Vip' = c("GALNTL6", "LRP1B", "VIP", "GRM7", "KCNT2", "THSD7A", "ERBB4", "SYNPR", "ADARB2", "SLC24A3"), 
  'L23_IT' = c("CBLN2", "EPHA6", "LAMA2", "CNTN5", "PDZD2", "CUX2", "RASGRF2", "FAM19A1", "LINC01378", "CA10"), 
  'L5_ET' = c("COL5A2", "FAM19A1", "VAT1L", "COL24A1", "CBL2", "NRP1", "PTCHD1.AS", "NRG1", "HOMER1", "SLC35F3"), 
  'L5_IT' = c("FSTL4", "CNTN5", "RORB", "FSTL5", "IL1RAPL2", "CHN2", "TOX", "CPNE4", "CADPS2", "POU6F2"), 
  'L56_NP' = c("TSHZ2", "NPSR1.AS1", "HTR2C", "ITGA8", "ZNF385D", "ASIC2", "CDH6", "CRYM", "NXPH2", "CPNE4"), 
  'L6_CT' = c("ADAMTSL1", "KIAA1217", "SORCS1", "HS3ST4", "TRPM3", "TOX", "SEMA3E", "EGFEM1P", "MEIS2", "SEMA5A"), 
  'L6_IT' = c("PTPRL", "PDZRN4", "CDH9", "THEMIS", "FSTL5", "CDH13", "CDH12", "CBLN2", "LY86.AS1", "MLIP"), 
  'L6IT_Car' = c("THEMIS", "RNF152", "NTNG2", "STK32B", "KCNMB2", "GAS2L3", "OLFML2B", "POSTN", "B3GAT2", "NR4A2"), 
  'L6b' = c("HS3ST4", "KCNMB2", "MDFIC", "C10orf11", "NTM", "CDH9", "MARCH1", "TLE4", "FOXP2", "KIAA1217")
)


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
  'Micro' = c("APBB1IP", "CX3CR1", "P2RY12", "TREM2", "CSF1R", "CD74", "C3"),
)

neuronal_genes <- list(
  'Lamp5' = c("LAMP5"), 
  'Pvalb' = c("PVALB"), 
  'Sncg' = c("SNCG"), 
  'Sst' = c("SST"), 
  'Sst_Chodl' = c("CHODL"), 
  'Vip' = c("VIP"), 
  'L23_IT' = c("CUX2"), 
  'L5_ET' = c("COL5A2", "NRG1", "FEZF2"), 
  'L5_IT' = c("CPNE4"), 
  'L56_NP' = c("NXPH2"), 
  'L6_CT' = c("MEIS2"), 
  'L6_IT' = c("THEMIS"), 
  'L6IT_Car' = c("NTNG2"), 
  'L6b' = c("HS3ST4")
)

marker_gene_lst <- c(unlist(marker_genes, use.names = F))
neuron_gene_lst <- c(unlist(neuronal_genes, use.names = F))

#' 
#' ### 02: Dotplot
## ---- fig.height=10, fig.width=8----------------------------------------------
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

p2 <- DotPlot(
  so,
  features = neuron_gene_lst,
  group.by = "seurat_clusters",
  dot.scale = 10
  ) &
  labs("Neuronal markers") &
  theme(
    axis.text.y=element_text(hjust=0, size = 15),
    axis.title.y=element_blank(), 
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size = 15)
  )

# p <- p1 + p2
# p

ggsave(
  file.path(outDir, paste0("MarkerGenesDotPlot.pdf")),
  plot = p1,
  height = 10, width = 10)

ggsave(
  file.path(outDir, paste0("NeuronalGenesDotPlot.pdf")),
  plot = p2,
  height = 10, width = 10)

#' 
#' ### 03: Featureplot
## ---- fig.height=10, fig.width=8----------------------------------------------
# Featureplot
p_sumList <- list()

for(i in seq_along(marker_genes)){
  # i = 1
  genes <- marker_genes[[i]]
  celltype <- names(marker_genes)[[i]]
  
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
s <- subset(so, subset = seurat_clusters == 23)
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
  file.path(outDir, "BetsCellNumberSubjectID.pdf"),
  plot = p,
  height = 5,
  width = 10)

#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../04_FindDEMarkers.Rmd"), documentation = 2)

