#' ---
#' title: "01 - Case vs Control DEG Analysis using DESeq2"
#' author: "ertakeuchi"
#' date: "22 August 2023"
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
library("Signac")
library("tidyverse")
# use_python("/usr/bin/python3")
library("patchwork")
library("gridExtra")
library("ComplexHeatmap")
library("pheatmap")
library("RColorBrewer")
library("data.table")
library("scater")
library("SingleCellExperiment")
library("Matrix.utils")
library("magrittr")
library("purrr")
library("reshape2")
library("S4Vectors")
library("apeglm")
library("DESeq2")
library("RColorBrewer")
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

tissue <- "spinalcord"
directory <- "230822_experiment1/out_1"

outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)
resolutions <- c(0.2, 0.4, 0.6, 0.8)

#' 
#' 
#' ################################################################################
#' # Analysis 01: Pseudobulk DE Analysis using DESeq2
#' ################################################################################
#' 
#' ## 1. Data Preprocessing
#' 
## -----------------------------------------------------------------------------
so.m <- readRDS(
    file = file.path(
        rdataDir, paste0(tissue, "_som_HVG3000_PrimaryAnnotationSymphony_ThirdAnnotation.rds")
        )
)
# clusters <- unique(so.m$ThirdAnnotation)

head(so.m@meta.data)
so.m$Sex <- factor(so.m$Sex)
so.m$Condition <- factor(so.m$Condition)
so.m$Onset <- factor(so.m$Onset)

# so.m@meta.data$detailed_annotation <- so.m@meta.data$ThirdAnnotation

# Extract Project RN01
# so.m <- subset(so.m, subset = Project == "RN01", invert = TRUE)
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(object = so.m, slot = "counts", assay = "RNA")
metadata <- so.m@meta.data %>%
            dplyr::select(c("Age", "Sex", "Condition",
                "PrimaryAnnotation", "ThirdAnnotation", "SubjectID"))
head(metadata)

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(so.m$ThirdAnnotation)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts),
                        colData = metadata)

## Preparing for Pseudobulk analysis
# cluster_idのvector
kids <- levels(as.factor(colData(sce)$cluster_id))
kids
# write.table(kids, file = file.path(outDir, "ThirdAnnotation.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# Total number of clusters
nk <- length(kids)
nk

sids <- levels(as.factor(colData(sce)$SubjectID))
# Total number of samples 
ns <- length(sids)
ns

## Determine the number of cells per sample
# table(sce$SampleID)
table(sce$SubjectID)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$SubjectID))
## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$SubjectID)
## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                n_cells, row.names = NULL) %>% 
                select(-"cluster_id")
ei
# # Perform QC if not already performed
# dim(sce)

# # Calculate quality control (QC) metrics
# sce <- calculateQCMetrics(sce)

# # Get cells w/ few/many detected genes
# sce$is_outlier <- isOutlier(
#         metric = sce$total_features_by_counts,
#         nmads = 2, type = "both", log = TRUE)

# # Remove outlier cells
# sce <- sce[, !sce$is_outlier]
# dim(sce)

# ## Remove lowly expressed genes which have less than 10 cells with any counts
# sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

# dim(sce)
# Aggregating counts to the sample level for each cluster
groups <- colData(sce)[, c("cluster_id", "SubjectID")]
# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(
    t(counts(sce)), 
    groupings = groups,
    fun = "sum") 

class(pb)
dim(pb)
pb[1:6,1:6]
# aggr_counts[1:6, 1:6]
# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_", n = 3), function(x) paste(x[1], x[2], sep = "_"))
splitf

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                    factor(splitf)) %>%
        lapply(function(u) 
                # set_colnames(t(u), gsub("^(?:[^_]*_){3}([^_]+_[^_]+).*", "\\1", rownames(u))
                set_colnames(t(u), sapply(stringr::str_split(rownames(u), pattern = "_", n = 4), function(x) paste(x[3], x[4], sep = "_")
                ))
        ) 

class(pb)

# Explore the different components of list
str(pb)

# Print out the table of cells in each cluster-sample group
options(width = 100)
# cluster_id = detailed_annotation
table(sce$cluster_id, sce$SubjectID)

# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
        pb[[x]] %>%
                colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
        unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
        rep(names(pb)[x], 
            each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
        unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    SubjectID = de_samples)

gg_df <- left_join(gg_df, ei[, c("SubjectID", "Condition", "Age", "Sex")]) 
head(gg_df)

metadata <- gg_df %>%
        dplyr::select(cluster_id, SubjectID, Condition, Age, Sex) 
        
metadata$cluster_id <- factor(metadata$cluster_id)
metadata$group <- paste0(metadata$cluster_id, "_", metadata$Condition) %>%
    factor()
head(metadata)

# Generate vector of cluster IDs
clusters <- levels(metadata$cluster_id)
# clusters

#' 
#' 
#' ## 2. Liklihood ratio test; LRT Function settings
#' 
## -----------------------------------------------------------------------------

get_dds_LRTresults <- function(x){
        # x = 1
        # padj.cutoff <- 0.2
        # fc <- 2 ^ 0.25

        cluster_metadata <- metadata[which(metadata$cluster_id == clusters[x]), ]
        rownames(cluster_metadata) <- cluster_metadata$SubjectID
        cluster_metadata$Condition <- factor(cluster_metadata$Condition, levels = c("HC", "ALS"))
        counts <- pb[[clusters[x]]]
        cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
        
        #all(rownames(cluster_metadata) == colnames(cluster_counts))        

        dds <- DESeqDataSetFromMatrix(cluster_counts,
                                colData = cluster_metadata,
                                design = ~ Age + Sex + Condition)

        dds_lrt <- DESeq2::DESeq(
            object = dds,
            test = "LRT",
            fitType = "glmGamPoi",
            useT = TRUE,
            minmu = 1e-6,
            minReplicatesForReplace = Inf,
            reduced = ~ Age + Sex)

        res_LRT <- DESeq2::results(dds_lrt)

        # Create a tibble for LRT results
        res_LRT_tb <- res_LRT %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>% 
        as_tibble()

        # Save all results
        write.csv(res_LRT_tb,
                    file = file.path(
                    outDir, paste0(clusters[x], "_LRT_AllGenes.csv")),
                    quote = FALSE,
                    row.names = FALSE) 

        # Subset to return genes with padj < 0.05
        sigLRT_genes <- res_LRT_tb %>% 
            filter(padj < padj.cutoff)
        sigLRT_genes

        # Save all results
        write.csv(sigLRT_genes,
                    file = file.path(
                    outDir, paste0(clusters[x], "_padj_", padj.cutoff, "_LRT_SigGenes.csv")),
                    quote = FALSE,
                    row.names = FALSE) 

        # Get number of significant genes
        nrow(sigLRT_genes)

        res_LRT_filtered <- subset(res_LRT_tb, padj < padj.cutoff & abs(log2FoldChange) > log2(fc)) # Select
        # table(sign(res_LRT_tb$log2FoldChange)) # N. of genes Down, Up
        res_LRT_filtered <- res_LRT_filtered[order(-res_LRT_filtered$log2FoldChange), ] #sort
        head(res_LRT_filtered) #top upregulated
        tail(res_LRT_filtered) #top downregulated
        res_LRT_filtered_up <- as.data.frame(res_LRT_filtered[res_LRT_filtered$log2FoldChange > 0, ])
        head(res_LRT_filtered_up)
        res_LRT_filtered_down <- as.data.frame(res_LRT_filtered[res_LRT_filtered$log2FoldChange < 0, ])
        head(res_LRT_filtered_down)

        write.csv(res_LRT_filtered_up,
                    file = file.path(
                    outDir, paste0(clusters[x], "_padj", padj.cutoff, "L2FC", log2(fc) ,"_LRT_SigGenesUp.csv")),
                    quote = FALSE,
                    row.names = FALSE)
        write.csv(res_LRT_filtered_down,
                    file = file.path(
                    outDir, paste0(clusters[x], "_padj", padj.cutoff, "L2FC", log2(fc) ,"_LRT_SigGenesDown.csv")),
                    quote = FALSE,
                    row.names = FALSE)
               
        # # Transform counts for data visualization
        # rld <- rlog(dds_lrt, blind=TRUE)
        
        # # Extract the rlog matrix from the object and compute pairwise correlation values
        # rld_mat <- assay(rld)
        # rld_cor <- cor(rld_mat)
        
        
        # # Obtain rlog values for those significant genes
        # cluster_rlog <- rld_mat[sigLRT_genes$gene, ]
        
        # cluster_meta_sig <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
        
        # # # Remove samples without replicates
        # # cluster_rlog <- cluster_rlog[, -1]
        # # cluster_metadata <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
        
        
        # # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
        # cluster_groups <- degPatterns(cluster_rlog, metadata = cluster_meta_sig, time = "group_id", col=NULL)
        # ggsave(paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.png"))
        
        # # Let's see what is stored in the `df` component
        # write.csv(cluster_groups$df,
        #         paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.csv"),
        #         quote = FALSE, 
        #         row.names = FALSE)
        
        # saveRDS(cluster_groups, paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.rds"))
        # save(dds_lrt, cluster_groups, res_LRT, sigLRT_genes, file = paste0("DESeq2/lrt/", clusters[x], "_all_LRTresults.Rdata"))
        
}

#' 
#' ## 3. Run LRT Function
#' 
## -----------------------------------------------------------------------------
padj_list <- c(0.2, 0.1, 0.05, 0.01)
fc_list <- c(2 ^ 0.25, 2 ^ 0.5, 2 ^ 1)

for (i in seq_along(padj_list)){
    # i = 1
    padj.cutoff <- padj_list[i]
    for (j in seq_along(fc_list)){
        # j = 1
        fc <- fc_list[j]
        result_list <- map(1:length(clusters), get_dds_LRTresults)
    }
}

#' 
## -----------------------------------------------------------------------------

# #-----
#         ## Extract normalized counts from dds object
#         normalized_counts <- counts(dds, normalized = TRUE)
        
#         ## Extract top 20 DEG from resLFC (make sure to order by padj)
#         top20_sig_genes <- res_up %>%
#             dplyr::arrange(padj) 
#         top20_sig_genes <- rownames(top20_sig_genes)
            
#         # top20_sig_genesがNULLでない場合のみ処理を実行
#         if (length(top20_sig_genes) > 0) {
#             ## Extract matching normalized count values from matrix
#             top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
#                 # top20_sig_genesが1つの場合のみ行名を修正
#                 if (length(top20_sig_genes) == 1) {
#                     top20_sig_counts <- data.frame(t(top20_sig_counts), row.names = top20_sig_genes)

#                 }
#             ## Convert wide matrix to long data frame for ggplot2
#             top20_sig_df <- data.frame(top20_sig_counts)
#             top20_sig_df$gene <- rownames(top20_sig_counts)
            
#             top20_sig_df <- melt(setDT(top20_sig_df), 
#                                 id.vars = c("gene"),
#                                 variable.name = "cluster_sample_id") %>% 
#                 data.frame()
            
#             ## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
#             top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
#             top20_sig_df$cluster_sample_id <- gsub("\\  ", "+ ", top20_sig_df$cluster_sample_id)
            
#             ## Join counts data frame with metadata
#             top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
#                                         by = "cluster_sample_id")
            
#             ## Generate plot
#             ggplot(top20_sig_df, aes(y = value, x = Condition, col = Condition, label = SubjectID)) +
#                 geom_jitter(height = 0, width = 0.15) +
#                 scale_y_continuous(trans = 'log10') +
#                 ylab("log10 of normalized expression level") +
#                 xlab("condition") +
#                 ggtitle(paste0(clustx, " Top 20 Significant DE Genes up without RN01")) +
#                 theme(plot.title = element_text(hjust = 0.5)) +
#                 facet_wrap(~ gene)
            
#             ggsave(file = file.path(
#                 outDir, "img", tissue, "DESeq2/subset/pairwise", paste0(
#                     tissue, "_", clustx, "_DESeq2_Top_20_Significant_DE_Genes_up_withoutRN01.pdf")), height = 20, width = 10)
        
#             } else {
#                 print(paste0(clustx, " has No significant genes."))
#             }

#         top20_sig_genes <- res_down %>%
#             dplyr::arrange(padj) 
#         top20_sig_genes <- rownames(top20_sig_genes)
            
#         # top20_sig_genesがNULLでない場合のみ処理を実行
#         if (length(top20_sig_genes) > 0) {
#             ## Extract matching normalized count values from matrix
#             top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
#                 # top20_sig_genesが1つの場合のみ行名を修正
#                 if (length(top20_sig_genes) == 1) {
#                     top20_sig_counts <- data.frame(t(top20_sig_counts), row.names = top20_sig_genes)

#                 }
#             ## Convert wide matrix to long data frame for ggplot2
#             top20_sig_df <- data.frame(top20_sig_counts)
#             top20_sig_df$gene <- rownames(top20_sig_counts)
            
#             top20_sig_df <- melt(setDT(top20_sig_df),
#                                 id.vars = c("gene"),
#                                 variable.name = "cluster_sample_id") %>%
#                 data.frame()
            
#             ## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
#             top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
#             top20_sig_df$cluster_sample_id <- gsub("\\  ", "+ ", top20_sig_df$cluster_sample_id)
            
#             ## Join counts data frame with metadata
#             top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
#                                         by = "cluster_sample_id")
            
#             ## Generate plot
#             ggplot(top20_sig_df, aes(y = value, x = Condition, col = Condition, label = SubjectID)) +
#                 geom_jitter(height = 0, width = 0.15) +
#                 scale_y_continuous(trans = 'log10') +
#                 ylab("log10 of normalized expression level") +
#                 xlab("condition") +
#                 ggtitle(paste0(clustx, " Top 20 Significant DE Genes down withoutRN01")) +
#                 theme(plot.title = element_text(hjust = 0.5)) +
#                 facet_wrap(~ gene)
            
#             ggsave(file = file.path(
#                 outDir, "img", tissue, "DESeq2/subset/pairwise", paste0(
#                     tissue, "_", clustx, "_DESeq2_Top_20_Significant_DE_Genes_down_withoutRN01.pdf")), height = 20, width = 10)
        
#             } else {
#                 print(paste0(clustx, " has No significant genes."))
#             }

#         # }

#         # map(cluster_names, get_dds_resultsAvsB, A = "ALS", B = "CTRL")
# # }


#' 
#' ## 1.1 VennDiagram MAST DEGs and DESeq2 DEGs
## -----------------------------------------------------------------------------
# FC <- 1.2 # Fold-change cutoff
# FDR <- 0.2 # FDR cutoff
# alpha <- FDR 
# tissue <- "spinalcord"


# # upかdownか
# for (i in seq_along(c("up", "down"))){
#     # i = 1
#     up_down <- c("up", "down")[i]
#     # DEGsの読込
#     df <- read.csv(file = file.path(
#         outDir, "table", tissue,
#         paste0(tissue, "_DEGs_MAST_LFC025_Padj005_withoutRN01_", up_down, ".csv")))

#     celltype_list <- unique(df$detailed_annotation)

#     # celltypeのどれか 
#     for(k in seq_along(celltype_list)){
#         # k = 24
#         clustx <- celltype_list[[k]]
#         df_sub <- subset(df, detailed_annotation == clustx)
#         genes <- df_sub$gene

#         # significant genes from DESeq2 pseudo-bulk analysis    
#         if(up_down == "up"){

#             deseq2_genes_up <- read.csv(file = file.path(
#                             outDir, "table", tissue, "DESeq2/subset",
#                             paste0(clustx, "_res_FDR_", FDR, "_FC_", FC, "_withoutRN01_up.csv"))) %>%
#                             dplyr::pull(X)

#             if (length(deseq2_genes_up) >= 1) {

#                 overlap_genes <- intersect(genes, deseq2_genes_up)
#                 num_overlap <- length(overlap_genes)

#                 # ベン図を作成する
#                 p.venn <- VennDiagram::venn.diagram(
#                     x = list(genes, deseq2_genes_up),
#                     category.names = c(paste0(clustx, "_", up_down, "_MAST"), "DESeq2"),
#                     filename = file.path(outDir, "img/spinalcord/DESeq2/subset",
#                         paste0(tissue, "_", clustx, "_", up_down, "_genes_vs_DESeq2_FC", FC, "_padj_", FDR,  "_withoutRN01_venn_diagramm.png")),
#                     # height = 500,
#                     # width = 500,
#                     col = c(custom_colors$discrete[5],custom_colors$discrete[2]),
#                     fill = c(custom_colors$discrete[5],custom_colors$discrete[2]),
#                     fontfamily = "sans",
#                     cat.fontfamily = "sans",
#                     cat.default.pos = "outer",
#                     cat.pos = c(-15, 15)
#                 )

#                 # ベン図に表示する遺伝子名の取得
#                 venn_genes <- overlap_genes
#                 write.table(venn_genes, file = file.path(outDir, "img/spinalcord/DESeq2/subset",
#                     paste0(tissue, "_", clustx, "_", up_down, "_intersect_MAST_DESeq2_genes_withoutRN01.txt")),
#                     sep = " ",
#                     quote = FALSE, row.names = FALSE, col.names = FALSE)
#             }
#         } else {

#                 deseq2_genes_down <- read.csv(file = file.path(
#                             outDir, "table", tissue, "DESeq2/subset",
#                             paste0(clustx, "_res_FDR_", FDR, "_FC_", FC, "_withoutRN01_down.csv"))) %>%
#                             dplyr::pull(X)
                
#                 if (length(deseq2_genes_down) >= 1) {

#                     overlap_genes <- intersect(genes, deseq2_genes_down)
#                     num_overlap <- length(overlap_genes)

#                     p.venn <- VennDiagram::venn.diagram(
#                             x = list(genes, deseq2_genes_down),
#                             category.names = c(paste0(clustx, "_", up_down, "_MAST"), "DESeq2"),
#                             filename = file.path(outDir, "img/spinalcord/DESeq2/subset",
#                                 paste0(tissue, "_", clustx, "_", up_down, "_genes_vs_DESeq2_FC", FC, "_padj_", FDR,  "_withoutRN01_venn_diagramm.png")),
#                             # height = 500,
#                             # width = 500,
#                             col = c(custom_colors$discrete[13],custom_colors$discrete[2]),
#                             fill = c(custom_colors$discrete[13],custom_colors$discrete[2]),
#                             fontfamily = "sans",
#                             cat.fontfamily = "sans",
#                             cat.default.pos = "outer",
#                             cat.pos = c(-15, 15)
#                     )

#                     # ベン図に表示する遺伝子名の取得
#                     venn_genes <- overlap_genes
#                     write.table(venn_genes, file = file.path(outDir, "img/spinalcord/DESeq2/subset",
#                         paste0(tissue, "_", clustx, "_", up_down, "_intersect_genes_withoutRN01.txt")),
#                         sep = " ",
#                         quote = FALSE, row.names = FALSE, col.names = FALSE)
#                     }
#             }
        
#     }
# }


#' 
#' 
#' ## SessionInfo
## -----------------------------------------------------------------------------
sessionInfo()

#' 
#' 
## -----------------------------------------------------------------------------
knitr::purl(file.path(outDir, "../DESeq2_LRT.Rmd"), documentation = 2)

