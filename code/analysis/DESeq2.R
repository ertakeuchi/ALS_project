#' ---
#' title: "01 - Case vs Control DEG Analysis using DESeq2"
#' author: "ertakeuchi"
#' date: "24 August 2023"
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
library("Signac")
library("tidyverse")
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
library("DEGreport")
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
# directory <- "230831_experiment3/out_3"

outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)
resolutions <- c(0.2, 0.4, 0.6, 0.8)

#' 
#' 
#' ################################################################################
#' # Analysis 01: Pseudobulk DE Analysis using DESeq2
#' ################################################################################
#' 
## -----------------------------------------------------------------------------
# 00: Set environment
# |
# |-- 01: Data Preprocessing 
# |
# |-- 02: Wald test in each cluster
# |     |
# |     |-- 01: set the function
# |     |-- 02: run it on all clusters comparing case to control
# |
# |-- 03: Liklihood ratio test; LRT Function settings
# |     |
# |     |-- 01: set the function
# |     |-- 02: run it on all clusters comparing case to control

#' 
#' ## 01. Data Preprocessing
#' 
## -----------------------------------------------------------------------------
so.m <- readRDS(
    file = file.path(
        rdataDir, paste0(tissue, "_som_HVG3000_PrimaryAnnotationSymphony_ThirdAnnotation.rds")
        )
)
# head(so.m@meta.data)
so.m$group_id <- factor(so.m$Condition)
so.m$Onset <- factor(so.m$Onset)

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(object = so.m, slot = "counts", assay = "RNA")
metadata <- so.m@meta.data %>%
            dplyr::select(c("Age", "Sex", "group_id",
                "PrimaryAnnotation", "ThirdAnnotation", "SubjectID"))

# head(metadata)

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(so.m$ThirdAnnotation)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts),
                        colData = metadata)

## Preparing for Pseudobulk analysis

# kids is cluster_ids
kids <- levels(as.factor(colData(sce)$cluster_id))
# kids

# Total number of clusters
nk <- length(kids)
nk

sids <- levels(as.factor(colData(sce)$SubjectID))
# Total number of samples 
ns <- length(sids)
ns


## Determine the number of cells per sample
table(sce$SubjectID)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$SubjectID))
# n_cells
## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$SubjectID)
## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                n_cells, row.names = NULL) %>% 
                select(-"cluster_id")
ei

# Aggregating counts to the sample level for each cluster
groups <- colData(sce)[, c("cluster_id", "SubjectID")]
# groups

# Aggregate raw counts across cluster-sample groups
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

# Turn into a list and split the list into components for each cluster and transform
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

# # Visualize the number of cells per cluster-sample group
# p <- ggplot(d, aes(x = SubjectID, y = cluster_id, fill = n_cells)) +
#         geom_tile() +
#         geom_text(aes(label = n_cells), color = "black", size = 2) +
#         scale_fill_gradient(low = "white", high = "green") +
#         labs(title = "Number of cells by Cell Type", 
#                 x = "gene", y = "SubjectID") +
#         coord_equal() +
#         theme(axis.text.y=element_text(hjust=0), axis.title.y=element_blank(), 
#                 axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), 
#                 plot.title=element_text(hjust=0.5), legend.position="none") 
# p
# ggsave(
#         file.path(outDir, "Number_of_cells_by_Cell_Type.pdf"),
#         plot = p, height = 10, width = 8
#         )

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

# Create a data frame with the sample IDs, cluster IDs and group_id

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    SubjectID = de_samples)

gg_df <- left_join(gg_df, ei[, c("SubjectID", "group_id", "Age", "Sex")]) 
head(gg_df)

metadata <- gg_df %>%
        dplyr::select(cluster_id, SubjectID, group_id, Age, Sex) 
        
metadata$cluster_id <- factor(metadata$cluster_id)
metadata$group <- paste0(metadata$cluster_id, "_", metadata$group_id) %>%
    factor()
head(metadata)

# Generate vector of cluster IDs
clusters <- levels(metadata$cluster_id)
# clusters

#' 
#' 
#' ## 02. Wald test in each cluster
#' 
#' ### 01. set the function
#' 
## -----------------------------------------------------------------------------
# Function to run DESeq2 and get results for all clusters
## x is index of cluster in clusters vector on which to run function
## A is the sample group to compare
## B is the sample group to compare against (base level)
# A = 2
# B = 1
# ThirdAnnotation <- "Astro_1"
# index <- which(clusters == ThirdAnnotation)
# x <- index

get_dds_resultsAvsB <- function(x, A, B){
  
        ThirdAnnotation <- clusters[x]
        cluster_metadata <- metadata[which(metadata$cluster_id == ThirdAnnotation), ]
        # head(cluster_metadata)
        rownames(cluster_metadata) <- cluster_metadata$SubjectID
        # head(cluster_metadata)

        cluster_metadata$group_id <- factor(cluster_metadata$group_id, levels = c("HC", "ALS"))

        # subset specific cluster counts
        counts <- pb[[ThirdAnnotation]]

        cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
        # check metadata rownames and counts colnames are the same
        all(rownames(cluster_metadata) == colnames(cluster_counts))  


        # Create DESeq2 object
        dds <- DESeqDataSetFromMatrix(cluster_counts,
                                        colData = cluster_metadata,
                                        design = ~ Age + Sex + group_id
                                        )


        # quality check
        # データ可視化のためにカウントを変換
        rld <-  rlog(dds, blind=TRUE)

        # PCA をプロット
        pdf(file.path(outDir, paste0(ThirdAnnotation, "_DESeq2_PCA.png")))
        DESeq2::plotPCA(rld, intgroup = "group_id")
        dev.off()

        # オブジェクトから rlog 行列を抽出し、一対の相関値を計算
        rld_mat <- assay(rld)
        rld_cor <- cor(rld_mat)

        # ヒートマップをプロット
        breaksList = seq(0.8, 1, by = 0.01)
        p <- pheatmap(
          rld_cor,
          annotation  = cluster_metadata[, c("group_id"), drop=F],
          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
          breaks = breaksList,
          )
        ggsave(file.path(outDir, paste0(ThirdAnnotation, "_DESeq2_heatmap.pdf")), plot = p, height = 10, width = 10)

        dds <- DESeq(dds)

        # 分散推定値をプロットする
        pdf(file.path(outDir, paste0(ThirdAnnotation, "_DESeq2_dispersion.png")))
        plotDispEsts(dds)
        dev.off()

        contrast <-  c("group_id", levels(cluster_metadata$group_id)[2], levels(cluster_metadata$group_id)[1])
        
        # check the coefficient names
        # resultsNames(dds) 
        res <- results(
          dds, 
          contrast  = contrast
        #   alpha =  0.05
          )
        res <-  lfcShrink(dds,
                        coef = 4,
                        res = res
                        )

        # 01.全遺伝子の結果表
        # Tidyverse関数で使用するために結果オブジェクトをtibbleにする
        res_tbl <- res  %>%
                data.frame()  %>%
                rownames_to_column(var="gene")  %>%
                as_tibble()

        # 結果出力をチェック
        res_tbl

        # 全ての結果をファイルに書き込む
        write.csv(res_tbl,
                file.path(
                        outDir,
                        paste0(
                                ThirdAnnotation, "_",
                                levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1],
                                "_all_genes.csv"
                                )
                        ),
                quote  =  FALSE, 
                row.names =  FALSE)


        # 02. 有意な遺伝子の表 
        # Set thresholds
        # padj_cutoff <- 0.05

        sig_res <-  dplyr::filter(res_tbl, padj < padj_cutoff)  %>%
                dplyr::arrange(padj)

        if(nrow(sig_res)>=1){

                write.csv(sig_res,
                        file.path(
                                outDir,
                                paste0(
                                        ThirdAnnotation, "_", 
                                        levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1],
                                        "_padj_", padj_cutoff, "_sig_genes.csv"
                                        )
                                ),
                        quote  =  FALSE, 
                        row.names =  FALSE)


                # 03. もっとも重要な上位20遺伝子の正規化発現の散布図
                ## ggplot of top genes
                normalized_counts <- counts(dds, normalized = TRUE)

                ## Order results by padj values
                top20_sig_genes <- sig_res  %>%
                        dplyr::arrange(padj)  %>%
                        dplyr::pull(gene) %>%
                        head(n=20)


                top20_sig_norm <- data.frame(normalized_counts)  %>%
                        rownames_to_column(var  = "gene")  %>%
                        dplyr::filter(gene  %in% top20_sig_genes)

                gathered_top20_sig <- top20_sig_norm  %>%
                        gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key  = "samplename", value  = "normalized_counts")
                        
                gathered_top20_sig <- inner_join(ei[,  c("SubjectID", "group_id")], gathered_top20_sig, by  =  c("SubjectID"  = "samplename"))

                ## plot using ggplot2
                p <- ggplot(gathered_top20_sig)  +
                        geom_point(aes(x = gene,
                                y  = normalized_counts, 
                                color  = group_id), 
                                position=position_jitter(w=0.1,h=0))  +
                        scale_y_log10()  +
                        xlab("Genes")  +
                        ylab("log10 Normalised Counts")  +
                        ggtitle("Top 20 Significant DE Genes")  +
                        theme_bw()  +
                        theme(axis.text.x  = element_text(angle  = 45,  hjust  =  1))  + 
                        theme(plot.title  = element_text(hjust =  0.5))

                # p

                ggsave(
                file.path(outDir, paste0(ThirdAnnotation, "_padj_", padj_cutoff, "_DESeq2_top20_sig_genes.pdf")),
                plot = p,
                height = 10, width = 10
                )

                # 04. 全有意遺伝子のヒートマップ
                sig_norm <- data.frame(normalized_counts) %>%
                        rownames_to_column(var = "gene") %>%
                        dplyr::filter(gene %in% sig_res$gene)
                
                if(nrow(sig_norm)>=2){
                        # カラーパレットの設定
                        heat_colors <- brewer.pal(6, "YlOrRd")

                        # アノテーションのメタデータデータフレームを用いて pheatmap を実行
                        p <- pheatmap(sig_norm[ ,  2:length(colnames(sig_norm))], 
                        colors  = heat_colors, 
                        cluster_rows  =  T, 
                        show_rownames  =  F,
                        annotation  = cluster_metadata[, c("group_id", "cluster_id")], 
                        border_color  = NA, 
                        fontsize  = 10, 
                        scale  = "row", 
                        fontsize_row  = , 
                        height = 20)        

                        # p
                        ggsave(file.path(outDir, paste0(ThirdAnnotation, "_padj_", padj_cutoff, "_DESeq2_Sigheatmap.pdf")), plot = p, height = 10, width = 10)
                }

                # 05. 結果のボルケーノプロット
                ## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
                res_table_thres <- res_tbl %>% 
                                mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.25)
                                
                ## Volcano plot
                p <- ggplot(res_table_thres) +
                geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
                ggtitle("Volcano plot of stimulated B cells relative to control") +
                xlab("log2 fold change") + 
                ylab("-log10 adjusted p-value") +
                scale_y_continuous(limits = c(0,50)) +
                theme(legend.position = "none",
                        plot.title = element_text(size = rel(1.5), hjust = 0.5),
                        axis.title = element_text(size = rel(1.25))) 
                
                # p
                ggsave(file.path(outDir, paste0(ThirdAnnotation, "_padj_", padj_cutoff, "_DESeq2_volcano.pdf")), plot = p, height = 10, width = 10)

        }else{print("No significant genes")}
        
}

#' 
#' 
#' ### 02. run it on all clusters comparing case to control
## -----------------------------------------------------------------------------


padj_list <- c(0.2, 0.1, 0.05, 0.01)
# padj_list <- c(0.2)

for (i in seq_along(padj_list)){
        padj_cutoff <- padj_list[i]
        # Run the script on all clusters comparing case to control 
        map(1:length(clusters), get_dds_resultsAvsB, A = 2, B = 1)
}

#' 
#' 
#' 
#' 
# ' ## 04. Liklihood ratio test
# ' 
# ' ### 01. set the function
# ' 
# -----------------------------------------------------------------------------

ThirdAnnotation <- "Opc_1"
index <- which(clusters == ThirdAnnotation)
x <- index
padj_cutoff <- 0.01

get_dds_LRTresults <- function(x){

        ThirdAnnotation <- clusters[x]
        cluster_metadata <- metadata[which(metadata$cluster_id == ThirdAnnotation), ]
        rownames(cluster_metadata) <- cluster_metadata$SubjectID
        cluster_metadata$group_id <- factor(cluster_metadata$group_id, levels = c("HC", "ALS"))
        counts <- pb[[ThirdAnnotation]]
        cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
        
        #all(rownames(cluster_metadata) == colnames(cluster_counts))        

        dds <- DESeqDataSetFromMatrix(cluster_counts,
                                colData = cluster_metadata,
                                design = ~ Age + Sex + group_id
                                )

        dds_lrt <- DESeq2::DESeq(
            object = dds,
            test = "LRT",
        #     fitType = "glmGamPoi",
        #     useT = TRUE,
        #     minmu = 1e-6,
        #     minReplicatesForReplace = Inf,
            reduced = ~ Age + Sex
            )

        # Extract results from a DESeq2 analysis
        res_LRT <- DESeq2::results(dds_lrt)

        # Create a tibble for LRT results
        res_LRT_tb <- res_LRT %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>% 
                as_tibble()
        head(res_LRT_tb)

        # Save all results
        write.csv(res_LRT_tb,
                    file = file.path(
                    outDir, paste0(ThirdAnnotation, "_LRT_AllGenes.csv")),
                    quote = FALSE,
                    row.names = FALSE) 

        # Subset to return genes with padj cutoff
        # padj_cutoff <- 0.01
        sigLRT_genes <- res_LRT_tb %>% 
        filter(padj < padj_cutoff)

        sigLRT_genes
        nrow(sigLRT_genes)

        if(nrow(sigLRT_genes) > 15){
                
                # Save all results
                write.csv(sigLRT_genes,
                        file = file.path(
                        outDir, paste0(ThirdAnnotation, "_padj_", padj_cutoff, "_LRT_SigGenes.csv")),
                        quote = FALSE,
                        row.names = FALSE) 

                
                # Transform counts for data visualization
                rld <- rlog(dds_lrt, blind=TRUE)
                
                # Extract the rlog matrix from the object and compute pairwise correlation values
                rld_mat <- assay(rld)
                rld_cor <- cor(rld_mat)
                
                
                # Obtain rlog values for those significant genes
                cluster_rlog <- rld_mat[sigLRT_genes$gene, ]

                cluster_meta_sig <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
                
                result <- tryCatch({
                        # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
                        cluster_groups <- degPatterns(cluster_rlog, metadata = cluster_meta_sig, time = "group_id", col=NULL)
                        ggsave(file.path(outDir, paste0(ThirdAnnotation, "_padj_", padj_cutoff, "_LRT_DEgene_groups.png")))

                        # Let's see what is stored in the `df` component
                        write.csv(cluster_groups$df,
                                file.path(outDir, paste0(ThirdAnnotation, "_padj_", padj_cutoff, "_LRT_DEgene_groups.csv")),
                                quote = FALSE, 
                                row.names = FALSE)

                }, error = function(e) {
                        print("Error in degPatterns function")
                })

                if(inherits(result, "try-error")){
                        next
                }


                
                # saveRDS(cluster_groups, file.path(outDir, paste0(ThirdAnnotation, "_LRT_DEgene_groups.rds")))
                # save(dds_lrt, cluster_groups, res_LRT, sigLRT_genes, file = file.path(paste0(ThirdAnnotation, "_all_LRTresults.Rdata")))
        }
}

#' 
#' ### 02. run it on all clusters
## -----------------------------------------------------------------------------

padj_list <- c(0.05, 0.01, 0.005)

for (i in seq_along(padj_list)){
        padj_cutoff <- padj_list[i]
        map(1:length(clusters), get_dds_LRTresults)
}

# ' 
# ' 
# ' ## 1.1 VennDiagram MAST DEGs and DESeq2 DEGs
# -----------------------------------------------------------------------------
FC <- 1.2 # Fold-change cutoff
FDR <- 0.2 # FDR cutoff
alpha <- FDR 
tissue <- "spinalcord"


# upかdownか
for (i in seq_along(c("up", "down"))){
    # i = 1
    up_down <- c("up", "down")[i]
    # DEGsの読込
    df <- read.csv(file = file.path(
        outDir, "table", tissue,
        paste0(tissue, "_DEGs_MAST_LFC025_Padj005_withoutRN01_", up_down, ".csv")))

    celltype_list <- unique(df$detailed_annotation)

    # celltypeのどれか 
    for(k in seq_along(celltype_list)){
        # k = 24
        clustx <- celltype_list[[k]]
        df_sub <- subset(df, detailed_annotation == clustx)
        genes <- df_sub$gene

        # significant genes from DESeq2 pseudo-bulk analysis    
        if(up_down == "up"){

            deseq2_genes_up <- read.csv(file = file.path(
                            outDir, "table", tissue, "DESeq2/subset",
                            paste0(clustx, "_res_FDR_", FDR, "_FC_", FC, "_withoutRN01_up.csv"))) %>%
                            dplyr::pull(X)

            if (length(deseq2_genes_up) >= 1) {

                overlap_genes <- intersect(genes, deseq2_genes_up)
                num_overlap <- length(overlap_genes)

                # ベン図を作成する
                p.venn <- VennDiagram::venn.diagram(
                    x = list(genes, deseq2_genes_up),
                    category.names = c(paste0(clustx, "_", up_down, "_MAST"), "DESeq2"),
                    filename = file.path(outDir, "img/spinalcord/DESeq2/subset",
                        paste0(tissue, "_", clustx, "_", up_down, "_genes_vs_DESeq2_FC", FC, "_padj_", FDR,  "_withoutRN01_venn_diagramm.png")),
                    # height = 500,
                    # width = 500,
                    col = c(custom_colors$discrete[5],custom_colors$discrete[2]),
                    fill = c(custom_colors$discrete[5],custom_colors$discrete[2]),
                    fontfamily = "sans",
                    cat.fontfamily = "sans",
                    cat.default.pos = "outer",
                    cat.pos = c(-15, 15)
                )

                # ベン図に表示する遺伝子名の取得
                venn_genes <- overlap_genes
                write.table(venn_genes, file = file.path(outDir, "img/spinalcord/DESeq2/subset",
                    paste0(tissue, "_", clustx, "_", up_down, "_intersect_MAST_DESeq2_genes_withoutRN01.txt")),
                    sep = " ",
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
            }
        } else {

                deseq2_genes_down <- read.csv(file = file.path(
                            outDir, "table", tissue, "DESeq2/subset",
                            paste0(clustx, "_res_FDR_", FDR, "_FC_", FC, "_withoutRN01_down.csv"))) %>%
                            dplyr::pull(X)
                
                if (length(deseq2_genes_down) >= 1) {

                    overlap_genes <- intersect(genes, deseq2_genes_down)
                    num_overlap <- length(overlap_genes)

                    p.venn <- VennDiagram::venn.diagram(
                            x = list(genes, deseq2_genes_down),
                            category.names = c(paste0(clustx, "_", up_down, "_MAST"), "DESeq2"),
                            filename = file.path(outDir, "img/spinalcord/DESeq2/subset",
                                paste0(tissue, "_", clustx, "_", up_down, "_genes_vs_DESeq2_FC", FC, "_padj_", FDR,  "_withoutRN01_venn_diagramm.png")),
                            # height = 500,
                            # width = 500,
                            col = c(custom_colors$discrete[13],custom_colors$discrete[2]),
                            fill = c(custom_colors$discrete[13],custom_colors$discrete[2]),
                            fontfamily = "sans",
                            cat.fontfamily = "sans",
                            cat.default.pos = "outer",
                            cat.pos = c(-15, 15)
                    )

                    # ベン図に表示する遺伝子名の取得
                    venn_genes <- overlap_genes
                    write.table(venn_genes, file = file.path(outDir, "img/spinalcord/DESeq2/subset",
                        paste0(tissue, "_", clustx, "_", up_down, "_intersect_genes_withoutRN01.txt")),
                        sep = " ",
                        quote = FALSE, row.names = FALSE, col.names = FALSE)
                    }
            }
        
    }
}


# ' 
# ' 
# ' ## SessionInfo
# -----------------------------------------------------------------------------
sessionInfo()

#' 
#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../DESeq2_LRT.Rmd"), documentation = 2)

