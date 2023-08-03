#' ---
#' title: "DESeq2 pseudo-bulk Wald analysis"
#' author: "Eriko Takeuchi"
#' date: "20230801"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)

#' 
#' # Packages
## ---- include=FALSE-------------------------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
# library(DEGreport)
# BiocManager::install("DEGreport")
set.seed(1234)

#' 
#' # inputs
#' ## set directory
## -------------------------------------------------------------------------------------------------
mainDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq/R/out/final_out/spinalcord")
outDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq/experimental_record/230801_experiment1/out_1")
tissue <- "spinalcord"

#' ## Color setting
## -------------------------------------------------------------------------------------------------
custom_colors <- list()

colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)

colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

colors_indian <- c(
   '#FEA47F', '#25CCF7', '#EAB543', '#55E6C1', '#CAD3C8',
   '#F97F51', '#1B9CFC', '#F8EFBA', '#58B19F', '#2C3A47',
   '#B33771', '#3B3B98', '#FD7272', '#9AECDB', '#D6A2E8',
   '#6D214F', '#182C61', '#FC427B', '#BDC581', '#82589F'
)

custom_colors$discrete <- c(colors_dutch, colors_spanish, colors_indian)

#' 
#' 
#' ################################################################################
#' # Step 1: Pseudobulk DE Analysis using DESeq2
#' ################################################################################
#' 
#' ## 1. Data Preprocessing
#' 
## -------------------------------------------------------------------------------------------------
so.m <- readRDS(
    file = file.path(
        mainDir, "som_spinalcord_3rdannotation_LogNormalized.rds"
        )
)

# sexを数字に置き換えておく
so.m$Sex <- ifelse(so.m$Sex == "Male", 0, 1)
# Celltype_3rdを"EX_1"みたいな簡易版にかえる as DetailedAnnotation
so.m@meta.data$detailed_annotation <- sapply(stringr::str_split(so.m@meta.data$CellType_3rd, pattern = "_", n = 4), 
    function(x) paste(x[1], x[2], sep = "_")) 
head(so.m@meta.data)

# Extract Project RN01
so.m <- subset(so.m, subset = Project == "RN01", invert = TRUE)
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(object = so.m, slot = "counts", assay = "RNA")
metadata <- so.m@meta.data %>%
            dplyr::select(c("Age", "Sex", "Onset", "Condition",
                "CellType_1st", "detailed_annotation", "SampleID"))
head(metadata)

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(so.m$detailed_annotation)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts),
                        colData = metadata)

## Preparing for Pseudobulk analysis
# cluster_idのvector
kids <- levels(as.factor(colData(sce)$cluster_id))
# Total number of clusters
nk <- length(kids)
nk

sids <- levels(as.factor(colData(sce)$SampleID))
# Total number of samples 
ns <- length(sids)
ns

## Determine the number of cells per sample
table(sce$SampleID)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$SampleID))
## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$SampleID)
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
groups <- colData(sce)[, c("cluster_id", "SampleID")]
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
# splitf <- sapply(stringr::str_split(rownames(pb), 
#                                     pattern = "_",  
#                                     n = 2), 
#                 `[`, 2)
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
table(sce$cluster_id, sce$SampleID)

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
                    SampleID = de_samples)

gg_df <- left_join(gg_df, ei[, c("SampleID", "Condition", "Age", "Sex")]) 
head(gg_df)

metadata <- gg_df %>%
        dplyr::select(cluster_id, SampleID, Condition, Age, Sex) 
        
metadata$cluster_id <- factor(metadata$cluster_id)
metadata$group <- paste0(metadata$cluster_id, "_", metadata$Condition) %>%
    factor()
head(metadata)

# Generate vector of cluster IDs
clusters <- levels(metadata$cluster_id)
# clusters

#' 
#' 
#' ## 2. Wald test (for comparison)
## -------------------------------------------------------------------------------------------------
# Function to run DESeq2 and get results for all clusters
## x is index of cluster in clusters vector on which to run function
## A is the sample group to compare
## B is the sample group to compare against (base level)

get_dds_resultsAvsB <- function(x, A, B){
        # A = 2
        # B = 1
        # x = 1
        cluster_metadata <- metadata[which(metadata$cluster_id == clusters[x]), ] 
        rownames(cluster_metadata) <- cluster_metadata$sample_id
        counts <- pb[[clusters[x]]]
        cluster_counts <- data.frame(counts[, which(colnames(counts) %in% cluster_metadata$SampleID)])
        
        #all(rownames(cluster_metadata) == colnames(cluster_counts))        
        
        dds <- DESeqDataSetFromMatrix(cluster_counts,
                                colData = cluster_metadata,
                                design = ~ Age + Sex + Condition)
        
        # Transform counts for data visualization
        rld <- rlog(dds, blind=TRUE)
        
        # Plot PCA
        
        DESeq2::plotPCA(rld, intgroup = "Condition")
        ggsave(file.path(outDir, paste0(clusters[x], "_specific_PCAplot.png")))
        
        
        # Extract the rlog matrix from the object and compute pairwise correlation values
        rld_mat <- assay(rld)
        rld_cor <- cor(rld_mat)
        
        # Plot heatmap
        

        my_sample_col = cluster_metadata[,c("SampleID", "Condition", "Sex", "Age")]
        rownames(my_sample_col) <- my_sample_col$SampleID
        my_sample_col <- my_sample_col[,c("Condition", "Sex", "Age"), drop=F]
        head(my_sample_col)
        
        png(file.path(outDir, paste0(clusters[x], "_specific_heatmap.png")))

        pheatmap(
            rld_cor,
            annotation_col = my_sample_col)
        dev.off()
        
        # Run DESeq2 differential expression analysis
        # dds <- DESeq(dds)
        dds_wld <- DESeq2::DESeq(
            object = dds,
            test = "Wald",
            # fitType = "glmGamPoi", # glmGamPoi dispersion estimator should be used in combination with a LRT and not a Wald test.
            # useT = TRUE,
            minmu = 1e-6,
            minReplicatesForReplace = Inf)
        
        # Plot dispersion estimates
        png(file.path(outDir, paste0(clusters[x], "_dispersion_plot.png")))
        plotDispEsts(dds_wld)
        dev.off()

        # Output results of Wald test for contrast for A vs B
        cluster_metadata$Condition <- factor(cluster_metadata$Condition)
        contrast <- c("Condition", levels(cluster_metadata$Condition)[A], levels(cluster_metadata$Condition)[B])
        contrast        
        contrast_name <- paste0("Condition_", levels(cluster_metadata$Condition)[A], "_vs_", levels(cluster_metadata$Condition)[B])
        contrast_name
        # resultsNames(dds)
        res <- results(dds_wld, 
                       contrast = contrast,
                       alpha = 0.2,
                       pAdjustMethod = "BH" # default param
                       )
        
        res <- lfcShrink(dds_wld, 
                         coef = contrast_name,
                         res=res)
        # Set thresholds
        padj_cutoff <- 0.2
        
        # Turn the results object into a tibble for use with tidyverse functions
        res_tbl <- res %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()
        
        write.csv(res_tbl,
                  file.path(outDir, paste0(clusters[x], "_", contrast_name, "_WaldAllGenes.csv")),
                  quote = FALSE, 
                  row.names = FALSE)
                               
        # Subset the significant results
        sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
                dplyr::arrange(padj)
        sig_res

        sig_res_up <- as.data.frame(sig_res[sig_res$log2FoldChange > 0, ])
        head(sig_res_up)
        sig_res_down <- as.data.frame(sig_res[sig_res$log2FoldChange < 0, ])
        head(sig_res_down)

        write.csv(sig_res_up,
                  file.path(outDir, paste0(clusters[x], "_", contrast_name, "_WaldSigGenesUp.csv")),
                  quote = FALSE, 
                  row.names = FALSE)
        write.csv(sig_res_down,
                  file.path(outDir, paste0(clusters[x], "_", contrast_name, "_WaldSigGenesDown.csv")),
                  quote = FALSE, 
                  row.names = FALSE)
        
        ## ggplot of top genes
        normalized_counts <- counts(dds_wld, 
                                    normalized = TRUE)
        
        ## Order results by padj values
        top20_sig_genes <- sig_res %>%
                dplyr::arrange(padj) %>%
                dplyr::pull(gene) %>%
                head(n=20)
        
        
        top20_sig_norm <- data.frame(normalized_counts) %>%
                rownames_to_column(var = "gene") %>%
                dplyr::filter(gene %in% top20_sig_genes)
        
        gathered_top20_sig <- top20_sig_norm %>%
                gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
        
        gathered_top20_sig <- inner_join(ei[, c("SampleID", "Condition" )], gathered_top20_sig, by = c("SampleID" = "samplename"))
        
        ## plot using ggplot2
        ggplot(gathered_top20_sig) +
                geom_point(aes(x = gene, 
                               y = normalized_counts, 
                               color = Condition), 
                           position=position_jitter(w=0.1,h=0)) +
                scale_y_log10() +
                xlab("Genes") +
                ylab("log10 Normalized Counts") +
                ggtitle("Top 20 Significant DE Genes") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                theme(plot.title = element_text(hjust = 0.5))

        ggsave(file.path(outDir, paste0(clusters[x], "_", contrast_name, "_WaldTop20SigGenes.png")))
        
}

# Run the script on all clusters comparing stim condition relative to control condition
map(1:length(clusters), get_dds_resultsAvsB, A = 2, B = 1)



#' 
#' 
#' 
#' ## 2. Liklihood ratio test; LRT
#' 
## -------------------------------------------------------------------------------------------------

# get_dds_LRTresults <- function(x){
#         # x = 1
#         padj.cutoff <- 0.2
#         fc <- 2 ^ 0.25

#         cluster_metadata <- metadata[which(metadata$cluster_id == clusters[x]), ]
#         rownames(cluster_metadata) <- cluster_metadata$SampleID
#         cluster_metadata$Condition <- factor(cluster_metadata$Condition, levels = c("CTRL", "ALS"))
#         counts <- pb[[clusters[x]]]
#         cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
        
#         #all(rownames(cluster_metadata) == colnames(cluster_counts))        

#         dds <- DESeqDataSetFromMatrix(cluster_counts,
#                                 colData = cluster_metadata,
#                                 design = ~ Age + Sex + Condition)

#         dds_lrt <- DESeq2::DESeq(
#             object = dds,
#             test = "LRT",
#             fitType = "glmGamPoi",
#             useT = TRUE,
#             minmu = 1e-6,
#             minReplicatesForReplace = Inf,
#             reduced = ~ Age + Sex)

#         res_LRT <- DESeq2::results(dds_lrt)

#         # Create a tibble for LRT results
#         sig_res$<- res_LRT %>%
#         data.frame() %>%
#         rownames_to_column(var="gene") %>% 
#         as_tibble()

#         # Save all results
#         write.csv(sig_res$
#                     file = file.path(
#                     outDir, paste0(clusters[x], "_LRT_AllGenes.csv")),
#                     quote = FALSE,
#                     row.names = FALSE) 

#         # Subset to return genes with padj < 0.05
#         sigLRT_genes <- sig_res$%>% 
#             filter(padj < padj.cutoff)
#         sigLRT_genes

#         # Save all results
#         write.csv(sigLRT_genes,
#                     file = file.path(
#                     outDir, paste0(clusters[x], "_padj_", padj.cutoff, "_LRT_SigGenes.csv")),
#                     quote = FALSE,
#                     row.names = FALSE) 

#         # Get number of significant genes
#         nrow(sigLRT_genes)

#         sig_res$<- subset(sig_res$ padj < padj.cutoff & abs(log2FoldChange) > log2(fc)) # Select
#         # table(sign(sig_res$log2FoldChange)) # N. of genes Down, Up
#         sig_res$<- sig_res$order(-sig_res$log2FoldChange), ] #sort
#         head(sig_res$ #top upregulated
#         tail(sig_res$ #top downregulated
#         sig_res_up <- as.data.frame(sig_res$log2FoldChange > 0, ])
#         head(sig_res_up)
#         res_wld_tb_down <- as.data.frame(sig_res$log2FoldChange < 0, ])
#         head(res_wld_tb_down)

#         write.csv(sig_res_up,
#                     file = file.path(
#                     outDir, paste0(clusters[x], "_padj_", padj.cutoff, "_LRT_SigGenesUp.csv")),
#                     quote = FALSE,
#                     row.names = FALSE)
#         write.csv(res_wld_tb_down,
#                     file = file.path(
#                     outDir, paste0(clusters[x], "_padj_", padj.cutoff, "_LRT_SigGenesDown.csv")),
#                     quote = FALSE,
#                     row.names = FALSE)
        
#         # # Transform counts for data visualization
#         # rld <- rlog(dds_lrt, blind=TRUE)
        
#         # # Extract the rlog matrix from the object and compute pairwise correlation values
#         # rld_mat <- assay(rld)
#         # rld_cor <- cor(rld_mat)
        
        
#         # # Obtain rlog values for those significant genes
#         # cluster_rlog <- rld_mat[sigLRT_genes$gene, ]
        
#         # cluster_meta_sig <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
        
#         # # # Remove samples without replicates
#         # # cluster_rlog <- cluster_rlog[, -1]
#         # # cluster_metadata <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
        
        
#         # # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
#         # cluster_groups <- degPatterns(cluster_rlog, metadata = cluster_meta_sig, time = "group_id", col=NULL)
#         # ggsave(paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.png"))
        
#         # # Let's see what is stored in the `df` component
#         # write.csv(cluster_groups$df,
#         #         paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.csv"),
#         #         quote = FALSE, 
#         #         row.names = FALSE)
        
#         # saveRDS(cluster_groups, paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.rds"))
#         # save(dds_lrt, cluster_groups, res_LRT, sigLRT_genes, file = paste0("DESeq2/lrt/", clusters[x], "_all_LRTresults.Rdata"))
        
# }

# map(1:length(clusters), get_dds_LRTresults)


# ```

# ```{r}

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
#             ggplot(top20_sig_df, aes(y = value, x = Condition, col = Condition, label = SampleID)) +
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
#             ggplot(top20_sig_df, aes(y = value, x = Condition, col = Condition, label = SampleID)) +
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
## -------------------------------------------------------------------------------------------------
sessionInfo()

#' 
#' 
## -------------------------------------------------------------------------------------------------
knitr::purl(file.path(outDir, "../00_spinalcord_DESeq2_Wald.Rmd"), documentation = 2)

