#' ---
#' title: "Spinalcord DEGs Model Comparison"
#' author: "Eriko Takeuchi"
#' date: "20230801"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


#' 
#' # Packages
## ---- include=FALSE-----------------------------------------------------------
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(data.table)
set.seed(1234)

#' 
#' # inputs
#' ## set directory
## -----------------------------------------------------------------------------
mainDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq/R/out/final_out/spinalcord")
outDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq/experimental_record/230801_experiment3/out_3")
tissue <- "spinalcord"

#' 
#' ## Color setting
## -----------------------------------------------------------------------------
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
#' ###################################################################################
#' # Step 01: DEGs Model Comparison
#' ###################################################################################
#' 
#' ## 1-1.DEGs of GLM with latent.vars 
## -----------------------------------------------------------------------------
# ディレクトリのpathを指定
directory_path <- file.path(outDir, "../../230731_experiment2/out_2/GLMLatent")
# ディレクトリ内のファイル一覧を取得
file_list <- list.files(directory_path, pattern = "_glm_padj02_log2fc025_Up.csv", full.name = FALSE)
# ファイル名から接頭辞を取り除いたリストを作成
file_prefix_list <- lapply(file_list, function(file_name) {
  prefix <- gsub("_glm_padj02_log2fc025_Up.csv", "", file_name)
  return(prefix)
})

# データを格納する空のリストを作成
file_data <- list()

# ファイルを1つずつ処理してデータをリストに格納
for (i in seq_along(file_list)) {
  data <- read.csv(file.path(directory_path, file_list[i]), header = TRUE)
  if (nrow(data) >= 1) {
    data <- data[complete.cases(data), ]

     if (nrow(data) == 0){
        # 1行目としてNAを追加
        data[nrow(data)+1,] <- NA
        # geneとCellType_3rd列を追加
        data2 <- data.frame(CellType_3rd = NA)
        data <- merge(data, data2, all = TRUE)
     }
     
    # "CellType_3rd"列のみ、セルの値を変更
    data$CellType_3rd <- file_prefix_list[[i]]
    data$gene <- data$X
    data$X <- NULL
    } else {
    # 1行目としてNAを追加
    data[nrow(data)+1,] <- NA
    # geneとCellType_3rd列を追加
    data2 <- data.frame(CellType_3rd = NA)
    data <- merge(data, data2, all = TRUE)
    data$CellType_3rd <- file_prefix_list[[i]]
    data$gene <- data$X
    data$X <- NULL
  }
  file_data[[i]] <- data
}
# リスト内のデータフレームを1つのデータフレームに結合
combined_data <- do.call(rbind, file_data)
head(combined_data)

# CellType_3rdごとのgeneの合計を計算する
sum_genes_GlmLatent <- combined_data %>% 
  group_by(CellType_3rd) %>% 
  summarize(totalGenesGlmLatent = sum(!is.na(gene)))
print(sum_genes_GlmLatent, n = 43)

sum_genes_GlmLatent$CellType_3rd <- sapply(stringr::str_split(sum_genes_GlmLatent$CellType_3rd, pattern = "_", n = 4), 
    function(x) paste(x[1], x[2], sep = "_")) 
print(sum_genes_GlmLatent, n = 43)

write.table(
  sum_genes_GlmLatent,
  file = file.path(outDir, "sumGenesGlmLatent_Up.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

#' 
#' 
#' ## 1-2.DEGs of GLMM with latent.vars and random.vars 
## -----------------------------------------------------------------------------
# ディレクトリのpathを指定
directory_path <- file.path(outDir, "../../230731_experiment2/out_2/GLMM")
# ディレクトリ内のファイル一覧を取得
file_list <- list.files(directory_path, pattern = "_glmm_padj02_log2fc025_Up.csv", full.name = FALSE)
# ファイル名から接頭辞を取り除いたリストを作成
file_prefix_list <- lapply(file_list, function(file_name) {
  prefix <- gsub("_glmm_padj02_log2fc025_Up.csv", "", file_name)
  return(prefix)
})
# データを格納する空のリストを作成
file_data <- list()

# ファイルを1つずつ処理してデータをリストに格納
for (i in seq_along(file_list)) {
  data <- read.csv(file.path(directory_path, file_list[i]), header = TRUE)
  if (nrow(data) >= 1) {
    data <- data[complete.cases(data), ]

     if (nrow(data) == 0){
        # 1行目としてNAを追加
        data[nrow(data)+1,] <- NA
        # geneとCellType_3rd列を追加
        data2 <- data.frame(CellType_3rd = NA)
        data <- merge(data, data2, all = TRUE)
     }
     
    # "CellType_3rd"列のみ、セルの値を変更
    data$CellType_3rd <- file_prefix_list[[i]]
    data$gene <- data$X
    data$X <- NULL
    } else {
    # 1行目としてNAを追加
    data[nrow(data)+1,] <- NA
    # geneとCellType_3rd列を追加
    data2 <- data.frame(CellType_3rd = NA)
    data <- merge(data, data2, all = TRUE)
    data$CellType_3rd <- file_prefix_list[[i]]
    data$gene <- data$X
    data$X <- NULL
  }
  file_data[[i]] <- data
}
# リスト内のデータフレームを1つのデータフレームに結合
combined_data <- do.call(rbind, file_data)
head(combined_data)

# CellType_3rdごとのgeneの合計を計算する
sum_genes_GLMM <- combined_data %>% 
  group_by(CellType_3rd) %>% 
  summarize(totalGenesGlmm = sum(!is.na(gene)))
print(sum_genes_GLMM, n = 43)

sum_genes_GLMM$CellType_3rd <- sapply(stringr::str_split(sum_genes_GLMM$CellType_3rd, pattern = "_", n = 4), 
    function(x) paste(x[1], x[2], sep = "_")) 
print(sum_genes_GlmLatent, n = 43)

write.table(
  sum_genes_GLMM,
  file = file.path(outDir, "sumGenesGLMM_Up.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)


#' 
#' ## 1-3.DEGs of GLM without latent.vars 
## -----------------------------------------------------------------------------
# ディレクトリのpathを指定
directory_path <- file.path(outDir, "../../230731_experiment2/out_2/GLMNoLatent")
# ディレクトリ内のファイル一覧を取得
file_list <- list.files(directory_path, pattern = "_glm_padj02_log2fc025_Up.csv", full.name = FALSE)
# ファイル名から接頭辞を取り除いたリストを作成
file_prefix_list <- lapply(file_list, function(file_name) {
  prefix <- gsub("_glm_padj02_log2fc025_Up.csv", "", file_name)
  return(prefix)
})
file_prefix_list

# データを格納する空のリストを作成
file_data <- list()

# ファイルを1つずつ処理してデータをリストに格納
for (i in seq_along(file_list)) {
  data <- read.csv(file.path(directory_path, file_list[i]), header = TRUE)
  if (nrow(data) >= 1) {
    data <- data[complete.cases(data), ]

     if (nrow(data) == 0){
        # 1行目としてNAを追加
        data[nrow(data)+1,] <- NA
        # geneとCellType_3rd列を追加
        data2 <- data.frame(CellType_3rd = NA)
        data <- merge(data, data2, all = TRUE)
     }
     
    # "CellType_3rd"列のみ、セルの値を変更
    data$CellType_3rd <- file_prefix_list[[i]]
    data$gene <- data$X
    data$X <- NULL
    } else {
    # 1行目としてNAを追加
    data[nrow(data)+1,] <- NA
    # geneとCellType_3rd列を追加
    data2 <- data.frame(CellType_3rd = NA)
    data <- merge(data, data2, all = TRUE)
    data$CellType_3rd <- file_prefix_list[[i]]
    data$gene <- data$X
    data$X <- NULL
  }
  file_data[[i]] <- data
}
# リスト内のデータフレームを1つのデータフレームに結合
combined_data <- do.call(rbind, file_data)
head(combined_data)

# CellType_3rdごとのgeneの合計を計算する
sum_genes_GlmNoLatent <- combined_data %>% 
  group_by(CellType_3rd) %>% 
  summarize(totalGenesGlmNoLatent = sum(!is.na(gene)))
print(sum_genes_GlmNoLatent, n = 43)

sum_genes_GlmNoLatent$CellType_3rd <- sapply(stringr::str_split(sum_genes_GlmNoLatent$CellType_3rd, pattern = "_", n = 4), 
    function(x) paste(x[1], x[2], sep = "_")) 
print(sum_genes_GlmNoLatent, n = 43)

write.table(
  sum_genes_GlmNoLatent,
  file = file.path(outDir, "sumGenesGlmNoLatent_Up.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

#' 
#' ## 1-4.DEGs of DESeq2 (LRT) 
## -----------------------------------------------------------------------------
# ディレクトリのpathを指定
directory_path <- file.path(outDir, "../../230728_experiment1/out_1")
# ディレクトリ内のファイル一覧を取得
file_list <- list.files(directory_path, pattern = "_padj_0.2_LRT_SigGenesUp.csv", full.name = FALSE)
# ファイル名から接頭辞を取り除いたリストを作成
file_prefix_list <- lapply(file_list, function(file_name) {
  prefix <- gsub("_padj_0.2_LRT_SigGenesUp.csv", "", file_name)
  return(prefix)
})
file_prefix_list
# データを格納する空のリストを作成
file_data <- list()

# ファイルを1つずつ処理してデータをリストに格納
for (i in seq_along(file_list)) {
  # i = 1
  data <- read.csv(file.path(directory_path, file_list[i]), header = TRUE)
  head(data)
  
  if (nrow(data) >= 2) {

     if (nrow(data) == 0){
        # 1行目としてNAを追加
        data[nrow(data)+1,] <- NA
        # geneとCellType_3rd列を追加
        data2 <- data.frame(CellType_3rd = NA)
        data <- merge(data, data2, all = TRUE)
     }
     
    # "CellType_3rd"列のみ、セルの値を変更
    data$CellType_3rd <- file_prefix_list[[i]]
    } else {
    # 1行目としてNAを追加
    data[nrow(data)+1,] <- NA
    # geneとCellType_3rd列を追加
    data2 <- data.frame(CellType_3rd = NA)
    data <- merge(data, data2, all = TRUE)
    data$CellType_3rd <- file_prefix_list[[i]]
  }
  file_data[[i]] <- data
}

# リスト内のデータフレームを1つのデータフレームに結合
combined_data <- do.call(rbind, file_data)
head(combined_data)
nrow(combined_data)

# CellType_3rdごとのgeneの合計を計算する
sum_genes_DESeq2LRT <- combined_data %>% 
  group_by(CellType_3rd) %>% 
  summarize(totalGenesDESeq2LRT = sum(!is.na(gene)))
print(sum_genes_DESeq2LRT, n = 43)

write.table(
  sum_genes_DESeq2LRT,
  file = file.path(outDir, "sumGenesDESeq2LRT_Up.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

#' 
#' ## 1-5.DEGs of DESeq2 (Wald) 
## -----------------------------------------------------------------------------
# ディレクトリのpathを指定
directory_path <- file.path(outDir, "../../230801_experiment1/out_1")
# ディレクトリ内のファイル一覧を取得
file_list <- list.files(directory_path, pattern = "_Condition_CTRL_vs_ALS_WaldSigGenesUp.csv", full.name = FALSE)
# ファイル名から接頭辞を取り除いたリストを作成
file_prefix_list <- lapply(file_list, function(file_name) {
  prefix <- gsub("_Condition_CTRL_vs_ALS_WaldSigGenesUp.csv", "", file_name)
  return(prefix)
})
file_prefix_list
# データを格納する空のリストを作成
file_data <- list()

# ファイルを1つずつ処理してデータをリストに格納
for (i in seq_along(file_list)) {
  # i = 1
  data <- read.csv(file.path(directory_path, file_list[i]), header = TRUE)
  head(data)
  
  if (nrow(data) >= 2) {

     if (nrow(data) == 0){
        # 1行目としてNAを追加
        data[nrow(data)+1,] <- NA
        # geneとCellType_3rd列を追加
        data2 <- data.frame(CellType_3rd = NA)
        data <- merge(data, data2, all = TRUE)
     }
     
    # "CellType_3rd"列のみ、セルの値を変更
    data$CellType_3rd <- file_prefix_list[[i]]
    } else {
    # 1行目としてNAを追加
    data[nrow(data)+1,] <- NA
    # geneとCellType_3rd列を追加
    data2 <- data.frame(CellType_3rd = NA)
    data <- merge(data, data2, all = TRUE)
    data$CellType_3rd <- file_prefix_list[[i]]
  }
  file_data[[i]] <- data
}

# リスト内のデータフレームを1つのデータフレームに結合
combined_data <- do.call(rbind, file_data)
head(combined_data)
nrow(combined_data)

# CellType_3rdごとのgeneの合計を計算する
sum_genes_DESeq2Latent <- combined_data %>% 
  group_by(CellType_3rd) %>% 
  summarize(totalGenesDESeq2Wald = sum(!is.na(gene)))
print(sum_genes_DESeq2Latent, n = 43)

write.table(
  sum_genes_DESeq2Latent,
  file = file.path(outDir, "sumGenesDESeq2Wald_Up.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

#' 
#' 
#' # 2. combine DEGｓ dataframes
## -----------------------------------------------------------------------------
# ディレクトリのpathを指定
directory_path <- outDir
# ファイルをフィルタリングするための正規表現パターンを指定
pattern <- "^sum.*\\_Up.tsv$"
# ディレクトリ内のファイル一覧を取得
file_list <- list.files(directory_path, pattern = pattern, full.name = FALSE)
file_list

combined_list <- lapply(X = file_list, function(data){
  # data = "sumGenesGlmLatent_Up.tsv"
  data <- read.table(file.path(directory_path, data), header = TRUE, na.strings = 0)
})

# 全てのデータフレームを"CellType_3rd"でマージ
merged_df <- reduce(combined_list, left_join, by = "CellType_3rd")
head(merged_df)

merged_df[is.na(merged_df)] <- 0
head(merged_df)
# save
write.table(
  merged_df,
  file = file.path(outDir, "mergedSumGenes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

#' 
#' 
#' # 3. Visualization: Heat Plot 
## -----------------------------------------------------------------------------
merged_df <- read.table(file.path(outDir, "mergedSumGenes.tsv"), header = TRUE)
head(merged_df)

# wide to long
set_colnames <- as.vector(colnames(merged_df)[-1])
mergedLong_df <- pivot_longer(data = merged_df, cols = set_colnames)
head(mergedLong_df)

# ヒートマップを描画する
p <- ggplot(mergedLong_df, aes(x = name, y = CellType_3rd, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_void() +
  labs(title = "Number of Significant Genes by Cell Type", 
  x = "gene", y = "") +
  coord_equal() +
  theme(axis.text.y=element_text(hjust=0), axis.title.y=element_blank(), 
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        plot.title=element_text(hjust=0.5), legend.position="none") 

p <- p + theme(plot.margin = unit(c(1, 0.5, 1, 0.5), "cm"))
p

ggsave(
  plot = p,
  file = file.path( 
                outDir,
                paste0(tissue, "_Lfc025Padj02_DEGsComparisonHeatmap.pdf")),
  width = 8, height = 20
)

#' 
#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
#' 
## -----------------------------------------------------------------------------
knitr::purl(file.path(outDir, "../00_spinalcordDEGsModelComparison.Rmd"), documentation = 2)

