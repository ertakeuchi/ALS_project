#' ---
#' title: "01 - Case vs Control DEG Analysis"
#' author: "ertakeuchi"
#' date: "28 August 2023"
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
#' # 01. Set environment
#' 
#' ## 0-0. Load libraries
## ---- include=FALSE-----------------------------------------------------------
library("Seurat")
library("tidyverse")
library("patchwork")
library("gridExtra")
library("ComplexHeatmap")
library("pheatmap")
library("RColorBrewer")
library("data.table")
set.seed(1234)


#' 
#' 
#' ## 0-1. Set a common directory (FOR CONTAINER)
## -----------------------------------------------------------------------------
homeDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq")
# homeDir <- file.path(".")
# NGSDir <- file.path(homeDir, "../../NGS_original")
metadataDir <- file.path(homeDir, "metadata")

# # Get command-line arguments
# args <- commandArgs(trailingOnly = TRUE)
# # Check if the correct number of arguments is provided
# if (length(args) != 2) {
#   stop("Usage: Rscript test.R tissue output_directory", call. = FALSE)
# }
# Extract the tissue argument
# tissue <- args[1]
# directory <- args[2]

tissue <- "spinalcord"
directory <- "230828_experiment1/out_1"
case_name <- "ALS"
control_name <- "HC"

outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)
resolutions <- c(0.2, 0.4, 0.6, 0.8)

#' 
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
   '#6D21https://vscode-remote+ssh-002dremote-002bfuji.vscode-resource.vscode-cdn.net/tmp/RtmpSxGaGx/vscode-R/plot.png?version%3D16922571488504F', '#182C61', '#FC427B', '#BDC581', '#82589F'
)

custom_colors$discrete <- c(colors_dutch, colors_spanish, colors_indian)

#' 
#' 
#' ###################################################################################
#' # DEG Analysis 00: DEGs Model Comparison
#' ###################################################################################
#' 
## -----------------------------------------------------------------------------
# 01: Set environment
# |
# |-- 02: Campared to each model 
# |   |
# |   |-- 01: Load MAST result
# |   |
# |   |-- 02: DESeq2 LRT
# |   |
# |   |-- 03: DESeq2 Wald
# |   
# |-- 03: Compared to each model


#' 
#' # 02. Campared to each model
#' 
#' ## 01. Load MAST result
## -----------------------------------------------------------------------------
# pattern <- paste0("_L2FC", l2fc, "_Padj", padj, "_DEResultsCount.txt")
pattern <- "DEResultsCount.txt"
dir <- file.path(outDir, "../..")
directory_list <- c(file.path(dir, "230820_experiment2/out_2"), file.path(dir, "230820_experiment3/out_3"))
directory_list

file_lists <- list()
for (i in seq_along(directory_list)){
  l <- list.files(directory_list[i], pattern = pattern, full.name = TRUE) 
  file_lists[[i]] <- l
}
file_lists <- unlist(file_lists)
# file_lists

files_filtered <- file_lists[grepl(paste0(case_name, "vs", control_name), file_lists)]
files_filtered

combined_list <- lapply(X = files_filtered, function(data){
  data <- read.table(file.path(data), header = TRUE, na.strings = 0)
})
l1 <- combined_list[[1]]
l2 <- combined_list[[2]]

merged_df <- left_join(l1, l2, by = "ThirdAnnotation")
head(merged_df)

merged_df[is.na(merged_df)] <- 0
head(merged_df)

# save
write.table(
  merged_df,
  file = file.path(outDir, "MAST_SumGenes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

#' 
#' ## 02. Load and modify DESeq2 Data
## -----------------------------------------------------------------------------
my_model <- "DESeq2_LRT"
directory_path <- file.path(outDir, "../../230824_experiment1/out_1")

clusters <- read.table(file.path(directory_path, "clusters.txt"))$V1
clusters

file_list <- list.files(directory_path, pattern = "LRT_SigGenes.csv", full.name = FALSE)
file_list

padj_values <- c(0.01, 0.05, 0.1, 0.2)
patterns <- padj_values

countsUp <- list()
countsDown <- list()

for (i in seq_along(patterns)){
  # i = 4
  pattern <- patterns[i]

  filtered_list <- file_list[grepl(pattern, file_list)]
  # filtered_list
  file_prefix_list <- lapply(filtered_list, function(filename) {
    elements <- unlist(strsplit(filename, "_"))
    prefix <- paste(elements[1], elements[2], sep = "_")
    return(prefix)
  })
  # file_prefix_list

  dataUp <- list()
  dataDown <- list()

    for (j in seq_along(clusters)){
        
      # j = 43
      ThirdAnnotation <- clusters[j]
      ThirdAnnotation

      if(ThirdAnnotation %in% file_prefix_list){

        data <- read.csv(file.path(directory_path, filtered_list[grepl(ThirdAnnotation, filtered_list)]), header = TRUE)
        head(data)

        dataU <- data %>% filter(log2FoldChange > 0)
        dataD <- data %>% filter(log2FoldChange < 0)
        # dataD
        
        datas <- list(dataU, dataD)

        data_list <- lapply(datas, FUN = function(x){

          if (nrow(x) >= 2) {
            x$ThirdAnnotation <- ThirdAnnotation
            } else {
            x[nrow(x)+1,] <- NA
            x2 <- data.frame(ThirdAnnotation = NA)
            x <- merge(x, x2, all = TRUE)
            x$ThirdAnnotation <- ThirdAnnotation
          }
          return(x)
        })
        data_list
        
        dataUp[[j]] <- data_list[[1]]
        dataDown[[j]] <- data_list[[2]]

      } else {

        data <- data.frame(
          gene = NA, baseMean = NA, log2FoldChange = NA,
          lfcSE = NA, stat = NA, pvalue = NA,
          padj = NA, ThirdAnnotation = ThirdAnnotation
          )
        dataUp[[j]] <- data
        dataDown[[j]] <- data
      }

    }

  dataUp <- do.call(rbind, dataUp)
  dataDown <- do.call(rbind, dataDown)
  
  head(dataUp)
  head(dataDown)

  for(l in seq_along(list(dataUp, dataDown))){

    data <- list(dataUp, dataDown)[[l]]
    sum_genes <- data %>% 
      group_by(ThirdAnnotation) %>% 
      summarize(!!paste0(my_model, pattern) := sum(!is.na(gene)))
      # summarize(GLM_old := sum(!is.na(gene)))
    # print(sum_genes, n = 43)
    head(sum_genes) 

    strings <- sum_genes$ThirdAnnotation
    strings <- strings[grepl("EX", strings)]
    numeric_part <- as.numeric(gsub("[^0-9]+", "", strings))
    sorted_strings <- strings[order(numeric_part)]
    indices <- which(clusters %in% sorted_strings)
    clusters[indices] <- sorted_strings
    sum_genes <- sum_genes[order(match(sum_genes$ThirdAnnotation, clusters)),]
    # head(sum_genes)
  
    if(l == 1){
      countsUp[[i]] <- sum_genes
    } else if(l ==2){
      countsDown[[i]] <- sum_genes
    }

  }
}

countsUp <- do.call(cbind, countsUp)
countsDown <- do.call(cbind, countsDown)

head(countsUp)
head(countsDown)

write.table(
  countsUp,
  file = file.path(outDir, paste0(my_model, "_AllSumGenesUp.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)

write.table(
  countsDown,
  file = file.path(outDir, paste0(my_model, "_AllSumGenesDown.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)

  for(num in seq_along(list(countsUp, countsDown))){

    if(num == 1){
      color <- "red"
      status <- "up"
      counts_df <- countsUp
    } else if(num == 2){
      color <- "blue"
      status <- "down"
      counts_df <- countsDown
    }

    df_long <- pivot_longer(data = counts_df, cols = -ThirdAnnotation)
    head(df_long)

    # Vizualize only LRT model
    p <- ggplot(df_long, aes(x = name, y = ThirdAnnotation, fill = value)) +
      geom_tile() +
      geom_text(aes(label = value), color = "black", size = 4) +
      scale_fill_gradient(low = "white", high = color) +
      theme_void() +
      labs(title = "Number of Significant Genes by Cell Type", 
      x = "gene", y = "") +
      coord_equal() +
      theme(axis.text.y=element_text(hjust=0), axis.title.y=element_blank(), 
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), 
            panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            plot.title=element_text(hjust=0.5), legend.position="none") 

    p <- p + theme(plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"))
    # p

    ggsave(
      plot = p,
      file = file.path( 
                    outDir,
                    paste0(my_model, "AllDEGsComparison", status, ".pdf")),
      width = 8, height = 20
    )
  }



#' 
#' 
#' ## 03. Load and modify DESeq2 Wald Data
## -----------------------------------------------------------------------------
my_model <- "DESeq2_Wald"
directory_path <- file.path(outDir, "../../230824_experiment3/out_3")

clusters <- read.table(file.path(directory_path, "../../230824_experiment1/out_1/clusters.txt"))$V1
# clusters

file_list <- list.files(directory_path, pattern = "_sig_genes.csv", full.name = FALSE)
file_list

padj_values <- c(0.01, 0.05, 0.1, 0.2)
patterns <- padj_values

countsUp <- list()
countsDown <- list()

for (i in seq_along(patterns)){
  # i = 4
  pattern <- patterns[i]

  filtered_list <- file_list[grepl(pattern, file_list)]
  filtered_list
  file_prefix_list <- lapply(filtered_list, function(filename) {
    elements <- unlist(strsplit(filename, "_"))
    prefix <- paste(elements[1], elements[2], sep = "_")
    return(prefix)
  })
  file_prefix_list

  dataUp <- list()
  dataDown <- list()

    for (j in seq_along(clusters)){
        
      # j = 39
      ThirdAnnotation <- clusters[j]
      ThirdAnnotation

      if(ThirdAnnotation %in% file_prefix_list){
        # print(TRUE)
        data <- read.csv(file.path(directory_path, filtered_list[grepl(ThirdAnnotation, filtered_list)]), header = TRUE)
        head(data)

        dataU <- data %>% filter(log2FoldChange > 0)
        dataD <- data %>% filter(log2FoldChange < 0)
        # dataU
        # dataD                  
        datas <- list(dataU, dataD)

        data_list <- lapply(datas, FUN = function(x){

          if (nrow(x) > 0) {
            x$ThirdAnnotation <- ThirdAnnotation
            } else {
            x[nrow(x)+1,] <- NA
            x2 <- data.frame(ThirdAnnotation = NA)
            x <- merge(x, x2, all = TRUE)
            x$ThirdAnnotation <- ThirdAnnotation
          }
          return(x)
        })
        data_list
        
        dataUp[[j]] <- data_list[[1]]
        dataDown[[j]] <- data_list[[2]]

     
      } else {
        # print(FALSE)
        data <- data.frame(
          gene = NA, baseMean = NA, log2FoldChange = NA,
          lfcSE = NA, pvalue = NA,
          padj = NA, ThirdAnnotation = ThirdAnnotation
          )
        data_list <- list(data, data)
        dataUp[[j]] <- data_list[[1]]
        dataDown[[j]] <- data_list[[2]]
      }

    }

  dataUp <- do.call(rbind, dataUp)
  dataDown <- do.call(rbind, dataDown)
  
  head(dataUp)
  head(dataDown)

  dataUp_short <- dataUp[,c("gene", "log2FoldChange", "padj", "ThirdAnnotation")]
  dataDown_short <- dataDown[,c("gene", "log2FoldChange", "padj", "ThirdAnnotation")]

  write.table(dataUp_short, file = file.path(outDir, paste0(my_model, "_SigGenesUp.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(dataDown_short, file = file.path(outDir, paste0(my_model, "_SigGenesDown.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  for(l in seq_along(list(dataUp, dataDown))){

    data <- list(dataUp, dataDown)[[l]]
    # l = 1
    sum_genes <- data %>% 
      group_by(ThirdAnnotation) %>% 
      summarize(!!paste0(my_model, pattern) := sum(!is.na(gene)))
      # summarize(GLM_old := sum(!is.na(gene)))
    # print(sum_genes, n = 43)
    head(sum_genes) 

    strings <- sum_genes$ThirdAnnotation
    strings <- strings[grepl("EX", strings)]
    numeric_part <- as.numeric(gsub("[^0-9]+", "", strings))
    sorted_strings <- strings[order(numeric_part)]
    indices <- which(clusters %in% sorted_strings)
    clusters[indices] <- sorted_strings
    sum_genes <- sum_genes[order(match(sum_genes$ThirdAnnotation, clusters)),]
    # head(sum_genes)
  
    if(l == 1){
      countsUp[[i]] <- sum_genes
    } else if(l == 2){
      countsDown[[i]] <- sum_genes
    }

  }
}

countsUp <- do.call(cbind, countsUp)
countsDown <- do.call(cbind, countsDown)

head(countsUp)
head(countsDown)

write.table(
  countsUp,
  file = file.path(outDir, paste0(my_model, "_AllSumGenesUp.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)

write.table(
  countsDown,
  file = file.path(outDir, paste0(my_model, "_AllSumGenesDown.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)

  for(num in seq_along(list(countsUp, countsDown))){

    if(num == 1){
      color <- "red"
      status <- "up"
      counts_df <- countsUp
    } else if(num == 2){
      color <- "blue"
      status <- "down"
      counts_df <- countsDown
    }

    df_long <- pivot_longer(data = counts_df, cols = -ThirdAnnotation)
    head(df_long)

    # Vizualize only LRT model
    p <- ggplot(df_long, aes(x = name, y = ThirdAnnotation, fill = value)) +
      geom_tile() +
      geom_text(aes(label = value), color = "black", size = 4) +
      scale_fill_gradient(low = "white", high = color) +
      theme_void() +
      labs(title = "Number of Significant Genes by Cell Type", 
      x = "gene", y = "") +
      coord_equal() +
      theme(axis.text.y=element_text(hjust=0), axis.title.y=element_blank(), 
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), 
            panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            plot.title=element_text(hjust=0.5), legend.position="none") 

    p <- p + theme(plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"))
    # p

    ggsave(
      plot = p,
      file = file.path( 
                    outDir,
                    paste0(my_model, "AllDEGsComparison", status, ".pdf")),
      width = 8, height = 20
    )
  }


#' 
#' 
#' # 03. Compared to each model
## -----------------------------------------------------------------------------
directory_path <- outDir
pattern <- "SumGenes"

file_list <- list.files(directory_path, pattern = pattern, full.name = FALSE)
file_list <- file_list[!grepl("mergedSumGenes", file_list)]
file_list

filter_up <- file_list[!grepl("Down", file_list)]
filter_down <- file_list[grepl("Down", file_list)]

up_list <- lapply(X = filter_up, function(data){
  data <- read.table(file.path(directory_path, data), header = TRUE, na.strings = 0)
})
down_list <- lapply(X = filter_down, function(data){
  data <- read.table(file.path(directory_path, data), header = TRUE, na.strings = 0)
})
up_list

for (i in seq_along(list(up_list, down_list))){

  # i = 1
  list_status <- list(up_list, down_list)[[i]]

  if(i == 1){
    color <- "red"
    status <- "up"
  } else if(i == 2){
    color <- "blue"
    status <- "down"
    }

  merged_df <- purrr::reduce(list_status, left_join, by = "ThirdAnnotation")
  head(merged_df)
  merged_df <- merged_df[, !grepl("^ThirdAnnotation\\.", names(merged_df))]

  merged_df[is.na(merged_df)] <- 0
  head(merged_df)
  # save
  write.table(
    merged_df,
    file = file.path(outDir, paste0(status, "_mergedSumGenes.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE)


  merged_long <- pivot_longer(data = merged_df, cols = -ThirdAnnotation)
  head(merged_long)

  # wide to long
  set_colnames <- as.vector(colnames(merged_df)[-1])
  mergedLong_df <- pivot_longer(data = merged_df, cols = set_colnames)
  head(mergedLong_df)
  # set_colnames <- c("GLMM", "GLM", "DESeq2_LRT")
  mergedLong_df$name <- factor(mergedLong_df$name)

  # ヒートマップを描画する
  p <- ggplot(merged_long, aes(x = name, y = ThirdAnnotation, fill = value)) +
    geom_tile() +
    geom_text(aes(label = value), color = "black", size = 4) +
    scale_fill_gradient(low = "white", high = color) +
    theme_void() +
    labs(title = paste0("Number of Significant Genes by Cell Type ", status), 
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
                  paste0(status, "_AllDEGsComparisonHeatmap.pdf")),
    width = 8, height = 20
  )
}


#' ## SessionInfo
## -----------------------------------------------------------------------------
sessionInfo()

#' 
#' 
## -----------------------------------------------------------------------------
knitr::purl(file.path(outDir, "../DESeq2_LRT.Rmd"), documentation = 2)

