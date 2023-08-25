#' ---
#' title: "02 - SCTransform"
#' author: "ertakeuchi"
#' date: "17 August 2023"
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
library("Signac")
library("tidyverse")
library("reticulate")
library("ggrepel")
use_python("/usr/bin/python3")
library("patchwork")
library("leiden")
library("Azimuth")
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

tissue <- "brain"
directory <- "230817_experiment3/out_3"
outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)


#' 
#' 
#' ################################################################################
#' # Preprocessing 02: SCTransform and DoubletFinder
#' ################################################################################
#' 
#' ## 1. Set a directory and prepare sample information
## -----------------------------------------------------------------------------
# 00: Set environment
# |
# |-- 01: Prepare sample information
# |   
# |-- 02: Add sample information and make a Seurat object list
# |   |
# |   |-- 02-1: Add chrY genes percent information
# |   
# |-- 03: Nomalize using SCTransform
# |   
# |-- 04: DoubletFinder and remove doublets
# |   |
# |   |-- 04-1: Check how many removed by DoubletFinder
# |   |-- 04-2: Save RData

# get a Seurat object file list
file_list <- list.files(rdataDir, pattern = "_so.rds", full.name = FALSE)
file_list <- file_list[!grepl("RN01", file_list)]
# file_list
# get a project name list
project_name_list <- lapply(file_list, function(file_name) {
  prefix <- gsub("_so.rds", "", file_name)
  return(prefix)
})
# read sample information
meta_data <- read.table(file.path(metadataDir, "SampleInformation.tsv"), header = TRUE, sep = "\t")
head(meta_data)
# meta_data$SubjectID <- as.factor(meta_data$SubjectID)


#' 
#' ## 2. Add sample information as meta_data and merge each Seurat object to a list
## -----------------------------------------------------------------------------
# Make chrY gene list
ref <-  rtracklayer::import(file.path(homeDir, "../../reference/Gencode_v37/gencode.v37.annotation.gtf")) %>% 
  as.data.frame()
# head(ref)
genes_Y_table <- ref[ref$seqnames == "chrY", c("seqnames", "start", "gene_name")] %>%
  distinct(gene_name, .keep_all = TRUE)
# head(genes_Y_table)

# Add sample information on every single Seurat object
so_list <- list()

for (i  in seq_along(file_list)){
  # i = 4
  # read each Seurat object
  so <- readRDS(file = file.path(rdataDir, file_list[i]))

  # add chrY gene information
  genes_Y_selected <- genes_Y_table[genes_Y_table$gene_name %in% rownames(so), ]
  chrY.gene = genes_Y_selected$gene_name
  # calculate the percentage of chrY genes in each cell : Summed chrY read count / total raw read count
  so@meta.data$pct_chrY <- 
    Matrix::colSums(so@assays$RNA@counts[chrY.gene,])/colSums(so@assays$RNA@counts)
 
  # add Seurat object to a list
  so_list <- append(so_list, so)
}

rm(ref, genes_Y_table)

#' 
#' 
#' ## 3. Nomalize using SCTransform
## ---- echo=TRUE, results='hide', warning=FALSE, message=FALSE-----------------
# SCTransform: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

reference <- LoadReference(path = file.path(homeDir, "../../reference/Azimuth/Human_M1_cortex"))

sct_list <- lapply(X = so_list, FUN = function(x){
  # add meta data
  x <- PercentageFeatureSet(x, "^RP[SL]", col.name = "percent.ribo")
  # do SCTransform
  x <- SCTransform(
    x,
    assay = "RNA",
    vst.flavor = "v2",
    new.assay.name = "refAssay",
    residual.features = rownames(x = reference$map),
    reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
    min_cells = 0, # default = 5
    method = "glmGamPoi",
    vars.to.regress = c("percent.mt", "pct_chrY"),
    verbose = FALSE,
    ncells = 2000, # default = 5000
    do.correct.umi = FALSE, # default = TRUE
    do.scale = FALSE, # default = FALSE
    do.center = TRUE # default = TRUE
    )
})

sct_list
# Merge Seurat objects
prefixes <- lapply(
  sct_list, function(x) unique(x@meta.data$Project)
) %>% unlist()

query_obj <- merge(
        sct_list[[1]],
        y = sct_list[2:length(sct_list)],
        project = "Azimuth",
        add.cell.ids = prefixes,
        merge.data = TRUE
        )

saveRDS(query_obj, file = file.path(outDir, "query_obj.rds"))

#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
#' 
## -----------------------------------------------------------------------------
knitr::purl(file.path(outDir, "../02.1_SCTransformForAzimuth.Rmd"), documentation = 2)

