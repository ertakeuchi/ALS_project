#' ---
#' title: "02 - FOR SYMPHONY SCTransform and DoubletFinder"
#' author: "ertakeuchi"
#' date: "14 August 2023"
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
#' ## 0-0. Load libraries
## ---- include=FALSE-----------------------------------------------------------
library("Seurat")
library("Signac")
library("tidyverse")
# library("DoubletFinder")
library("reticulate")
library("ggrepel")
use_python("/usr/bin/python3")
library("patchwork")
# library("leiden")
# library("future")  
# library("leiden")
set.seed(1234)


#' 
#' 
#' ## 0-1. Set a common directory
## -----------------------------------------------------------------------------
# homeDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq")
# NGSDir <- file.path(homeDir, "../../NGS_original")
# metadataDir <- file.path(homeDir, "metadata")
tissue <- "spinalcord"

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check if the correct number of arguments is provided
if (length(args) != 1) {
  stop("Usage: Rscript test.R tissue output_directory", call. = FALSE)
}
# Extract the tissue argument
rdataDir <- args[1]

# directory <- "230814_experiment4/out_4"
# outDir <- file.path(homeDir, "experimental_record", directory)
# rdataDir <- file.path(homeDir, "R/RData", tissue)


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
# |-----------------------------------------------------------------
# |-- 03: Nomalize using SCTransform                          <----- here
# |-----------------------------------------------------------------
# |-- 04: DoubletFinder and remove doublets
# |   |
# |   |-- 04-1: Check how many removed by DoubletFinder
# |   |-- 04-2: Save RData


#' 
#' 
#' ## Save Here for collaboraters
#' 
#' ## 3. Nomalize using SCTransform
## ---- echo=TRUE, results='hide', warning=FALSE, message=FALSE-----------------
# SCTransform: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
file_list <- list.files(path = rdataDir, pattern = "soListForSymphony", full.names = TRUE)
# file_list

soSCTList <- list()

for (i in seq_along(file_list)) {

  filename <- file_list[i]

  if (grepl("withoutRN01", filename)){
    name <- "withoutRN01"
  } else {
    name <- "includedRN01"
  }

  # print(name)

  so_list <- readRDS(file = file.path(filename))

  sct_list <- lapply(X = so_list, FUN = function(x){
    # add meta data
    x <- PercentageFeatureSet(x, "^RP[SL]", col.name = "percent.ribo")
    # do SCTransform
    x <- SCTransform(
      x,
      assay = "RNA",
      vst.flavor = "v2",
      min_cells = 0, # default = 5
      method = "glmGamPoi",
      vars.to.regress = c("percent.mt", "pct_chrY"),
      verbose = FALSE,
      return.only.var.genes = FALSE # This is important
      )
    x <- RunPCA(x, assay = "SCT", npcs = 30, verbose = FALSE)
  })

  soSCTList[[i]] <- sct_list
  names(soSCTList)[i] <- name
  
}

saveRDS(soSCTList, file = file.path(rdataDir, "soSCTList_ForSymphony.rds"))

#' 
#' 
## -----------------------------------------------------------------------------

