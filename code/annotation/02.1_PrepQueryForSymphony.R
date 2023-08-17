#' ---
#' title: "02 - FOR SYMPHONY SCTransform and DoubletFinder WITHOUT RN01"
#' author: "ertakeuchi"
#' date: "16 August 2023"
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
library("DoubletFinder")
library("reticulate")
library("ggrepel")
use_python("/usr/bin/python3")
library("patchwork")
library("leiden")
library("future")
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
# rdataDir <- args[3]
# outputDir <- args[4]

# tissue <- "spinalcord"
# directory <- "230816_experiment2/out_2"
gse <- "GSE190442"
outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)

#' 
#' # 3. Load query data and MapQuery
## -----------------------------------------------------------------------------
# Load query data
soSCTList <- readRDS(file = file.path(rdataDir, "soSCTList_ForSymphony.rds"))
# names(soSCTList)
sct_list <- soSCTList[["withoutRN01"]]
rm(soSCTList)

so.features <- SelectIntegrationFeatures(
        object.list = sct_list
        ) # default HVG = 2000

sct_list <- PrepSCTIntegration(
        object.list = sct_list,
        anchor.features = so.features
        )

query_obj <- merge(
        sct_list[[1]],
        y = sct_list[2:length(sct_list)],
        project = gse,
        merge.data = TRUE
        )

saveRDS(query_obj, file = file.path(outDir, "query_obj.rds"))

#' 
#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../02_SCT_DoubletFinder.Rmd"), documentation = 2)

