#' ---
#' title: "02 - SCTransform and DoubletFinder"
#' author: "ertakeuchi"
#' date: "04 August 2023"
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
library("tidyverse")
library("DoubletFinder")
library("reticulate")
library("ggrepel")
use_python("/usr/bin/python3")
library("patchwork")
library("leiden")
set.seed(1234)

#' 
#' 
#' ## 0-1. Set a common directory
## -----------------------------------------------------------------------------
homeDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq")
outDir <- file.path(homeDir, "experimental_record/230804_experiment2/out_2")
NGSDir <- file.path(homeDir, "../../NGS_original")
metadataDir <- file.path(homeDir, "metadata")

#' 
#' 
#' ################################################################################
#' # Preprocessing 02: SCTransform and DoubletFinder
#' ################################################################################
#' 
#' ## 1. Set a directory and prepare sample information
## -----------------------------------------------------------------------------
tissue <- "spinalcord"

rdataDir <- file.path(homeDir, "R/RData", tissue)

# get a Seurat object file list
file_list <- list.files(rdataDir, pattern = "_so.rds", full.name = FALSE)
# get a project name list
project_name_list <- lapply(file_list, function(file_name) {
  prefix <- gsub("_so.rds", "", file_name)
  return(prefix)
})
# read sample information
meta_data <- read.table(file.path(metadataDir, "SampleInformation.tsv"), header = TRUE, sep = "\t")
head(meta_data)
# meta_data$SubjectID <- as.factor(meta_data$SubjectID)


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

  # prepare sample information
  so@meta.data$Project <- project_name_list[[i]]
  meta <- (meta_data[meta_data$Project == project_name_list[[i]],])
  meta$Project <- NULL

  # add sample information
  meta_add <- so@meta.data %>%
    rownames_to_column(var = "rowname") %>%
    left_join(meta, by = "SubjectID") %>%
    column_to_rownames(var = "rowname")

  so <- AddMetaData(so, meta_add)

  # add chrY gene information
  genes_Y_selected <- genes_Y_table[genes_Y_table$gene_name %in% rownames(so), ]
  chrY.gene = genes_Y_selected$gene_name
  # calculate the percentage of chrY genes in each cell : Summed chrY read count / total raw read count
  so@meta.data$pct_chrY <- 
    Matrix::colSums(so@assays$RNA@counts[chrY.gene,])/colSums(so@assays$RNA@counts)
 
  # add Seurat object to a list
  so_list <- append(so_list, so)
}


#' 
#' 
#' ## 3. Nomalize using SCTransform
## ---- echo=TRUE, results='hide', warning=FALSE, message=FALSE-----------------
# SCTransform: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

sct_list <- lapply(X = so_list, FUN = function(x){
  # add meta data
  x <- PercentageFeatureSet(x, "^RP[SL]", col.name = "percent.ribo")
  # do SCTransform
  x <- SCTransform(
    x,
    assay = "RNA",
    vst.flavor = "v2",
    method = "glmGamPoi",
    vars.to.regress = c("percent.mt", "pct_chrY"),
    verbose = FALSE)
  x <- RunPCA(x, assay = "SCT", npcs = 30, verbose = FALSE)
  x <- RunUMAP(x, dims = 1:30, verbose = FALSE, reduction.name = "umap.rna")
})

so_list <- sct_list

#' 
#' ## 4. DoubletFinder and remove doublets
## ---- echo=TRUE, include=FALSE------------------------------------------------
# Count per cell :https://satijalab.org/costpercell/
# DoubletFinder: https://github.com/chris-mcginnis-ucsf/DoubletFinder

so_list_filtered <- list()
nExp_l <- list()
mpK_l <- list()

for (i in seq_along(so_list)){
  # i = 4
  so <- so_list[[i]]

  # preprocessing
  project <- so@meta.data$Project[i]
  m <- meta_data[meta_data$Project == project,]
  dbl <- as.numeric(unique(m$ExpectDoublet))

  # pK identification (no ground-truth)
  sweep.res.list <- paramSweep_v3(so, PCs = 1:30, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  nExp <- round(ncol(so) * dbl)  # expect "dbl[i]"% doublets
  nExp_l <- append(nExp_l, nExp)
  mpK_l <- append(mpK_l, mpK)

  so <- doubletFinder_v3(
    seu = so,
    PCs = 1:30,
    pN = 0.25, # default
    pK = mpK,
    nExp = nExp,
    sct = TRUE
    )

  so_list_filtered <- append(so_list_filtered, so) 
}

so_list <- list()
sod_list <- list()

for (i in seq_along(so_list_filtered)){

  # i = 4
  so <- so_list_filtered[[i]]
  nExp <- as.character(nExp_l[[i]])
  mpK <- mpK_l[[i]]

  doublet_state <- FetchData(so, vars = paste0("DF.classifications_0.25_", mpK, "_", nExp))
  so <- so[, which(x = doublet_state == "Singlet")]
  sod <- so[, which(x = doublet_state == "Doublet")]

  so_list <- append(so_list, so)
  sod_list <- append(sod_list, sod)
}

# Check how many removed by DoubletFinder
s_cells <- list()
m_cells <- list()

for (i in seq_along(so_list)){
  s_cell <- length(Cells(so_list[[i]]))
  m_cell <- length(Cells(sod_list[[i]]))
  s_cells <- append(s_cells, s_cell)
  m_cells <- append(m_cells, m_cell)
}


df <- data.frame(matrix(unlist(s_cells), length(s_cells), byrow=TRUE))
colnames(df) <- "Singlet"
df$Doublet <- unlist(m_cells)
df$ID <- rownames(df)

df <- gather(data = df, key = "doublet_status", value = "CellNumber", c("Singlet", "Doublet"))
print(df)

# 可視化
p <- ggplot(
  df, aes(x = ID, y = CellNumber, fill=doublet_status)) + 
  geom_bar(stat = 'identity', position = 'dodge')
# p

ggsave(
  file.path(outDir, paste0(tissue, "_AfterDoubletFinder.pdf")),
  height = 10, width = 21,
  plot = p)

# save RData
saveRDS(so_list, file = file.path(rdataDir, paste0(tissue, "_singlet_list.rds")))

#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../02_SCT_DoubletFinder.Rmd"), documentation = 2)

