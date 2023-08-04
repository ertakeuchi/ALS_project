params <-
list(tissue_name = "spinalcord")

#' ---
#' title: "01 - cellranger count -> seurat object"
#' author: "ertakeuchi"
#' date: "04 August 2023"
#' output: 
#'   html_document:
#'     # code_folding: "hide"
#'     # fig_caption: "true"
#'   theme: 
#'     sketchy: "true"
#'     highlight: "tango"
#'   warnings: "false"
#' params: 
#'   tissue_name: 
#'     value: spinalcord
#'     choices:
#'       - spinalcord
#'       - brain
#' ---
#' 
## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)

#' 
#' ## 0-0. Load libraries
## ---- include=FALSE-----------------------------------------------------------
library("Signac")
library("Seurat")
library("EnsDb.Hsapiens.v86")
library("BSgenome.Hsapiens.UCSC.hg38")
library("tidyverse")
library("cowplot")
library("patchwork")
set.seed(1234)

#' 
#' ## 0-1. Set a common directory
## -----------------------------------------------------------------------------
homeDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq")
outDir <- file.path(homeDir, "experimental_record/230804_experiment1/out_1")
NGSDir <- file.path(homeDir, "../../NGS_original")
metadataDir <- file.path(homeDir, "metadata")

#' 
#' ## 0-2. Get gene annotations for hg38
## -----------------------------------------------------------------------------
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
annotation <- renameSeqlevels(annotation, mapSeqlevels(seqlevels(annotation), "UCSC"))
genome(annotation) <- "hg38"

#' 
#' ################################################################################
#' # Preprocessing 01: cellranger count -> seurat object
#' ################################################################################
#' 
#' ## 1. set directory and metadata
## -----------------------------------------------------------------------------
tissue <- "spinalcord"

rdataDir <- file.path(homeDir, "R/RData", tissue)
file_datas <- read.table(file.path(metadataDir, "FilePathInformation.tsv"), header = TRUE, sep = "\t")
file_data <- file_datas[file_datas$tissue == tissue, ]
head(file_data)

#' 
#' ## 2. cellranger count to seurat object
#' 
## ---- echo=FALSE, warning=FALSE, message=FALSE--------------------------------

# 00: Set environment
# |
# |-- 01: Data preprocessing on the FACS method
# |   |
# |   |-- 01-1: Load 10X data and create Seurat object
# |        |-- 01-1-1: Filtering with Hushtag Oligo
# |
# |-- 02: Data preprocessing on the density gradient method:Multiome
# |   |
# |   |-- 02-1: Load 10X data(Multiome) and create Seurat object
# |       |-- 02-1-1: DEMUXLET (RNA_POSTPROB > 0.99 & ATAC_POSTPROB > 0.99)
# |       |-- 02-1-2: skip DEMUXLET
# |
# |-- 03: Data preprocessing on the density gradient method:Gene Expression 
# |   |
# |   |-- 03-1: Load 10X data(Gene Expression) and create Seurat object
# |       |-- 03-1-1: DEMUXLET (RNA_POSTPROB > 0.99)
# |       |-- 03-1-2: skip DEMUXLET
# |
# |-- 04: Calculate Basic QC metrics
# |-- 05: Filter out low quality cells
# |   |
# |   |-- 05-1: Filter out on the Multiome data
# |   |-- 05-2: Filter out on the Gene Expression data
# |
# |-- 06: Visualize QC metrics and save RData


for (no in seq_along(file_data$project)){
  # no = 6
  print(paste0("Project ", file_data$project[no], " is being processed ..."))

  #------------------------------------------------------------------------------------------------#
  # 01: FACS METHOD - LOAD 10X DATA AND CREATE SEURAT OBJECT
  #------------------------------------------------------------------------------------------------#
  if(file_data$facs[no] == TRUE){

    print(paste0(file_data$project[no], " is the FACS data. Loading the FACS data ..."))

    data <- Read10X(data.dir = file.path(NGSDir, file_data$dir[no], "outs/filtered_feature_bc_matrix/"))
    data.hto <- data[["Antibody Capture"]]
    data.rna <- data[["Gene Expression"]]

    so <- CreateSeuratObject(counts = data.rna, min.cells = 0, min.features = 200)
    # create a new assay to store HTO information
    so[["HTO"]] <- CreateAssayObject(counts = data.hto[, colnames(so)])
    so <-  so %>% 
      NormalizeData(assay = "HTO", normalization.method = "CLR") %>%
      MULTIseqDemux(assay = "HTO", quantile = 0.7, qrange = seq)
    DefaultAssay(so) <- "HTO"
    # table(so$MULTI_classification)
    Idents(so) <- "MULTI_ID"
    RidgePlot(so, assay = "HTO", features = rownames(so[["HTO"]])[1:4], ncol = 2)
    ggsave(file = file.path(outDir, paste0(file_data$project[no], "_MULTIseqDemux.pdf")),
          height = 8, width = 8)
    # Filter out low possibility of MULTIseq
    multi_ids <- c("ALS-SPINALCORD-1", "HC-SPINALCORD-1")
    # multi_ids_brain <- c("ALS-BRAIN-1", "HC-BRAIN-1")
    so_filtered <- so %>% subset(MULTI_ID %in% multi_ids)
    DefaultAssay(so_filtered) <- "RNA"

    so_filtered@meta.data <- so_filtered@meta.data %>%
      dplyr::mutate(SubjectID = case_when(
        MULTI_ID == "HC-SPINALCORD-1" ~ "HC_6",
        MULTI_ID == "ALS-SPINALCORD-1" ~ "ALS_6",
        TRUE ~ as.character(MULTI_ID)  
      ))
    head(so_filtered@meta.data)
    so_filtered@meta.data$Project <- file_data$project[no]
    so_filtered@meta.data$SubjectID <- file_data$subject[no]
    so <- so_filtered
  } else{

    print(paste0(file_data$project[no], " is not the FACS data."))

  #------------------------------------------------------------------------------------------------#
  # 02: DENSITY GRADIENT METHOD - LOAD 10X DATA AND CREATE SEURAT OBJECT: MULTIOME
  #------------------------------------------------------------------------------------------------#

    if(file_data$multiome[no] == TRUE){

      print(paste0(file_data$project[no], " is the multiome data. Loading the multiome data ..."))

      counts <- Read10X_h5(file.path(NGSDir, file_data$dir[no], "outs/filtered_feature_bc_matrix.h5"))
      fragpath <- file.path(NGSDir, file_data$dir[no], "outs/atac_fragments.tsv.gz")

      # create a Seurat object containing the RNA data
      so_raw <- CreateSeuratObject(
        counts = counts$`Gene Expression`,
        assay = "RNA"
      )
      # "so" used later for filtering
      so_qc <- CreateSeuratObject(
        counts = counts$`Gene Expression`,
        assay = "RNA",
        min.cells = 0, min.features = 200
      )
      # only use peaks in standard chromosomes
      grange.counts <- StringToGRanges(rownames(counts$Peaks), sep = c(":", "-"))
      grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
      counts$Peaks <- counts$Peaks[as.vector(grange.use), ]

      # create ATAC assay and add it to the object
      so_raw[["ATAC"]] <-
        CreateChromatinAssay(
          counts = counts$Peaks,
          sep = c(":", "-"),
          min.cells = 10,
          fragments = fragpath,
          annotation = annotation
        )
    #------------------------------------------------------------------------------------------------#
    # 02-1: DENSITY GRADIENT METHOD - DEMUXLET: MULTIOME
    #------------------------------------------------------------------------------------------------#
      if (file_data$demuxlet[no] == TRUE){

        print(paste0(file_data$project[no], " has demuxlet output. Filtering by DEMUXLET ..."))
        # Read demuxlet data
        demuxlet_sng_rna <- read.table(
          file.path(NGSDir, file_data$dir[no], "outs", paste0(file_data$project[no], "_ALS_rna_demuxlet.single")),
          sep = "\t", header = TRUE)
        demuxlet_sng_atac <- read.table(
          file.path(NGSDir, file_data$dir[no], "outs", paste0(file_data$project[no], "_ALS_atac_demuxlet.single")),
          sep = "\t", header = TRUE)
        # head(demuxlet_sng_rna)
        demuxlet_highPB_rna <- demuxlet_sng_rna[which(demuxlet_sng_rna$POSTPRB > 0.99), ]
        demuxlet_highPB_atac <- demuxlet_sng_atac[which(demuxlet_sng_atac$POSTPRB > 0.99), ]
        demuxlet_highPB <- inner_join(demuxlet_highPB_rna, demuxlet_highPB_atac, by="BARCODE", suffix = c("_rna", "_atac"))
        # head(demuxlet_highPB)

        # Select cells which are matched between RNA and ATAC on DEMUXLET *SM_ID = SubjectID*
        meta_data <- data.frame(demuxlet_highPB[demuxlet_highPB$SM_ID_rna == demuxlet_highPB$SM_ID_atac, c("BARCODE", "SM_ID_rna")], row.names = 1)
        colnames(meta_data) <- "SubjectID"
        meta_data$Project <- file_data$project[no]
        head(meta_data)

        # Visualize the number of cells per subject
        cellnumber_df <- data.frame(table(meta_data$SubjectID))
        colnames(cellnumber_df) <- c("subject_id", "number")
        p <- ggplot(cellnumber_df, aes(x = subject_id, y = number), fill = subject_id) +
          geom_bar(stat = "identity") +
          theme(text = element_text(size = 30), 
            axis.text.x = element_text(angle = 90))
        ggsave(file = file.path(outDir, paste0(
          file_data$project[no], "_NumberOfCellsPerSubjectsAfterDemuxlet.pdf")),
            height = 14, width = 10,
            plot = p)

        so_raw <- AddMetaData(object = so_raw, metadata = meta_data)
        # head(so@meta.data)

        # Make a subset removed SubjectID == "NA"
        subject_ids <- unlist(strsplit(file_data$subject[no], ","))
        so_filtered <- so_raw %>% subset(SubjectID %in% subject_ids)
        so <- so_filtered

      } else {
        print(paste0(file_data$project[no], " does not have demuxlet output."))
        so_raw@meta.data$Project <- file_data$project[no]
        so_raw@meta.data$SubjectID <- file_data$subject[no]
        so <- so_raw
      }
    } else {
    #------------------------------------------------------------------------------------------------#
    # 03: DENSITY GRADIENT METHOD - LOAD 10X DATA AND CREATE SEURAT OBJECT: GENE EXPRESSION
    #------------------------------------------------------------------------------------------------#
      print(paste0(file_data$project[no], " is not the multiome data."))

      data <- Read10X(data.dir = file.path(NGSDir, file_data$dir[no], "outs/filtered_feature_bc_matrix/"))
      so <-  CreateSeuratObject(counts = data, min.cells = 0, min.features = 200)

    #------------------------------------------------------------------------------------------------#
    # 03-1: DENSITY GRADIENT METHOD - DEMUXLET: GENE EXPRESSION
    #------------------------------------------------------------------------------------------------#
      if (file_data$demuxlet[no] == TRUE){
        print(paste0(file_data$project[no], " has demuxlet output. Filtering by DEMUXLET ..."))
        # Read demuxlet data
        demuxlet_sng <- read.table(
          file.path(
            NGSDir, file_data$dir[no], "outs", paste0(
              file_data$project[no], "_ALS_rna_demuxlet.single")),
          sep = "\t", header = TRUE)

        # Select cells from the demuxlet output with high probability
        demuxlet_highPB <- demuxlet_sng[which(demuxlet_sng$POSTPRB > 0.99), ]
      
        # Add subject and project information to the Demuxlet table and filter cells with the demuxlet high probability
        barcode_df <- data.frame(BARCODE = colnames(so), Project = file_data$project[no])
        data_demuxlet <- right_join(barcode_df, demuxlet_highPB, by = "BARCODE")
        # change a column name "SM_ID" to "SubjectID"
        colnames(data_demuxlet)[colnames(data_demuxlet) == "SM_ID"] <- "SubjectID"
        # head(data_demuxlet)
        
        # Visualize the number of cells per subject
        cellnumber_df <- data.frame(table(data_demuxlet$SubjectID))
        colnames(cellnumber_df) <- c("subject_id", "number")
        p <- ggplot(cellnumber_df, aes(x = subject_id, y = number), fill = subject_id) +
          geom_bar(stat = "identity") +
          theme(text = element_text(size = 30), 
            axis.text.x = element_text(angle = 90))
        ggsave(file = file.path(outDir, paste0(
          file_data$project[no], "_NumberOfCellsPerSubjectsAfterDemuxlet.pdf")),
            height = 14, width = 10,
            plot = p)

        # Add Demuxlet information to the Seurat object
        meta_data <- data.frame(data_demuxlet[, c("BARCODE", "SubjectID")], row.names = 1) 
        # head(meta_data)
        so <- AddMetaData(object = so, metadata = meta_data)
        # head(so@meta.data)
        # table(so@meta.data$SubjectID)
      } else {
        print(paste0(file_data$project[no], " does not have demuxlet output."))
        so@meta.data$Project <- file_data$project[no]
        so@meta.data$SubjectID <- file_data$subject[no]
      }

    }
  }
  #------------------------------------------------------------------------------------------------#
  # 04: CALCULATE BASIC QC METRICS 
  #------------------------------------------------------------------------------------------------#

    so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

    p1 <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, raster = TRUE) &
      theme(axis.text.x = element_text(angle = 90))
    p2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = TRUE) + NoLegend()
    p3 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = TRUE) + NoLegend() 
    p <- (p2 + p3) / p1  
    p <- p + plot_annotation(
        title = paste0(
          "Project ", file_data$project[no], " before filtering"))
    # p
    ggsave(file = file.path(outDir, paste0(file_data$project[no], "_QCBeforeFiltering.pdf")),
          height = 16, width = 12,
          plot = p)

  #------------------------------------------------------------------------------------------------#
  # 05: FILTER OUT LOW QUALITY CELLS
  #------------------------------------------------------------------------------------------------#

  if(file_data$multiome == TRUE){

    so.QC <- subset(
      x = so,
      cells = colnames(so),
      subset = nFeature_RNA > 200 &
      nFeature_RNA < 10000 &
      nCount_ATAC < 100000 &
      nCount_ATAC > 1000 &
      percent.mt < 5
    )
  } else {
      so.QC <- subset(
        x = so,
        subset = nFeature_RNA > 200 &
        nFeature_RNA < 10000 &
        percent.mt < 5
        )
  }
  #------------------------------------------------------------------------------------------------#
  # 06: VISUALIZE QC METRICS AND SAVE RDATA
  #------------------------------------------------------------------------------------------------#

    p1 <- VlnPlot(so.QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, raster = TRUE) &
      theme(axis.text.x = element_text(angle = 90))
    p2 <- FeatureScatter(so.QC, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = TRUE) + NoLegend()
    p3 <- FeatureScatter(so.QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = TRUE) + NoLegend() 
    p <- (p2 + p3) / p1  
    p <- p + plot_annotation(
        title = paste0(
          "Project ", file_data$project[no], " before filtering"))
    # p
    ggsave(file = file.path(outDir, paste0(file_data$project[no], "_QCAfterFiltering.pdf")),
          height = 16, width = 12,
          plot = p)

    saveRDS(so.QC, file = file.path(rdataDir, paste0(file_data$project[no], "_so.rds")))

}

#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../01_createSeuratObject.Rmd"), documentation = 2)

