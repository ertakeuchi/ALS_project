#' ---
#' title: "02 - Case vs Control DEG Analysis MAST: DOWN REGULATED"
#' author: "ertakeuchi"
#' date: "20 August 2023"
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
library("SeuratDisk")
library("Signac")
library("tidyverse")
use_python("/usr/bin/python3")
library("patchwork")
library("SingleCellExperiment")
library("Matrix")
library("SCOPfunctions")
library("lme4")
library("tidylo")
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
if (length(args) != 5) {
  stop("Usage: Rscript test.R tissue output_directory my_model case_name control_name", call. = FALSE)
}
# Extract the tissue argument
tissue <- args[1]
directory <- args[2]
my_model <- args[3]
case_name <- args[4]
control_name <- args[5]

# tissue <- "spinalcord"
# directory <- "230820_experiment2/out_2"

# case_name <- "HC"
# control_name <- "ALS"
# my_model <- "GLM"

outDir <- file.path(homeDir, "experimental_record", directory)
rdataDir <- file.path(homeDir, "R/RData", tissue)
resolutions <- c(0.2, 0.4, 0.6, 0.8)

l2fc <- 0.25
padj <- 0.2



#' 
#' ################################################################################
#' # DEG Analysis 02: MAST
#' ################################################################################
#' 
#' ## 1. Data Preprocessing
#' 
## -----------------------------------------------------------------------------
#' 
#' # 1. MAST with latent variables and random variables
#' 
#' ## 1-1. upregurated genes in case
## -----------------------------------------------------------------------------

subsets <- read.table(file = file.path(rdataDir, paste0(tissue, "_subsets.txt")), header = FALSE, stringsAsFactors = FALSE)$V1

soAllList <- list()
clustersList <- list()

for (i in seq_along(subsets)){
    # i = 1
    subPrimaryAnnotation <- subsets[[i]]

    so <- readRDS(file.path(rdataDir, paste0(subPrimaryAnnotation, "_so_withoutRN01_ThirdAnnotation.rds")))
    DefaultAssay(so) <- 'RNA'
    so <- NormalizeData(so)
    so$Sex <- factor(so$Sex)
    
    so$ThirdAnnotation_Condition <- paste0(so$ThirdAnnotation, "_", so$Condition)
    Idents(so) <- "ThirdAnnotation_Condition"
    cluster_id <- unique(so$ThirdAnnotation)

    soAllList[[i]] <- so
    names(soAllList)[[i]] <- subPrimaryAnnotation
    clustersList[[i]] <- cluster_id
    names(clustersList)[[i]] <- subPrimaryAnnotation
    

    }

# soAllList


#' 
#' 
## -----------------------------------------------------------------------------
#' Compute average log fold changes
#'
#' @param data feature * sample matrix
#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#' @param features Features to calculate fold change for.
#' If NULL, use all features
#' @importFrom Matrix rowSums
#' @rdname FoldChange
#' @concept differential_expression
#' @method FoldChange default
#' @references Copied from Seurat 3 source at https://github.com/satijalab/seurat/blob/master/R/differential_expression.R

.FoldChange.default <- function(
  data,
  cells.1,
  cells.2,
  mean.fxn,
  fc.name,
  features = NULL,
  ...
) {

  features <- if (is.null(features)) rownames(data) else features
  # Calculate percent expressed (各指定したクラスター内で、各遺伝子ごとの、発現している細胞の割合)
  thresh.min <- 0
  pct.1 <- round(
    x = rowSums(x = data[rownames(data) %in% features, colnames(data) %in% cells.1, drop = FALSE] > thresh.min) /
      length(x = cells.1),
    digits = 3
  )
  pct.2 <- round(
    x = rowSums(x = data[rownames(data) %in% features, colnames(data) %in% cells.2, drop = FALSE] > thresh.min) /
      length(x = cells.2),
    digits = 3
  )
  # Calculate fold change ("data"スロットの(LongNormalized data)平均値を比較)
  data.1 <- mean.fxn(data[rownames(data) %in% features, colnames(data) %in% cells.1, drop = FALSE])
  data.2 <- mean.fxn(data[rownames(data) %in% features, colnames(data) %in% cells.2, drop = FALSE])
  fc <- (data.1 - data.2)
  fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2))
  colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
  return(fc.results)
}

#' @param model Formula to use for MAST test
#' @param object Seurat >=3 object, with log-normalized data in the `data` slot
#' @param random_effect.vars Character vector of variables to add as random intersects in MAST model i.e. ~ ...  + (1|random_effect.var1) + (1|random_effect.var2)
#' @param ident.1 Identity class to define markers for
#' @param ident.2 A second identity class for comparison. Leave as NULL to compare with all other cells
#' @param cells.1 Vector of cell names belonging to group 1. Alternative way to specify ident.1
#' @param cells.2 Vector of cell names belonging to group 2. Alternative way to specify ident.2
#' @param logfc.threshold Only return results with a DE exceeding threshold
#' @param base base for log when computing mean and for output, default exp(1). NB: Seurat has changed to log2 in V4.
#' @param group.by Regroup cells into a different identity class prior to performing differential expression
#' @param assay Assay to pull data from; defaults to default assay
#' @param slot Slot to pull data from; defaults to "data" as MAST expects log-normalized data
#' @param features Genes to test. Default is to use all genes
#' @param min.pct Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
#' @param max.cells.per.ident Downsample each identity class to a max number. Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed to use for down sampling
#' @param latent.vars Variables to test
#' @param n_cores How many cores to use (temporarily sets options(mc.cores=n_cores))
#' @param verbose whether to print stdout from MAST functions
#' @param p.adjust.method passed to p.adjust. Note that n is all the genes in the object
#' @param \dots Additional parameters (other than formula, sca, method, ebayes, strictConvergence) to pass to MAST::zlm
#' @return data.frame with column names p_val, avg_log[base]FC, pct.1, pct.2, p_val_adj (identical format to Seurat::FindMarkers)
#' @export
#'
#' @references 
#' 
#' 

DE_MAST_RE_seurat = function(
  model,
  object,
  random_effect.vars = NULL,
  ident.1 = NULL,
  ident.2 = NULL,
  cells.1 = NULL,
  cells.2 = NULL,
  group.by = NULL,
  logfc.threshold = 0.25,
  base = exp(1),
  assay=NULL,
  slot="data",
  features = NULL,
  min.pct = 0.1,
  max.cells.per.ident = NULL,
  random.seed = 1,
  latent.vars = NULL,
  n_cores=NULL,
  verbose=TRUE,
  p.adjust.method="bonferroni",
  ...
) {


  if (!is.null(n_cores)) {
    # https://www.tidyverse.org/blog/2020/04/self-cleaning-test-fixtures/#the-onexit-pattern
    # set new value and capture old in op
    op <- options("mc.cores"=n_cores)
    on.exit(options(op), add = T)
  }

  if (is.null(c(ident.1, ident.2, cells.1, cells.2))) stop("at least one of ident.1, ident.2, cells.1, cells.2 must be provided")
  if (!is.null(ident.1) & !is.null(cells.1)) stop("Provide one of ident.1 or cells.1 but not both")
  if (!is.null(ident.2) & !is.null(cells.2)) stop("Provide one of ident.2 or cells.2 but not both")

  fc.name  <- if (base == exp(1)) "avg_logFC" else paste0("avg_log", base, "FC")
  #======== check inputs ========================================

  if(model == "GLMM"){
  stopifnot(!is.null(random_effect.vars))
  stopifnot(all(random_effect.vars %in% colnames(object@meta.data)))
  } else if(model == "GLM"){
    # do nothing
  }

  if (slot != "data") warning(paste0("MAST uses the logNormalised counts which are usually in the 'data' slot, but you are using ",slot))

  if (verbose) message(paste0("using ", round(base,2), " as log base. Make sure that the data slot has been log-transformed using this base"))
  #======== resolve idents and get cells ===================

  anyNA = if (!is.null(group.by)) any(is.na(object@meta.data[[group.by]])) else any(is.na(Seurat::Idents(object)))
  if (anyNA) stop("Some identities are NA, please check the metadata")

  logical.cells.1 = if (!is.null(cells.1)) colnames(object) %in% cells.1 else { if (!is.null(group.by)) {object@meta.data[[group.by]] == ident.1} else {Seurat::Idents(object) == ident.1}}
  if (sum(logical.cells.1)==0) {
    stop("no cells found matching ident.1. Did you forget to set Idents(object) or use group.by?")
  }
  if (is.null(cells.1)) cells.1 = colnames(object)[logical.cells.1]

  logical.cells.2 = if (!is.null(cells.2)) colnames(object) %in% cells.2 else {if (is.null(ident.2)) {!logical.cells.1} else {if (!is.null(group.by)) {object@meta.data[[group.by]] == ident.2} else {Seurat::Idents(object) == ident.2}}}
  if (sum(logical.cells.2)==0) {
    stop("no cells found matching ident.2. Did you forget to set Idents(object) or use group.by?")
  }
  if (is.null(cells.2)) cells.2 = colnames(object)[logical.cells.2]

  #======== features =================================

  features <- if (is.null(features)) rownames(object) else features

  if (any(!features %in% rownames(object))) {
    warning(paste0(paste(features[!features %in% rownames(object)], collapse = ", "), " not found in data!"))
  }

  features = features[features %in% rownames(object)]
  vec_logical_features = rep(T, length(features))

  data = Seurat::GetAssayData(object = object, assay=assay, slot=slot)

  densemat = utils_big_as.matrix(data, n_slices_init = 1, verbose=F)

  # compute average log fold change
  pseudocount.use = 1
  mean.fxn = function(x) {
    return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
  }
  # outputs a data.frame with columns fc.name, "pct.1", "pct.2"
  fc.results = .FoldChange.default(
      data = densemat,
      cells.1 = cells.1,
      cells.2 = cells.2,
      mean.fxn = mean.fxn,
      fc.name = fc.name,
      features = features
      )

  # ここでlog2fc > 0.25の遺伝子だけ論理型を逆転
  if (!is.null(logfc.threshold)) {
    vec_logical_features = vec_logical_features & fc.results[[fc.name]] > logfc.threshold  # need to convert back to non-log space to take mean
  }
  if (sum(vec_logical_features)<2) stop(paste0("Fewer than two features have a log fold change above ", logfc.threshold))
  # ここでmin.pctがどちらかの集団でマッチする遺伝子のみ論理型を逆転
  if (!is.null(min.pct)) {
    vec_logical_features = vec_logical_features & (fc.results$pct.1 > min.pct | fc.results$pct.2 > min.pct)
  }
  if (sum(vec_logical_features)<2) stop(paste0("Fewer than two features are expressed in ", min.pct, " of cells in any condition"))

  features <- features[vec_logical_features]　# ここでlfc.threshold, min.pctによる遺伝子の絞り込みが完了

　# Matrixを指定の遺伝子のみに縮小
  densemat = densemat[rownames(densemat) %in% features, , drop=F]

  #======== down sample cells ==============================

  idx.cells.1 = which(logical.cells.1)
  idx.cells.2 = which(logical.cells.2)

  if (!is.null(max.cells.per.ident)) {
    if (sum(logical.cells.1) > max.cells.per.ident) {
      set.seed(randomSeed)
      idx.cells.1 = sample(x = idx.cells.1, size = max.cells.per.ident, replace = F)
    }
    if (sum(!logical.cells.1) > max.cells.per.ident) {
      set.seed(randomSeed)
      idx.cells.2 = sample(x = idx.cells.2, size = max.cells.per.ident, replace = F)
    }
  }

  cells.1 = colnames(densemat)[idx.cells.1]
  cells.2 = colnames(densemat)[idx.cells.2]
  densemat = densemat[,c(idx.cells.1,idx.cells.2), drop=F]

  #======== filter out all zero features (if not already removed) ==

  densemat = densemat[apply(X = densemat, MARGIN=1, FUN = sum)>0,]

  #======== prep column (cell/sample) data =================

  df_coldata = data.frame(
    row.names = c(cells.1, cells.2),
    "wellKey" = c(cells.1, cells.2)
    ) # wellKey is hardcoded feature name in MAST

  # DE group factor
  df_coldata[cells.1, "group"] <- "Group1"
  df_coldata[cells.2, "group"] <- "Group2"
  df_coldata[, "group"] <- factor(x = df_coldata[, "group"])

  # latent vars
  if (!is.null(latent.vars)) {
    for (latent.var in latent.vars) {
      df_coldata[[latent.var]] <- object@meta.data[c(idx.cells.1,idx.cells.2), latent.var]
    }
  }
  
  if(model == "GLMM"){
  # add random vars to column data
  for (random_effect.var in random_effect.vars) {
    df_coldata[[random_effect.var]] <- as.factor(
        object@meta.data[c(idx.cells.1,idx.cells.2),random_effect.var]
      )
  }
  } else if(model == "GLM"){
    # do nothing
  }

  # check for NAs in column data
  for (colname in colnames(df_coldata)) {
    if (any(is.na(df_coldata[[colname]]))) {
      stop(paste0("variable '", colname, "' contains NAs"))
    }
  }

  if(model == "GLMM"){
  # check whether random var levels overlap entirely with test condition
  for (random_effect.var in random_effect.vars) {
    for (re_lvl in levels(df_coldata[[random_effect.var]])) {
      if (all(df_coldata[[random_effect.var]]==re_lvl & 1:ncol(densemat) %in% idx.cells.1)) {
        stop(paste0(random_effect.var, " level ", re_lvl, " overlaps entirely with ident.1. Increase ", max.cells.per.ident, " or check the data"))
      }
      else if (all(df_coldata[[random_effect.var]]==re_lvl & 1:ncol(densemat) %in% idx.cells.2)) {
        stop(paste0(random_effect.var, " level ", re_lvl, " overlaps entirely with ident.2. Increase ", max.cells.per.ident, " or check the data"))
      }
    }
  }
  } else if(model == "GLM"){
    # do nothing
  }

  #======== prep  feature data =========================

  df_feature = data.frame("primerid"=rownames(densemat)) # primerid is hardcoded feature name in MAST

  #======== prep SingleCellAssay object ============================

  expr = quote(MAST::FromMatrix(exprsArray=densemat,
                                 check_sanity = TRUE,
                                 cData=df_coldata,
                                 fData=df_feature))
  sca <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  SummarizedExperiment::colData(sca)$group <- relevel(SummarizedExperiment::colData(sca)$group, ref="Group2")

  #======== run MAST test ======================

  str_plus = if (!is.null(latent.vars)) " + " else ""

  # formula
  if(model == "GLMM"){
    fmla <- as.formula(
      object = paste0(" ~ group", str_plus, paste(latent.vars, collapse = " + "), paste(" + (1 |", random_effect.vars, ")", collapse=""))
    )
    print(fmla)
    # fit model parameters METHOD = glmer
    expr = quote(MAST::zlm(formula = fmla,
                    sca = sca,
                    method = 'glmer',
                    ebayes = F,　# if TRUE, regularize variance using empirical bayes method
                    strictConvergence = FALSE))
                    #  ...))

  } else if(model == "GLM"){
    fmla <- as.formula(
      object = paste0(" ~ group", str_plus, paste(latent.vars, collapse = " + "))
    )
    print(fmla)
    # fit model parameters METHOD = bayesglm
    expr = quote(MAST::zlm(formula = fmla,
                   sca = sca,
                   method = 'bayesglm'))
  }

  zlmCond <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  if (verbose) {
    print(zlmCond)
    message("Compute likelihoods and p-values")
  }
  # call a likelihood ratio test on the fitted object
  expr = quote(MAST::summary(object = zlmCond, doLRT = 'groupGroup1'))
  summaryCond <- if (verbose) eval(expr) else suppressWarnings(suppressMessages(eval(expr)))

  if ("character" %in% class(summaryCond)) stop("No differentially genes detected")

  de.results <- data.frame(
    "p_val" = summaryCond$datatable[contrast=='groupGroup1' & component=='H', `Pr(>Chisq)`],
    fc.results[vec_logical_features,])  #setDF(summaryCond$datatable[contrast=='groupGroup1' & component=='logFC', .(coef)])

  de.results$p_val_adj = p.adjust(de.results$p_val, method=p.adjust.method, n=length(vec_logical_features))

  de.results = de.results[order(de.results$p_val, -de.results[[fc.name]]),]

  return(de.results)

}

#' 
## -----------------------------------------------------------------------------
resultsList <- list()

for (i in seq_along(subsets)){
    # i = 1
    subPrimaryAnnotation <- subsets[[i]]
    so <- soAllList[[subPrimaryAnnotation]]

    result_list <- list()
    cluster_ids <- levels(as.factor(so$ThirdAnnotation))

    for (j in seq_along(cluster_ids)){
      # j = 2
      ThirdAnnotation <- cluster_ids[[j]]
      de.results <- DE_MAST_RE_seurat(
        model = my_model,
        object = so,
        # random_effect.vars=c("SubjectID"),
        ident.1 = paste(ThirdAnnotation, case_name, sep = "_"),
        ident.2 = paste(ThirdAnnotation, control_name, sep = "_"),
        group.by = "ThirdAnnotation_Condition",
        base = 2,
        latent.vars = c("Age", "Sex"),
        n_cores = 14
        ) 
      de.results$ThirdAnnotation <- ThirdAnnotation
      result_list[[j]] <- de.results

    }
    result_df <- do.call(rbind, result_list)
    resultsList[[i]] <- result_df
    names(resultsList)[[i]] <- subPrimaryAnnotation

}

resultsDf <- do.call(rbind, resultsList)

write.table(
  resultsDf,
  file = file.path(outDir, paste(my_model, case_name, "vs", control_name, "DEresults.txt", sep = "_")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  )


#' 
#' 
#' # 05. Count significant genes
## -----------------------------------------------------------------------------
resultsList <- list()
allClusters <- unlist(clustersList)

df <- read.table(
    file = file.path(
      outDir, paste(my_model, case_name, "vs", control_name, "DEresults.txt", sep = "_")),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
    )
head(df)

for (i in seq_along(allClusters)){
  # i = 1
  ThirdAnnotation <- allClusters[[i]]
  d <- df[df$ThirdAnnotation == ThirdAnnotation,] %>% na.omit()
  sig <- d[d$p_val_adj < padj & abs(d$avg_log2FC) >= l2fc,]
  print(sig)
  
  if(nrow(sig) == 0){
    print(paste0(ThirdAnnotation, " No significant genes"))
    sig[1, ] <- rep(NA, ncol(sig))
    sig$ThirdAnnotation <- ThirdAnnotation
  } else{
    rownames(sig) <- NULL
  }
  resultsList[[i]] <- sig
  }

data <- do.call(rbind, resultsList)
# head(data)

# ThirdAnnotationごとに行数をカウント
result <- data %>%
  group_by(ThirdAnnotation) %>%
  summarize(!!my_model := sum(!is.na(p_val_adj)), .groups = "keep")
# print(result)

write.table(
  result,
  file = file.path(outDir, paste0(my_model, "_", case_name, "vs", control_name, "_L2FC", l2fc, "_Padj",padj, "_DEResultsCount.txt")),
  sep = "\t",
  quote = FALSE
  )


#' 
#' ## SessionInfo
## -----------------------------------------------------------------------------
sessionInfo()

#' 
#' 
## -----------------------------------------------------------------------------
# knitr::purl(file.path(outDir, "../02.2_spinalcord_MAST_GLM.Rmd"), documentation = 2)

