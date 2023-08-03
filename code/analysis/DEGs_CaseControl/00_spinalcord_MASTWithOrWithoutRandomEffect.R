#' ---
#' title: "Spinalcord MAST with/without random effect"
#' author: "Eriko Takeuchi"
#' date: "20230731"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' # Packages
## ---- include=FALSE-----------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
library(SCOPfunctions)
library(lme4)
library(tidylo)
# devtools::install_github("CBMR-Single-Cell-Omics-Platform/SCOPfunctions")
set.seed(1234)

#' 
#' # inputs
#' ## set directory
## -----------------------------------------------------------------------------
mainDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq/R/out/final_out/spinalcord")
outDir <- file.path("/mnt/home/etakeuchi/bioinformatics/ALS_snRNAseq/experimental_record/230731_experiment2/out_2")
tissue <- "spinalcord"

#' 
#' # 1. MAST with latent variables and random variables
#' 
#' ## 1-1. upregurated genes in case
## -----------------------------------------------------------------------------
sub_list <- c("neuron", "others", "astro", "micro", "oligos")

for (i in seq_along(sub_list)){
    # i = 1
    celltype <- sub_list[[i]]

    so.sub <- readRDS(
        file = file.path(
            mainDir,
            "subset",
            paste0("sosub_", toString(celltype), "_3rd_layer_prepSCT.rds")
            )
        )
    so.sub <- subset(so.sub, subset = Project == "RN01", invert = TRUE)


    celltype3rd_list <- unique(so.sub$CellType_3rd)

        for (j in seq_along(celltype3rd_list)){ 
        # j = 3
        celltype_3rd <- celltype3rd_list[[j]]
        print(celltype_3rd)

        DefaultAssay(so.sub) <- 'RNA'
        so.sub <- NormalizeData(so.sub)


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
#  # compute average log fold change（seuratのFindMarkersに入っている関数　expm1($num)はexp($num)-1）
#   pseudocount.use = 1
#   mean.fxn = function(x) {
#     return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
#   }

#' MAST differential expression test with random effect for Seurat object
#'
#' Uses Seurat::FindMarkers MAST test hurdle model with added random effects variables
#'
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
#' @references Zimmerman, K.D., Espeland, M.A. & Langefeld, C.D.
#' . A practical solution to pseudoreplication bias in single-cell studies.
#' . Nat Commun 12, 738 (2021). https://doi.org/10.1038/s41467-021-21038-1
#'
#' . McDavid A, Finak G, Yajima M (2020). MAST: Model-based Analysis of Single Cell Transcriptomics. R package version 1.16.0, https://github.com/RGLab/MAST/.
#'
#' . Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019)
#'
#' . Seurat 3 source at https://github.com/satijalab/seurat/blob/master/R/differential_expression.R

###-------------test--------------------------###
object = so.sub
random_effect.vars=c("SampleID")
# random_effect.vars = NULL
ident.1 = paste(celltype_3rd, "ALS", sep = "_")
ident.2 = paste(celltype_3rd, "CTRL", sep = "_")
group.by = "CellType_3rd_Condition"
base = 2
latent.vars = c("Age", "Sex")
# latent.vars = NULL
n_cores = 12

cells.1 = NULL
cells.2 = NULL
logfc.threshold = 0.25
assay=NULL
slot="data"
features = NULL
min.pct = 0.1
max.cells.per.ident = NULL
random.seed = 1
verbose=TRUE
# p.adjust.method="fdr"
p.adjust.method="bonferroni" # seuratのFindMarkersのデフォルト
###-------------test--------------------------###


# DE_MAST_RE_seurat = function(
#   object,
#   random_effect.vars,
#   ident.1 = NULL,
#   ident.2 = NULL,
#   cells.1 = NULL,
#   cells.2 = NULL,
#   group.by = NULL,
#   logfc.threshold = 0.25,
#   base = exp(1),
#   assay=NULL,
#   slot="data",
#   features = NULL,
#   min.pct = 0.1,
#   max.cells.per.ident = NULL,
#   random.seed = 1,
#   latent.vars = NULL,
#   n_cores=NULL,
#   verbose=TRUE,
#   p.adjust.method="fdr",
#   ...
# ) {

  # require(Seurat)
  # require(MAST)

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

  stopifnot(!is.null(random_effect.vars))
  stopifnot(all(random_effect.vars %in% colnames(object@meta.data)))

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
      data=densemat,
      cells.1=cells.1,
      cells.2=cells.2,
      mean.fxn=mean.fxn,
      fc.name=fc.name,
      features = features)

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

  # add random vars to column data
  for (random_effect.var in random_effect.vars) {
    df_coldata[[random_effect.var]] <- as.factor(
        object@meta.data[c(idx.cells.1,idx.cells.2),random_effect.var]
      )
  }

  # check for NAs in column data
  for (colname in colnames(df_coldata)) {
    if (any(is.na(df_coldata[[colname]]))) {
      stop(paste0("variable '", colname, "' contains NAs"))
    }
  }

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

  #======== prep  feature data =========================

  df_feature = data.frame("primerid"=rownames(densemat)) # primerid is hardcoded feature name in MAST

  #======== prep SingleCellAssay object ============================

  expr = quote(MAST::FromMatrix(exprsArray=densemat,
                                 check_sanity = TRUE,
                                 cData=df_coldata,
                                 fData=df_feature))
  sca <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  # set group as a factor and make group1 the reference
  #cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  # in Seurat source this is Group1!?
  # see https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
  SummarizedExperiment::colData(sca)$group <- relevel(SummarizedExperiment::colData(sca)$group, ref="Group2")
  # cond <- relevel(x = cond, ref = "Group2")
  # SummarizedExperiment::colData(sca)$condition <- cond

  # add the centered number of expressed genes to the model
  # cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
  # SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)

  #======== run MAST test ======================

  str_plus = if (!is.null(latent.vars)) " + " else ""

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

  # # METHOD = bayesglm
  # expr = quote(MAST::zlm(formula = fmla,
  #                  sca = sca,
  #                  method = 'bayesglm'))

  zlmCond <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  if (verbose) {
    print(zlmCond)
    message("Compute likelihoods and p-values")
  }
  # call a likelihood ratio test on the fitted object
  expr = quote(MAST::summary(object = zlmCond, doLRT = 'groupGroup1'))
  summaryCond <- if (verbose) eval(expr) else suppressWarnings(suppressMessages(eval(expr)))

  if ("character" %in% class(summaryCond)) stop("No differentially genes detected")
  # print(summaryCond,n=4)
  # Fitted zlm with top 4 genes per contrast:
  #   ( log fold change Z-score )
  # primerid           conditionGroup1 cngeneson
  # Akr1c20              -4.1            13.2*
  # Bclaf3               18.5*            4.8
  # Cmss1                -0.3            -9.9*
  # Cux2                -16.5*            4.6
  # ENSMUSG00000115022  -19.5*            7.6
  # Pot1b                13.5*            5.3
  # Siah1a                4.6            17.6*
  # Sipa1l2               1.8            10.5*
  #
  #
  # summaryCond$datatable has colnames
  # "primerid" feature name
  # "component" if no latent.vars, "C" (continuous) "D" (discrete)  "H" (hurdle) "S" (RE?) "logFC"
  # "contrast" if no latent.vars,  conditionGroup1 (Intercept) cngeneson
  # "ci.hi" confidence interval (high)
  # "ci.lo" confidence interval (low)
  # "coef" coefficient estimate
  # "z" z-scores

  # uses
  # https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#differential-expression-using-a-hurdle-model
  # and
  # https://github.com/kdzimm/PseudoreplicationPaper/blob/master/Type_1_Error/Type%201%20-%20MAST%20RE.Rmd

  de.results <- data.frame(
    "p_val" = summaryCond$datatable[contrast=='groupGroup1' & component=='H', `Pr(>Chisq)`],
    fc.results[vec_logical_features,])  #setDF(summaryCond$datatable[contrast=='groupGroup1' & component=='logFC', .(coef)])

  de.results$p_val_adj = p.adjust(de.results$p_val, method=p.adjust.method, n=length(vec_logical_features))

  de.results = de.results[order(de.results$p_val, -de.results[[fc.name]]),]

  # return(de.results)
  head(de.results)

  test <- de.results[(de.results$avg_log2FC >= 0.25 & de.results$p_val_adj < 0.2), ]
  # head(test, n = 20)
  # nrow(test)

  write.csv(
    test,
    file = file.path(outDir, "GLMM", paste0(celltype_3rd, "_glmm_padj02_log2fc025_Up.csv")),
    quote = FALSE
    )
  }
}


#' 
#' 
#' ## 1-2. downregurated genes in case
## -----------------------------------------------------------------------------
sub_list <- c("neuron", "others", "astro", "micro", "oligos")

for (i in seq_along(sub_list)){
    i = 4
    celltype <- sub_list[[i]]

    so.sub <- readRDS(
        file = file.path(
            mainDir,
            "subset",
            paste0("sosub_", toString(celltype), "_3rd_layer_prepSCT.rds")
            )
        )
    so.sub <- subset(so.sub, subset = Project == "RN01", invert = TRUE)


    celltype3rd_list <- unique(so.sub$CellType_3rd)
    celltype3rd_list

        for (j in seq_along(celltype3rd_list)){ 
        j = 6
        celltype_3rd <- celltype3rd_list[[j]]
        print(celltype_3rd)

        DefaultAssay(so.sub) <- 'RNA'
        so.sub <- NormalizeData(so.sub)


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
#  # compute average log fold change（seuratのFindMarkersに入っている関数　expm1($num)はexp($num)-1）
#   pseudocount.use = 1
#   mean.fxn = function(x) {
#     return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
#   }

#' MAST differential expression test with random effect for Seurat object
#'
#' Uses Seurat::FindMarkers MAST test hurdle model with added random effects variables
#'
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
#' @references Zimmerman, K.D., Espeland, M.A. & Langefeld, C.D.
#' . A practical solution to pseudoreplication bias in single-cell studies.
#' . Nat Commun 12, 738 (2021). https://doi.org/10.1038/s41467-021-21038-1
#'
#' . McDavid A, Finak G, Yajima M (2020). MAST: Model-based Analysis of Single Cell Transcriptomics. R package version 1.16.0, https://github.com/RGLab/MAST/.
#'
#' . Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019)
#'
#' . Seurat 3 source at https://github.com/satijalab/seurat/blob/master/R/differential_expression.R

###-------------test--------------------------###
object = so.sub
random_effect.vars=c("SampleID")
# random_effect.vars = NULL
# Downregulated genesのために、ident.1とident.2を逆にする
ident.1 = paste(celltype_3rd, "CTRL", sep = "_")
ident.2 = paste(celltype_3rd, "ALS", sep = "_")
group.by = "CellType_3rd_Condition"
base = 2
latent.vars = c("Age", "Sex")
# latent.vars = NULL
n_cores = 12

cells.1 = NULL
cells.2 = NULL
logfc.threshold = 0.25
assay=NULL
slot="data"
features = NULL
min.pct = 0.1
max.cells.per.ident = NULL
random.seed = 1
verbose=TRUE
# p.adjust.method="fdr"
p.adjust.method="bonferroni" # seuratのFindMarkersのデフォルト
###-------------test--------------------------###


# DE_MAST_RE_seurat = function(
#   object,
#   random_effect.vars,
#   ident.1 = NULL,
#   ident.2 = NULL,
#   cells.1 = NULL,
#   cells.2 = NULL,
#   group.by = NULL,
#   logfc.threshold = 0.25,
#   base = exp(1),
#   assay=NULL,
#   slot="data",
#   features = NULL,
#   min.pct = 0.1,
#   max.cells.per.ident = NULL,
#   random.seed = 1,
#   latent.vars = NULL,
#   n_cores=NULL,
#   verbose=TRUE,
#   p.adjust.method="fdr",
#   ...
# ) {

  # require(Seurat)
  # require(MAST)

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

  stopifnot(!is.null(random_effect.vars))
  stopifnot(all(random_effect.vars %in% colnames(object@meta.data)))

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
      data=densemat,
      cells.1=cells.1,
      cells.2=cells.2,
      mean.fxn=mean.fxn,
      fc.name=fc.name,
      features = features)

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

  # add random vars to column data
  for (random_effect.var in random_effect.vars) {
    df_coldata[[random_effect.var]] <- as.factor(
        object@meta.data[c(idx.cells.1,idx.cells.2),random_effect.var]
      )
  }

  # check for NAs in column data
  for (colname in colnames(df_coldata)) {
    if (any(is.na(df_coldata[[colname]]))) {
      stop(paste0("variable '", colname, "' contains NAs"))
    }
  }

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

  #======== prep  feature data =========================

  df_feature = data.frame("primerid"=rownames(densemat)) # primerid is hardcoded feature name in MAST

  #======== prep SingleCellAssay object ============================

  expr = quote(MAST::FromMatrix(exprsArray=densemat,
                                 check_sanity = TRUE,
                                 cData=df_coldata,
                                 fData=df_feature))
  sca <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  # set group as a factor and make group1 the reference
  #cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  # in Seurat source this is Group1!?
  # see https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
  SummarizedExperiment::colData(sca)$group <- relevel(SummarizedExperiment::colData(sca)$group, ref="Group2")
  # cond <- relevel(x = cond, ref = "Group2")
  # SummarizedExperiment::colData(sca)$condition <- cond

  # add the centered number of expressed genes to the model
  # cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
  # SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)

  #======== run MAST test ======================

  str_plus = if (!is.null(latent.vars)) " + " else ""

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

  # # METHOD = bayesglm
  # expr = quote(MAST::zlm(formula = fmla,
  #                  sca = sca,
  #                  method = 'bayesglm'))

  zlmCond <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  if (verbose) {
    print(zlmCond)
    message("Compute likelihoods and p-values")
  }
  # call a likelihood ratio test on the fitted object
  expr = quote(MAST::summary(object = zlmCond, doLRT = 'groupGroup1'))
  summaryCond <- if (verbose) eval(expr) else suppressWarnings(suppressMessages(eval(expr)))

  if ("character" %in% class(summaryCond)) stop("No differentially genes detected")
  # print(summaryCond,n=4)
  # Fitted zlm with top 4 genes per contrast:
  #   ( log fold change Z-score )
  # primerid           conditionGroup1 cngeneson
  # Akr1c20              -4.1            13.2*
  # Bclaf3               18.5*            4.8
  # Cmss1                -0.3            -9.9*
  # Cux2                -16.5*            4.6
  # ENSMUSG00000115022  -19.5*            7.6
  # Pot1b                13.5*            5.3
  # Siah1a                4.6            17.6*
  # Sipa1l2               1.8            10.5*
  #
  #
  # summaryCond$datatable has colnames
  # "primerid" feature name
  # "component" if no latent.vars, "C" (continuous) "D" (discrete)  "H" (hurdle) "S" (RE?) "logFC"
  # "contrast" if no latent.vars,  conditionGroup1 (Intercept) cngeneson
  # "ci.hi" confidence interval (high)
  # "ci.lo" confidence interval (low)
  # "coef" coefficient estimate
  # "z" z-scores

  # uses
  # https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#differential-expression-using-a-hurdle-model
  # and
  # https://github.com/kdzimm/PseudoreplicationPaper/blob/master/Type_1_Error/Type%201%20-%20MAST%20RE.Rmd

  de.results <- data.frame(
    "p_val" = summaryCond$datatable[contrast=='groupGroup1' & component=='H', `Pr(>Chisq)`],
    fc.results[vec_logical_features,])  #setDF(summaryCond$datatable[contrast=='groupGroup1' & component=='logFC', .(coef)])

  de.results$p_val_adj = p.adjust(de.results$p_val, method=p.adjust.method, n=length(vec_logical_features))

  de.results = de.results[order(de.results$p_val, -de.results[[fc.name]]),]

  # return(de.results)
  head(de.results)

  test <- de.results[(de.results$avg_log2FC >= 0.25 & de.results$p_val_adj < 0.2), ]
  # head(test, n = 20)
  # nrow(test)

  write.csv(
    test,
    file = file.path(outDir, "GLMM", paste0(celltype_3rd, "_glmm_padj02_log2fc025_Down.csv")),
    quote = FALSE
    )
  }
}


#' 
#' 
#' 
#' # 2. MAST with only latent variables
#' 
#' ## 2-1. Upregurated genes in case
## -----------------------------------------------------------------------------
sub_list <- c("neuron", "others", "astro", "micro", "oligos")

for (i in seq_along(sub_list)){
    # i = 1
    celltype <- sub_list[[i]]

    so.sub <- readRDS(
        file = file.path(
            mainDir,
            "subset",
            paste0("sosub_", toString(celltype), "_3rd_layer_prepSCT.rds")
            )
        )
    so.sub <- subset(so.sub, subset = Project == "RN01", invert = TRUE)


    celltype3rd_list <- unique(so.sub$CellType_3rd)

        for (j in seq_along(celltype3rd_list)){ 
        # j = 3
        celltype_3rd <- celltype3rd_list[[j]]
        print(celltype_3rd)

        DefaultAssay(so.sub) <- 'RNA'
        so.sub <- NormalizeData(so.sub)


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
#  # compute average log fold change（seuratのFindMarkersに入っている関数　expm1($num)はexp($num)-1）
#   pseudocount.use = 1
#   mean.fxn = function(x) {
#     return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
#   }

#' MAST differential expression test with random effect for Seurat object
#'
#' Uses Seurat::FindMarkers MAST test hurdle model with added random effects variables
#'
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
#' @references Zimmerman, K.D., Espeland, M.A. & Langefeld, C.D.
#' . A practical solution to pseudoreplication bias in single-cell studies.
#' . Nat Commun 12, 738 (2021). https://doi.org/10.1038/s41467-021-21038-1
#'
#' . McDavid A, Finak G, Yajima M (2020). MAST: Model-based Analysis of Single Cell Transcriptomics. R package version 1.16.0, https://github.com/RGLab/MAST/.
#'
#' . Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019)
#'
#' . Seurat 3 source at https://github.com/satijalab/seurat/blob/master/R/differential_expression.R

###-------------test--------------------------###
object = so.sub
random_effect.vars = NULL
ident.1 = paste(celltype_3rd, "ALS", sep = "_")
ident.2 = paste(celltype_3rd, "CTRL", sep = "_")
group.by = "CellType_3rd_Condition"
base = 2
latent.vars = c("Age", "Sex")
n_cores = 12

cells.1 = NULL
cells.2 = NULL
logfc.threshold = 0.25
assay=NULL
slot="data"
features = NULL
min.pct = 0.1
max.cells.per.ident = NULL
random.seed = 1
verbose=TRUE
p.adjust.method="bonferroni" # seuratのFindMarkersのデフォルト
###-------------test--------------------------###


# DE_MAST_RE_seurat = function(
#   object,
#   random_effect.vars,
#   ident.1 = NULL,
#   ident.2 = NULL,
#   cells.1 = NULL,
#   cells.2 = NULL,
#   group.by = NULL,
#   logfc.threshold = 0.25,
#   base = exp(1),
#   assay=NULL,
#   slot="data",
#   features = NULL,
#   min.pct = 0.1,
#   max.cells.per.ident = NULL,
#   random.seed = 1,
#   latent.vars = NULL,
#   n_cores=NULL,
#   verbose=TRUE,
#   p.adjust.method="fdr",
#   ...
# ) {

  # require(Seurat)
  # require(MAST)

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

  # stopifnot(!is.null(random_effect.vars))
  # stopifnot(all(random_effect.vars %in% colnames(object@meta.data)))

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
      data=densemat,
      cells.1=cells.1,
      cells.2=cells.2,
      mean.fxn=mean.fxn,
      fc.name=fc.name,
      features = features)

  # ここでlog2fc > 0.25の遺伝子だけ論理型を逆転
  if (!is.null(logfc.threshold)) {
    vec_logical_features = vec_logical_features & fc.results[[fc.name]]>logfc.threshold  # need to convert back to non-log space to take mean
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

  # add random vars to column data
  for (random_effect.var in random_effect.vars) {
    df_coldata[[random_effect.var]] <- as.factor(
        object@meta.data[c(idx.cells.1,idx.cells.2),random_effect.var]
      )
  }

  # check for NAs in column data
  for (colname in colnames(df_coldata)) {
    if (any(is.na(df_coldata[[colname]]))) {
      stop(paste0("variable '", colname, "' contains NAs"))
    }
  }

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

  #======== prep  feature data =========================

  df_feature = data.frame("primerid"=rownames(densemat)) # primerid is hardcoded feature name in MAST

  #======== prep SingleCellAssay object ============================

  expr = quote(MAST::FromMatrix(exprsArray=densemat,
                                 check_sanity = TRUE,
                                 cData=df_coldata,
                                 fData=df_feature))
  sca <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  # set group as a factor and make group1 the reference
  #cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  # in Seurat source this is Group1!?
  # see https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
  SummarizedExperiment::colData(sca)$group <- relevel(SummarizedExperiment::colData(sca)$group, ref="Group2")
  # cond <- relevel(x = cond, ref = "Group2")
  # SummarizedExperiment::colData(sca)$condition <- cond

  # add the centered number of expressed genes to the model
  # cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
  # SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)

  #======== run MAST test ======================

  str_plus = if (!is.null(latent.vars)) " + " else ""

  fmla <- as.formula(
    object = paste0(" ~ group", str_plus, paste(latent.vars, collapse = " + "))
  )
    
  print(fmla)

  # METHOD = bayesglm
  expr = quote(MAST::zlm(formula = fmla,
                   sca = sca,
                   method = 'bayesglm'))

  zlmCond <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  if (verbose) {
    print(zlmCond)
    message("Compute likelihoods and p-values")
  }
  # call a likelihood ratio test on the fitted object
  expr = quote(MAST::summary(object = zlmCond, doLRT = 'groupGroup1'))
  summaryCond <- if (verbose) eval(expr) else suppressWarnings(suppressMessages(eval(expr)))

  if ("character" %in% class(summaryCond)) stop("No differentially genes detected")
  # print(summaryCond,n=4)
  # Fitted zlm with top 4 genes per contrast:
  #   ( log fold change Z-score )
  # primerid           conditionGroup1 cngeneson
  # Akr1c20              -4.1            13.2*
  # Bclaf3               18.5*            4.8
  # Cmss1                -0.3            -9.9*
  # Cux2                -16.5*            4.6
  # ENSMUSG00000115022  -19.5*            7.6
  # Pot1b                13.5*            5.3
  # Siah1a                4.6            17.6*
  # Sipa1l2               1.8            10.5*
  #
  #
  # summaryCond$datatable has colnames
  # "primerid" feature name
  # "component" if no latent.vars, "C" (continuous) "D" (discrete)  "H" (hurdle) "S" (RE?) "logFC"
  # "contrast" if no latent.vars,  conditionGroup1 (Intercept) cngeneson
  # "ci.hi" confidence interval (high)
  # "ci.lo" confidence interval (low)
  # "coef" coefficient estimate
  # "z" z-scores

  # uses
  # https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#differential-expression-using-a-hurdle-model
  # and
  # https://github.com/kdzimm/PseudoreplicationPaper/blob/master/Type_1_Error/Type%201%20-%20MAST%20RE.Rmd

  de.results <- data.frame(
    "p_val" = summaryCond$datatable[contrast=='groupGroup1' & component=='H', `Pr(>Chisq)`],
    fc.results[vec_logical_features,])  #setDF(summaryCond$datatable[contrast=='groupGroup1' & component=='logFC', .(coef)])

  de.results$p_val_adj = p.adjust(de.results$p_val, method=p.adjust.method, n=length(vec_logical_features))

  de.results = de.results[order(de.results$p_val, -de.results[[fc.name]]),]

  head(de.results)

  test <- de.results[(de.results$avg_log2FC >= 0.25 & de.results$p_val_adj < 0.2), ]
  # head(test, n = 20)
  # nrow(test)

  write.csv(
    test,
    file = file.path(outDir, "GLMLatent", paste0(celltype_3rd, "_glm_padj02_log2fc025_Up.csv")),
    quote = FALSE
    )
  }
}


#' 
#' ## 2-2. Downregurated genes in case
## -----------------------------------------------------------------------------
sub_list <- c("neuron", "others", "astro", "micro", "oligos")

for (i in seq_along(sub_list)){
    # i = 1
    celltype <- sub_list[[i]]

    so.sub <- readRDS(
        file = file.path(
            mainDir,
            "subset",
            paste0("sosub_", toString(celltype), "_3rd_layer_prepSCT.rds")
            )
        )
    so.sub <- subset(so.sub, subset = Project == "RN01", invert = TRUE)


    celltype3rd_list <- unique(so.sub$CellType_3rd)

        for (j in seq_along(celltype3rd_list)){ 
        # j = 3
        celltype_3rd <- celltype3rd_list[[j]]
        print(celltype_3rd)

        DefaultAssay(so.sub) <- 'RNA'
        so.sub <- NormalizeData(so.sub)


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
#  # compute average log fold change（seuratのFindMarkersに入っている関数　expm1($num)はexp($num)-1）
#   pseudocount.use = 1
#   mean.fxn = function(x) {
#     return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
#   }

#' MAST differential expression test with random effect for Seurat object
#'
#' Uses Seurat::FindMarkers MAST test hurdle model with added random effects variables
#'
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
#' @references Zimmerman, K.D., Espeland, M.A. & Langefeld, C.D.
#' . A practical solution to pseudoreplication bias in single-cell studies.
#' . Nat Commun 12, 738 (2021). https://doi.org/10.1038/s41467-021-21038-1
#'
#' . McDavid A, Finak G, Yajima M (2020). MAST: Model-based Analysis of Single Cell Transcriptomics. R package version 1.16.0, https://github.com/RGLab/MAST/.
#'
#' . Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019)
#'
#' . Seurat 3 source at https://github.com/satijalab/seurat/blob/master/R/differential_expression.R

###-------------test--------------------------###
object = so.sub
random_effect.vars = NULL
ident.1 = paste(celltype_3rd, "CTRL", sep = "_")
ident.2 = paste(celltype_3rd, "ALS", sep = "_")
group.by = "CellType_3rd_Condition"
base = 2
latent.vars = c("Age", "Sex")
n_cores = 12

cells.1 = NULL
cells.2 = NULL
logfc.threshold = 0.25
assay=NULL
slot="data"
features = NULL
min.pct = 0.1
max.cells.per.ident = NULL
random.seed = 1
verbose=TRUE
p.adjust.method="bonferroni" # seuratのFindMarkersのデフォルト
###-------------test--------------------------###


# DE_MAST_RE_seurat = function(
#   object,
#   random_effect.vars,
#   ident.1 = NULL,
#   ident.2 = NULL,
#   cells.1 = NULL,
#   cells.2 = NULL,
#   group.by = NULL,
#   logfc.threshold = 0.25,
#   base = exp(1),
#   assay=NULL,
#   slot="data",
#   features = NULL,
#   min.pct = 0.1,
#   max.cells.per.ident = NULL,
#   random.seed = 1,
#   latent.vars = NULL,
#   n_cores=NULL,
#   verbose=TRUE,
#   p.adjust.method="fdr",
#   ...
# ) {

  # require(Seurat)
  # require(MAST)

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

  # stopifnot(!is.null(random_effect.vars))
  # stopifnot(all(random_effect.vars %in% colnames(object@meta.data)))

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
      data=densemat,
      cells.1=cells.1,
      cells.2=cells.2,
      mean.fxn=mean.fxn,
      fc.name=fc.name,
      features = features)

  # ここでlog2fc > 0.25の遺伝子だけ論理型を逆転
  if (!is.null(logfc.threshold)) {
    vec_logical_features = vec_logical_features & fc.results[[fc.name]]>logfc.threshold  # need to convert back to non-log space to take mean
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

  # add random vars to column data
  for (random_effect.var in random_effect.vars) {
    df_coldata[[random_effect.var]] <- as.factor(
        object@meta.data[c(idx.cells.1,idx.cells.2),random_effect.var]
      )
  }

  # check for NAs in column data
  for (colname in colnames(df_coldata)) {
    if (any(is.na(df_coldata[[colname]]))) {
      stop(paste0("variable '", colname, "' contains NAs"))
    }
  }

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

  #======== prep  feature data =========================

  df_feature = data.frame("primerid"=rownames(densemat)) # primerid is hardcoded feature name in MAST

  #======== prep SingleCellAssay object ============================

  expr = quote(MAST::FromMatrix(exprsArray=densemat,
                                 check_sanity = TRUE,
                                 cData=df_coldata,
                                 fData=df_feature))
  sca <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  # set group as a factor and make group1 the reference
  #cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  # in Seurat source this is Group1!?
  # see https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
  SummarizedExperiment::colData(sca)$group <- relevel(SummarizedExperiment::colData(sca)$group, ref="Group2")
  # cond <- relevel(x = cond, ref = "Group2")
  # SummarizedExperiment::colData(sca)$condition <- cond

  # add the centered number of expressed genes to the model
  # cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
  # SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)

  #======== run MAST test ======================

  str_plus = if (!is.null(latent.vars)) " + " else ""

  fmla <- as.formula(
    object = paste0(" ~ group", str_plus, paste(latent.vars, collapse = " + "))
  )
    
  print(fmla)

  # METHOD = bayesglm
  expr = quote(MAST::zlm(formula = fmla,
                   sca = sca,
                   method = 'bayesglm'))

  zlmCond <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  if (verbose) {
    print(zlmCond)
    message("Compute likelihoods and p-values")
  }
  # call a likelihood ratio test on the fitted object
  expr = quote(MAST::summary(object = zlmCond, doLRT = 'groupGroup1'))
  summaryCond <- if (verbose) eval(expr) else suppressWarnings(suppressMessages(eval(expr)))

  if ("character" %in% class(summaryCond)) stop("No differentially genes detected")
  # print(summaryCond,n=4)
  # Fitted zlm with top 4 genes per contrast:
  #   ( log fold change Z-score )
  # primerid           conditionGroup1 cngeneson
  # Akr1c20              -4.1            13.2*
  # Bclaf3               18.5*            4.8
  # Cmss1                -0.3            -9.9*
  # Cux2                -16.5*            4.6
  # ENSMUSG00000115022  -19.5*            7.6
  # Pot1b                13.5*            5.3
  # Siah1a                4.6            17.6*
  # Sipa1l2               1.8            10.5*
  #
  #
  # summaryCond$datatable has colnames
  # "primerid" feature name
  # "component" if no latent.vars, "C" (continuous) "D" (discrete)  "H" (hurdle) "S" (RE?) "logFC"
  # "contrast" if no latent.vars,  conditionGroup1 (Intercept) cngeneson
  # "ci.hi" confidence interval (high)
  # "ci.lo" confidence interval (low)
  # "coef" coefficient estimate
  # "z" z-scores

  # uses
  # https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#differential-expression-using-a-hurdle-model
  # and
  # https://github.com/kdzimm/PseudoreplicationPaper/blob/master/Type_1_Error/Type%201%20-%20MAST%20RE.Rmd

  de.results <- data.frame(
    "p_val" = summaryCond$datatable[contrast=='groupGroup1' & component=='H', `Pr(>Chisq)`],
    fc.results[vec_logical_features,])  #setDF(summaryCond$datatable[contrast=='groupGroup1' & component=='logFC', .(coef)])

  de.results$p_val_adj = p.adjust(de.results$p_val, method=p.adjust.method, n=length(vec_logical_features))

  de.results = de.results[order(de.results$p_val, -de.results[[fc.name]]),]

  head(de.results)

  test <- de.results[(de.results$avg_log2FC >= 0.25 & de.results$p_val_adj < 0.2), ]
  # head(test, n = 20)
  # nrow(test)

  write.csv(
    test,
    file = file.path(outDir, "GLMLatent", paste0(celltype_3rd, "_glm_padj02_log2fc025_Down.csv")),
    quote = FALSE
    )
  }
}


#' 
#' 
#' # 3. MAST without latent variables
#' 
#' ## 3-1. Upregurated genes in case
## -----------------------------------------------------------------------------
sub_list <- c("neuron", "others", "astro", "micro", "oligos")

for (i in seq_along(sub_list)){
    # i = 1
    celltype <- sub_list[[i]]

    so.sub <- readRDS(
        file = file.path(
            mainDir,
            "subset",
            paste0("sosub_", toString(celltype), "_3rd_layer_prepSCT.rds")
            )
        )
    so.sub <- subset(so.sub, subset = Project == "RN01", invert = TRUE)


    celltype3rd_list <- unique(so.sub$CellType_3rd)

        for (j in seq_along(celltype3rd_list)){ 
        # j = 3
        celltype_3rd <- celltype3rd_list[[j]]
        print(celltype_3rd)

        DefaultAssay(so.sub) <- 'RNA'
        so.sub <- NormalizeData(so.sub)


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
#  # compute average log fold change（seuratのFindMarkersに入っている関数　expm1($num)はexp($num)-1）
#   pseudocount.use = 1
#   mean.fxn = function(x) {
#     return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
#   }

#' MAST differential expression test with random effect for Seurat object
#'
#' Uses Seurat::FindMarkers MAST test hurdle model with added random effects variables
#'
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
#' @references Zimmerman, K.D., Espeland, M.A. & Langefeld, C.D.
#' . A practical solution to pseudoreplication bias in single-cell studies.
#' . Nat Commun 12, 738 (2021). https://doi.org/10.1038/s41467-021-21038-1
#'
#' . McDavid A, Finak G, Yajima M (2020). MAST: Model-based Analysis of Single Cell Transcriptomics. R package version 1.16.0, https://github.com/RGLab/MAST/.
#'
#' . Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019)
#'
#' . Seurat 3 source at https://github.com/satijalab/seurat/blob/master/R/differential_expression.R

###-------------test--------------------------###
object = so.sub
random_effect.vars = NULL
ident.1 = paste(celltype_3rd, "ALS", sep = "_")
ident.2 = paste(celltype_3rd, "CTRL", sep = "_")
group.by = "CellType_3rd_Condition"
base = 2
latent.vars = NULL
n_cores = 12

cells.1 = NULL
cells.2 = NULL
logfc.threshold = 0.25
assay=NULL
slot="data"
features = NULL
min.pct = 0.1
max.cells.per.ident = NULL
random.seed = 1
verbose=TRUE
p.adjust.method="bonferroni" # seuratのFindMarkersのデフォルト
###-------------test--------------------------###


# DE_MAST_RE_seurat = function(
#   object,
#   random_effect.vars,
#   ident.1 = NULL,
#   ident.2 = NULL,
#   cells.1 = NULL,
#   cells.2 = NULL,
#   group.by = NULL,
#   logfc.threshold = 0.25,
#   base = exp(1),
#   assay=NULL,
#   slot="data",
#   features = NULL,
#   min.pct = 0.1,
#   max.cells.per.ident = NULL,
#   random.seed = 1,
#   latent.vars = NULL,
#   n_cores=NULL,
#   verbose=TRUE,
#   p.adjust.method="fdr",
#   ...
# ) {

  # require(Seurat)
  # require(MAST)

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

  # stopifnot(!is.null(random_effect.vars))
  # stopifnot(all(random_effect.vars %in% colnames(object@meta.data)))

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
      data=densemat,
      cells.1=cells.1,
      cells.2=cells.2,
      mean.fxn=mean.fxn,
      fc.name=fc.name,
      features = features)

  # ここでlog2fc > 0.25の遺伝子だけ論理型を逆転
  if (!is.null(logfc.threshold)) {
    vec_logical_features = vec_logical_features & fc.results[[fc.name]]>logfc.threshold  # need to convert back to non-log space to take mean
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

  # add random vars to column data
  for (random_effect.var in random_effect.vars) {
    df_coldata[[random_effect.var]] <- as.factor(
        object@meta.data[c(idx.cells.1,idx.cells.2),random_effect.var]
      )
  }

  # check for NAs in column data
  for (colname in colnames(df_coldata)) {
    if (any(is.na(df_coldata[[colname]]))) {
      stop(paste0("variable '", colname, "' contains NAs"))
    }
  }

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

  #======== prep  feature data =========================

  df_feature = data.frame("primerid"=rownames(densemat)) # primerid is hardcoded feature name in MAST

  #======== prep SingleCellAssay object ============================

  expr = quote(MAST::FromMatrix(exprsArray=densemat,
                                 check_sanity = TRUE,
                                 cData=df_coldata,
                                 fData=df_feature))
  sca <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  # set group as a factor and make group1 the reference
  #cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  # in Seurat source this is Group1!?
  # see https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
  SummarizedExperiment::colData(sca)$group <- relevel(SummarizedExperiment::colData(sca)$group, ref="Group2")
  # cond <- relevel(x = cond, ref = "Group2")
  # SummarizedExperiment::colData(sca)$condition <- cond

  # add the centered number of expressed genes to the model
  # cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
  # SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)

  #======== run MAST test ======================

  str_plus = if (!is.null(latent.vars)) " + " else ""

  fmla <- as.formula(
    object = paste0(" ~ group", str_plus, paste(latent.vars, collapse = " + "))
  )
    
  print(fmla)

  # METHOD = bayesglm
  expr = quote(MAST::zlm(formula = fmla,
                   sca = sca,
                   method = 'bayesglm'))

  zlmCond <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  if (verbose) {
    print(zlmCond)
    message("Compute likelihoods and p-values")
  }
  # call a likelihood ratio test on the fitted object
  expr = quote(MAST::summary(object = zlmCond, doLRT = 'groupGroup1'))
  summaryCond <- if (verbose) eval(expr) else suppressWarnings(suppressMessages(eval(expr)))

  if ("character" %in% class(summaryCond)) stop("No differentially genes detected")
  # print(summaryCond,n=4)
  # Fitted zlm with top 4 genes per contrast:
  #   ( log fold change Z-score )
  # primerid           conditionGroup1 cngeneson
  # Akr1c20              -4.1            13.2*
  # Bclaf3               18.5*            4.8
  # Cmss1                -0.3            -9.9*
  # Cux2                -16.5*            4.6
  # ENSMUSG00000115022  -19.5*            7.6
  # Pot1b                13.5*            5.3
  # Siah1a                4.6            17.6*
  # Sipa1l2               1.8            10.5*
  #
  #
  # summaryCond$datatable has colnames
  # "primerid" feature name
  # "component" if no latent.vars, "C" (continuous) "D" (discrete)  "H" (hurdle) "S" (RE?) "logFC"
  # "contrast" if no latent.vars,  conditionGroup1 (Intercept) cngeneson
  # "ci.hi" confidence interval (high)
  # "ci.lo" confidence interval (low)
  # "coef" coefficient estimate
  # "z" z-scores

  # uses
  # https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#differential-expression-using-a-hurdle-model
  # and
  # https://github.com/kdzimm/PseudoreplicationPaper/blob/master/Type_1_Error/Type%201%20-%20MAST%20RE.Rmd

  de.results <- data.frame(
    "p_val" = summaryCond$datatable[contrast=='groupGroup1' & component=='H', `Pr(>Chisq)`],
    fc.results[vec_logical_features,])  #setDF(summaryCond$datatable[contrast=='groupGroup1' & component=='logFC', .(coef)])

  de.results$p_val_adj = p.adjust(de.results$p_val, method=p.adjust.method, n=length(vec_logical_features))

  de.results = de.results[order(de.results$p_val, -de.results[[fc.name]]),]

  head(de.results)

  test <- de.results[(de.results$avg_log2FC >= 0.25 & de.results$p_val_adj < 0.2), ]
  # head(test, n = 20)
  # nrow(test)

  write.csv(
    test,
    file = file.path(outDir, "GLMNoLatent", paste0(celltype_3rd, "_glm_padj02_log2fc025_Up.csv")),
    quote = FALSE
    )
  }
}


#' 
#' ## 3-2. Downregurated genes in case
## -----------------------------------------------------------------------------
sub_list <- c("neuron", "others", "astro", "micro", "oligos")

for (i in seq_along(sub_list)){
    # i = 1
    celltype <- sub_list[[i]]

    so.sub <- readRDS(
        file = file.path(
            mainDir,
            "subset",
            paste0("sosub_", toString(celltype), "_3rd_layer_prepSCT.rds")
            )
        )
    so.sub <- subset(so.sub, subset = Project == "RN01", invert = TRUE)


    celltype3rd_list <- unique(so.sub$CellType_3rd)

        for (j in seq_along(celltype3rd_list)){ 
        # j = 3
        celltype_3rd <- celltype3rd_list[[j]]
        print(celltype_3rd)

        DefaultAssay(so.sub) <- 'RNA'
        so.sub <- NormalizeData(so.sub)


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
#  # compute average log fold change（seuratのFindMarkersに入っている関数　expm1($num)はexp($num)-1）
#   pseudocount.use = 1
#   mean.fxn = function(x) {
#     return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
#   }

#' MAST differential expression test with random effect for Seurat object
#'
#' Uses Seurat::FindMarkers MAST test hurdle model with added random effects variables
#'
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
#' @references Zimmerman, K.D., Espeland, M.A. & Langefeld, C.D.
#' . A practical solution to pseudoreplication bias in single-cell studies.
#' . Nat Commun 12, 738 (2021). https://doi.org/10.1038/s41467-021-21038-1
#'
#' . McDavid A, Finak G, Yajima M (2020). MAST: Model-based Analysis of Single Cell Transcriptomics. R package version 1.16.0, https://github.com/RGLab/MAST/.
#'
#' . Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019)
#'
#' . Seurat 3 source at https://github.com/satijalab/seurat/blob/master/R/differential_expression.R

###-------------test--------------------------###
object = so.sub
random_effect.vars = NULL
ident.1 = paste(celltype_3rd, "CTRL", sep = "_")
ident.2 = paste(celltype_3rd, "ALS", sep = "_")
group.by = "CellType_3rd_Condition"
base = 2
latent.vars = NULL
n_cores = 12

cells.1 = NULL
cells.2 = NULL
logfc.threshold = 0.25
assay=NULL
slot="data"
features = NULL
min.pct = 0.1
max.cells.per.ident = NULL
random.seed = 1
verbose=TRUE
p.adjust.method="bonferroni" # seuratのFindMarkersのデフォルト
###-------------test--------------------------###


# DE_MAST_RE_seurat = function(
#   object,
#   random_effect.vars,
#   ident.1 = NULL,
#   ident.2 = NULL,
#   cells.1 = NULL,
#   cells.2 = NULL,
#   group.by = NULL,
#   logfc.threshold = 0.25,
#   base = exp(1),
#   assay=NULL,
#   slot="data",
#   features = NULL,
#   min.pct = 0.1,
#   max.cells.per.ident = NULL,
#   random.seed = 1,
#   latent.vars = NULL,
#   n_cores=NULL,
#   verbose=TRUE,
#   p.adjust.method="fdr",
#   ...
# ) {

  # require(Seurat)
  # require(MAST)

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

  # stopifnot(!is.null(random_effect.vars))
  # stopifnot(all(random_effect.vars %in% colnames(object@meta.data)))

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
      data=densemat,
      cells.1=cells.1,
      cells.2=cells.2,
      mean.fxn=mean.fxn,
      fc.name=fc.name,
      features = features)

  # ここでlog2fc > 0.25の遺伝子だけ論理型を逆転
  if (!is.null(logfc.threshold)) {
    vec_logical_features = vec_logical_features & fc.results[[fc.name]]>logfc.threshold  # need to convert back to non-log space to take mean
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

  # add random vars to column data
  for (random_effect.var in random_effect.vars) {
    df_coldata[[random_effect.var]] <- as.factor(
        object@meta.data[c(idx.cells.1,idx.cells.2),random_effect.var]
      )
  }

  # check for NAs in column data
  for (colname in colnames(df_coldata)) {
    if (any(is.na(df_coldata[[colname]]))) {
      stop(paste0("variable '", colname, "' contains NAs"))
    }
  }

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

  #======== prep  feature data =========================

  df_feature = data.frame("primerid"=rownames(densemat)) # primerid is hardcoded feature name in MAST

  #======== prep SingleCellAssay object ============================

  expr = quote(MAST::FromMatrix(exprsArray=densemat,
                                 check_sanity = TRUE,
                                 cData=df_coldata,
                                 fData=df_feature))
  sca <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  # set group as a factor and make group1 the reference
  #cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  # in Seurat source this is Group1!?
  # see https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
  SummarizedExperiment::colData(sca)$group <- relevel(SummarizedExperiment::colData(sca)$group, ref="Group2")
  # cond <- relevel(x = cond, ref = "Group2")
  # SummarizedExperiment::colData(sca)$condition <- cond

  # add the centered number of expressed genes to the model
  # cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
  # SummarizedExperiment::colData(sca)$cngeneson <- scale(cdr2)

  #======== run MAST test ======================

  str_plus = if (!is.null(latent.vars)) " + " else ""

  fmla <- as.formula(
    object = paste0(" ~ group", str_plus, paste(latent.vars, collapse = " + "))
  )
    
  print(fmla)

  # METHOD = bayesglm
  expr = quote(MAST::zlm(formula = fmla,
                   sca = sca,
                   method = 'bayesglm'))

  zlmCond <- if (verbose) eval(expr) else suppressMessages(eval(expr))

  if (verbose) {
    print(zlmCond)
    message("Compute likelihoods and p-values")
  }
  # call a likelihood ratio test on the fitted object
  expr = quote(MAST::summary(object = zlmCond, doLRT = 'groupGroup1'))
  summaryCond <- if (verbose) eval(expr) else suppressWarnings(suppressMessages(eval(expr)))

  if ("character" %in% class(summaryCond)) stop("No differentially genes detected")
  # print(summaryCond,n=4)
  # Fitted zlm with top 4 genes per contrast:
  #   ( log fold change Z-score )
  # primerid           conditionGroup1 cngeneson
  # Akr1c20              -4.1            13.2*
  # Bclaf3               18.5*            4.8
  # Cmss1                -0.3            -9.9*
  # Cux2                -16.5*            4.6
  # ENSMUSG00000115022  -19.5*            7.6
  # Pot1b                13.5*            5.3
  # Siah1a                4.6            17.6*
  # Sipa1l2               1.8            10.5*
  #
  #
  # summaryCond$datatable has colnames
  # "primerid" feature name
  # "component" if no latent.vars, "C" (continuous) "D" (discrete)  "H" (hurdle) "S" (RE?) "logFC"
  # "contrast" if no latent.vars,  conditionGroup1 (Intercept) cngeneson
  # "ci.hi" confidence interval (high)
  # "ci.lo" confidence interval (low)
  # "coef" coefficient estimate
  # "z" z-scores

  # uses
  # https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#differential-expression-using-a-hurdle-model
  # and
  # https://github.com/kdzimm/PseudoreplicationPaper/blob/master/Type_1_Error/Type%201%20-%20MAST%20RE.Rmd

  de.results <- data.frame(
    "p_val" = summaryCond$datatable[contrast=='groupGroup1' & component=='H', `Pr(>Chisq)`],
    fc.results[vec_logical_features,])  #setDF(summaryCond$datatable[contrast=='groupGroup1' & component=='logFC', .(coef)])

  de.results$p_val_adj = p.adjust(de.results$p_val, method=p.adjust.method, n=length(vec_logical_features))

  de.results = de.results[order(de.results$p_val, -de.results[[fc.name]]),]

  head(de.results)

  test <- de.results[(de.results$avg_log2FC >= 0.25 & de.results$p_val_adj < 0.2), ]
  # head(test, n = 20)
  # nrow(test)

  write.csv(
    test,
    file = file.path(outDir, "GLMNoLatent", paste0(celltype_3rd, "_glm_padj02_log2fc025_Down.csv")),
    quote = FALSE
    )
  }
}


#' 
#' 
#' 
#' 
#' # 4. Errorに関するテスト
#' 
#' 
#' ## 4-1. Stack Overflowのサンプルデータ
## -----------------------------------------------------------------------------
data4<-structure(list(code = structure(1:10, .Label = c("10888", "10889", 
"10890", "10891", "10892", "10893", "10894", "10896", "10897", 
"10898", "10899", "10900", "10901", "10902", "10903", "10904", 
"10905", "10906", "10907", "10908", "10909", "10910", "10914", 
"10916", "10917", "10919", "10920", "10922", "10923", "10924", 
"10925", "10927"), class = "factor"), speed = c(0.0296315046039244, 
0.0366986630049636, 0.0294297725505692, 0.048316183511095, 0.0294275666501456, 
0.199924957584131, 0.0798850288176711, 0.0445886457047146, 0.0285993712316451, 
0.0715158276875623), meanflow = c(0.657410742496051, 0.608271363339857, 
0.663241108786611, 0.538259450171821, 0.666299529534762, 0.507156583629893, 
0.762448863636364, 37.6559178370787, 50.8557196935557, 31.6601587837838
), length = c(136, 157, 132, 140, 135, 134, 144, 149, 139, 165
), river = structure(c(2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L
), .Label = c("c", "f"), class = "factor")), .Names = c("code", 
"speed", "meanflow", "length", "river"), row.names = c(2L, 4L, 
6L, 8L, 10L, 12L, 14L, 16L, 18L, 20L), class = "data.frame")

data4

model1<-lmer(speed ~ river + length +(1|meanflow)+(1|code), data4)


model2 <- lmer(speed ~ river + length +(1|meanflow)+(1|code), data4,
    control=lmerControl(check.nobs.vs.nlev = "ignore",
     check.nobs.vs.rankZ = "ignore",
     check.nobs.vs.nRE="ignore"))

#' 
#' # 2-2. cell numbersの内訳
## -----------------------------------------------------------------------------
so.m <- readRDS(
        file = file.path(
            mainDir, "som_spinalcord_3rdannotation_LogNormalized.rds"
            )
    )
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(object = so.m, slot = "counts", assay = "RNA")
metadata <- so.m@meta.data %>%
            dplyr::select(c("Age", "Sex", "Onset", "Condition",
                "CellType_1st", "CellType_3rd", "SampleID"))

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(so.m$CellType_3rd)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts),
                        colData = metadata)

sids <- levels(as.factor(colData(sce)$SampleID))
cids <- levels(as.factor(colData(sce)$cluster_id))

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$SampleID))
n_cells_cluster <- data.frame(table(sce$cluster_id, sce$SampleID))
colnames(n_cells_cluster) <- c("cluster_id", "SampleID", "n_cells")
head(n_cells_cluster)

# ## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
# m <- match(sids, sce$SampleID)
# m_c <- match(cids, sce$cluster_id)
# ## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
# ei <- data.frame(colData(sce)[m, ], 
#                 n_cells, row.names = NULL) %>% 
#                 select(-"cluster_id")
# ei

# ei_c <- data.frame(colData(sce)[m_c],
#                   n_cells_cluster, row.names = NULL)

p <- ggplot(
    data = n_cells_cluster, aes(x = SampleID, y = cluster_id, fill = n_cells)) +
    geom_tile() +
    geom_text(aes(label = n_cells), color = "black", size = 4) +
    scale_fill_gradient(low = "white", high = "red") +
    labs(x = "subject_id", y = "cluster_id") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        text = element_text(size = 20)) + 
    theme_classic() 
p

ggsave(
    file = file.path(outDir, "NcellsHeatmap.pdf"),
    plot = p,
    height = 10, width = 10)
#----------- test
# # ggplotを使用してbarplotを描く
# ggplot(test_df, aes(x = name, y = mpg)) +
#   geom_bar(stat = "identity") +
#   labs(x = "Car Name", y = "Miles Per Gallon (mpg)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# data(mtcars)
# test_df <- mtcars
# head(mtcars)
# test_df$name <- rownames(test_df)
# head(test_df)


#' ```
#' 
#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
#' 
## -----------------------------------------------------------------------------
knitr::purl(file.path(outDir, "../01_Preprocessing_CellRangerCountToSeuratObject.Rmd"), documentation = 2)

