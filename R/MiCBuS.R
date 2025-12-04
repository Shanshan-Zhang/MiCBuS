#' Identify marker genes from scRNA-seq data for known cell types using Seurat
#'
#' This function identifies marker genes for specified known cell types
#' from a single-cell RNA-seq dataset. It leverages the Seurat framework
#' for normalization, clustering, and differential expression analysis.
#'
#' @param scdat A single-cell expression dataset, typically an object
#'        from the \pkg{Biobase} class (`ExpressionSet`) or any object
#'        containing accessible `exprs(scdat)` and `pData(scdat)` slots.
#' @param project Character string. Name of the Seurat project to create.
#'        Defaults to `"scMarker"`.
#' @param indent Character string. Metadata column name used as cell identity
#'        (e.g., `"cluster"` or `"celltype"`).
#' @param ct Character vector specifying the known cell types (or clusters)
#'        to include in marker detection.
#' @param refdat A numeric matrix or data frame representing the reference
#'        (generated from scdat) expression data for corresponding
#'        cell types. Rows correspond to genes, columns to cell types.
#'
#' @return A data frame containing marker genes identified by Seurat’s
#'         `FindAllMarkers()` along with FC_exp_score.
#'
#' @import Seurat
#' @export
scMarkerSeurat <- function(scdat, project = "scMarker", indent, ct, refdat) {
  # for R CMD check: avoid "no visible binding for global variable 'cluster'"
  cluster <- NULL

  # Initialize the Seurat object with the raw (non-normalized data)
  scdat.seu <- CreateSeuratObject(
    counts = Biobase::exprs(scdat),
    project = project,
    min.cells = 3,
    min.features = 200
  )
  scdat.seu <- AddMetaData(scdat.seu, Biobase::pData(scdat))
  scdat.seu <- NormalizeData(scdat.seu)
  scdat.seu <- SetIdent(scdat.seu, value = indent)
  scdat.seu <- subset(scdat.seu, cluster %in% ct)

  scDEGs <- FindAllMarkers(
    scdat.seu,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

  # Normalize reference dataset and map gene-level expression
  refdat.norm <- SECRET::normcount(refdat, "cpm")
  expression_values <- sapply(seq_len(nrow(scDEGs)), function(i) {
    gene <- scDEGs$gene[i]
    cluster <- scDEGs$cluster[i]
    refdat[gene, cluster]
  })

  # Add metrics to results
  scDEGs$refct_exp <- expression_values
  median50 <- stats::median(scDEGs$refct_exp) * 50
  scDEGs$refct_exp_adj <- pmin(scDEGs$refct_exp, median50)
  scDEGs$FC_exp_score <- (scDEGs$avg_log2FC)^2 * log10(scDEGs$refct_exp_adj + 1)

  return(scDEGs)
}




#' Reference-based deconvolution with SECRET used for p_initial()
#'
#' @param bulkdat A genes-by-samples bulk RNA-seq expression matrix (rownames = genes).
#' @param refdat A genes-by-celltypes reference expression matrix (rownames = genes, columns = cell types).
#' @param scDEGs A data frame containing at least a \code{"gene"} column of marker genes for the known cell types.
#' @param ct A character vector of cell-type names (must match columns of \code{refdat}) to include.
#'
#' @return The first element of \code{SECRET()} output (a matrix of estimated
#'   cell-type proportions).
#'
#' @export
secret_v2c <- function(bulkdat, refdat, scDEGs, ct) {
  # directly load output of SECRET generated reference
  refdat <- refdat[, ct]
  bulkdat <- SECRET::normcount(bulkdat, "cpm")
  refdat <- SECRET::normcount(refdat, "cpm")

  genelist <- scDEGs[, "gene", drop = FALSE]

  # count for gene frequency
  w <- as.data.frame(table(genelist$gene))
  colnames(w)[1] <- "gene"
  rownames(w) <- w$gene

  ### find the shared genes and subset the data
  glist <- Reduce(
    intersect,
    list(rownames(bulkdat), rownames(refdat), genelist$gene)
  )
  w <- w[glist, ]
  bulkdat2 <- bulkdat[rownames(bulkdat) %in% w$gene, ]
  refdat2  <- refdat[rownames(refdat) %in% w$gene, ]

  w$weight <- 1 / w$Freq
  wt <- w$weight

  est_prop <- SECRET::SECRET(
    bulkdat2,
    refdat2,
    withUnknown = FALSE,
    w = wt,
    yNorm = "cpm",
    bNorm = "cpm"
  )[[1]]

  return(est_prop)
}



#' Estimate initial cell-type proportions from bulk and single-cell data
#'
#' This function estimates initial cell-type proportions in bulk RNA-seq data
#' using reference-based deconvolution methods SECRET without unknown cell types.
#' It calculates mean across estimated proportions for each cell type.
#'
#' @param bulkdat A bulk RNA-seq expression matrix or `ExpressionSet` object.
#' @param refdat A single-cell RNA-seq reference matrix or `ExpressionSet` object.
#' @param scDEGs A data frame of marker genes (typically from single-cell DE analysis).
#' @param ct A vector of cell type names to include in the estimation.
#'
#' @return A data frame with estimated cell-type proportions (rows = cell types),
#'   including the mean and standard deviation across samples.
#'
#' @export
p_initial <- function(bulkdat, refdat, scDEGs, ct) {
  est_prop <- secret_v2c(bulkdat, refdat, scDEGs, ct)
  est_prop.t <- as.data.frame(t(est_prop))
  est_prop.t$mean <- apply(est_prop.t, 1, mean)
  est_prop.t$sd   <- apply(est_prop.t, 1, stats::sd)
  est_prop.t$cellType <- rownames(est_prop.t)
  return(est_prop.t)
}



#' Generate Dirichlet-based pseudo-bulk RNA-seq samples
#'
#' @param nsamples Number of Dirichlet-based pseudo-bulk RNA-seq samples to generate
#' @param prop Initial proportions of cell types
#' @param s A parameter to adjust the spread around initial proportions
#' @param ct Cell types existed in the scRNA-seq data
#' @param ds An object prepared for SimBu::simulate_bulk()
#' @param ncells Number of cells for each sample
#'
#' @return A list
#' @export
simulate_bulk_dir <- function(nsamples, prop, s = 10, ct, ds, ncells = 1000) {
  pdir <- rdirichlet(nsamples, alpha = prop * s)

  pure_scenario_dataframe <- data.frame(
    matrix(
      pdir,
      ncol = length(ct),
      dimnames = list(paste0("dirichlet", seq_len(nsamples)), ct)
    )
  )

  simulation.rd <- SimBu::simulate_bulk(
    data = ds,
    scenario = "custom",
    custom_scenario_data = pure_scenario_dataframe,
    scaling_factor = "NONE",
    ncells = ncells,
    nsamples = nsamples,
    BPPARAM = BiocParallel::MulticoreParam(workers = 4),
    run_parallel = TRUE
  )
  return(simulation.rd)
}



#' DESeq2 DE with LFC shrinkage for pseudo-dirichlet-bulk vs. bulk
#'
#' Runs a simple two-group differential expression using DESeq2 on a matrix
#' of read counts where the first \code{n} columns are baseline bulk samples
#' (labeled \code{"bulk"}) and the next \code{m} columns are Dirichlet-based
#' pseudo-bulk samples (labeled \code{"dir_bulk"}). Returns a \code{DESeqResults}
#' object with shrunken log2 fold changes.
#'
#' @param readcounts A genes-by-samples count matrix (or data frame coercible to
#'   matrix) of non-negative integer read counts. Columns are samples.
#' @param n Integer. Number of baseline bulk samples at the beginning of
#'   \code{readcounts}.
#' @param m Integer. Number of Dirichlet-based pseudo-bulk samples following the
#'   baseline samples in \code{readcounts}.
#'
#' @details
#' The design is \code{~ condition} with factor levels ordered as
#' \code{c("bulk", "dir_bulk")}, so the reported log2 fold change is
#' \code{bulk - dir_bulk} in the call to \code{lfcShrink(contrast = c("condition","bulk","dir_bulk"))}.
#'
#' @return A \code{DESeqResults} object with shrunken log2 fold-changes.
#'
#' @export
deseqShrunk_counts <- function(readcounts, n, m) {
  # Basic checks
  stopifnot(is.numeric(n), is.numeric(m), n >= 1, m >= 1)
  readcounts <- as.matrix(readcounts)
  if (ncol(readcounts) < (n + m)) {
    stop("ncol(readcounts) must be at least n + m.")
  }
  if (is.null(colnames(readcounts))) {
    colnames(readcounts) <- paste0("s", seq_len(ncol(readcounts)))
  }

  # Label the pseudo-bulk columns
  colnames(readcounts)[(n + 1):(n + m)] <- paste0("dir_bulk", seq_len(m))

  # Build sample metadata with explicit level order ("bulk" reference)
  sample_info <- data.frame(
    row.names = colnames(readcounts),
    condition = factor(
      c(rep("bulk", n), rep("dir_bulk", m)),
      levels = c("bulk", "dir_bulk")
    )
  )

  # DESeq2 pipeline
  DESeq.ds <- DESeq2::DESeqDataSetFromMatrix(
    countData = readcounts,
    colData   = sample_info,
    design    = ~ condition
  )
  DESeq.ds <- DESeq2::estimateSizeFactors(DESeq.ds)
  DESeq.ds <- DESeq2::DESeq(DESeq.ds)

  DGE.results <- DESeq2::results(DESeq.ds, alpha = 0.05)

  # Shrink LFCs (type = "normal" is broadly available; consider "apeglm" if installed)
  res.shrunk <- DESeq2::lfcShrink(
    dds      = DESeq.ds,
    contrast = c("condition", "bulk", "dir_bulk"),
    res      = DGE.results,
    type     = "normal"
  )

  return(res.shrunk)
}



#' MiCBuS: derive marker gene for unknown cell types using bulk and Dirichlet-pseudo bulk data.
#'
#' This function identifies pseudo marker genes for unknown cell types by
#' combining bulk RNA-seq data with Dirichlet-based pseudo-bulk samples
#' generated from single-cell RNA-seq.
#'
#' (1) validating bulk and single-cell inputs;
#' (2) constructing a cell-type reference;
#' (3) identifying scRNA-seq marker genes using \code{scMarkerSeurat()};
#' (4) estimating initial bulk cell-type proportions;
#' (5) preparing a \pkg{SimBu} dataset;
#' (6) generating Dirichlet-based pseudo-bulk RNA-seq samples;
#' (7) performing DESeq2 shrinkage-based differential expression; and
#' (8) extracting significant pseudo marker genes for unknown cell types.
#'
#' @param bulkdat A genes-by-samples integer count matrix for bulk RNA-seq.
#'   Row names must be gene IDs.
#' @param sc_eset A single-cell \code{ExpressionSet} with
#'   \code{exprs(sc_eset)} as counts and metadata variables
#'   \code{cluster} and \code{sample}.
#' @param m Integer. Number of Dirichlet-based pseudo-bulk samples.
#' @param seed Optional integer for reproducibility.
#'
#' @return A data frame containing the pseudo-bulk DE results (DEpseudo.df).
#'
#' @export
#' @import alabama
#' @examples
#' \dontrun{
#' DEpseudo.df <- MiCBuS(
#'   bulkdat = bulk_counts,
#'   sc_eset = sc_eset,
#'   m = 20,
#'   seed = 1234
#' )
#' }
MiCBuS <- function(bulkdat, sc_eset, m = 20, seed = NULL) {

  ## ----------------------
  ## 1. INPUT VALIDATION
  ## ----------------------

  ### bulkdat checks
  if (is.null(bulkdat) || !is.matrix(bulkdat))
    stop("bulkdat must be a genes-by-samples matrix.")

  if (is.null(rownames(bulkdat)))
    stop("bulkdat must have gene IDs as rownames.")

  if (!is.numeric(bulkdat))
    stop("bulkdat must contain numeric/integer counts.")

  ### sc_eset checks
  if (!inherits(sc_eset, "ExpressionSet"))
    stop("sc_eset must be an ExpressionSet.")

  pd <- Biobase::pData(sc_eset)
  if (!all(c("cluster", "sample") %in% colnames(pd)))
    stop("pData(sc_eset) must include 'cluster' and 'sample' columns.")

  if (is.null(Biobase::exprs(sc_eset)))
    stop("exprs(sc_eset) must contain count data.")

  ### seed control
  if (!is.null(seed)) set.seed(seed)


  ## ----------------------
  ## 2. SECRET reference
  ## ----------------------
  refdat <- SECRET::scRefer(
    input_data = sc_eset,
    ct_var = "cluster",
    sample_var = "sample"
  )

  ## ----------------------
  ## 3. scRNA-seq marker genes
  ## ----------------------
  ct <- unique(pd$cluster)

  scMarkers <- scMarkerSeurat(
    scdat   = sc_eset,
    project = "sc_eset",
    indent  = "cluster",
    ct      = ct,
    refdat  = refdat
  )

  ## ----------------------
  ## 4. Initial proportion estimation
  ## ----------------------
  est_prop.t <- p_initial(
    bulkdat = bulkdat,
    refdat  = refdat,
    scDEGs  = scMarkers,
    ct      = ct
  )

  prop1 <- est_prop.t$mean
  names(prop1) <- ct
  n <- ncol(bulkdat)

  ## ----------------------
  ## 5. Prepare SimBu dataset
  ## ----------------------
  counts <- Matrix::Matrix(Biobase::exprs(sc_eset), sparse = TRUE)
  tpm <- Matrix::t(1e6 * Matrix::t(counts) / Matrix::colSums(counts))

  annotation <- pd
  annotation$ID <- rownames(annotation)
  annotation$cell_type <- annotation$cluster

  ds2 <- SimBu::dataset(
    annotation   = annotation,
    count_matrix = counts,
    tpm_matrix   = tpm,
    name         = "MiCBuS_dataset"
  )

  ## ----------------------
  ## 6. Dirichlet-pseudo-bulk simulation
  ## ----------------------
  simulation.rd2 <- simulate_bulk_dir(
    nsamples = m,
    prop     = prop1,
    s        = 10,
    ct       = ct,
    ds       = ds2,
    ncells   = 1000
  )

  ## ----------------------
  ## 7. DESeq2 shrinkage
  ## ----------------------
  bulkdat.count <- bulkdat
  bulkdat.count.ps <- as.matrix(
    SummarizedExperiment::assays(simulation.rd2$bulk)[["bulk_counts"]]
  )

  common.gene <- intersect(
    rownames(bulkdat.count),
    rownames(bulkdat.count.ps)
  )

  readcounts <- cbind(
    bulkdat.count[common.gene, ],
    bulkdat.count.ps[common.gene, ]
  )

  res.shrunk <- deseqShrunk_counts(readcounts, n = n, m = m)

  ## ----------------------
  ## 8. Extract significant DE genes
  ## ----------------------
  res.shrunk.df <- as.data.frame(res.shrunk)
  DEpseudo.df <- res.shrunk.df[
    res.shrunk.df$padj < 0.05 &
      res.shrunk.df$log2FoldChange > 0,
  ]
  DEpseudo.df$gene <- rownames(DEpseudo.df)

  return(DEpseudo.df)
}
