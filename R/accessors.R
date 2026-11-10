# Internal helpers that hide the differences between SingleCellExperiment,
# Seurat and plain matrices. Every other function in the package goes through
# these, so supporting a new container only requires touching this file.

.object_type <- function(object) {
  if (inherits(object, "SingleCellExperiment")) return("sce")
  if (inherits(object, "Seurat")) return("seurat")
  if (is.matrix(object) || inherits(object, "Matrix")) return("matrix")
  stop("`object` must be a SingleCellExperiment, a Seurat object or a ",
       "cells x dimensions matrix; got <", class(object)[1], ">.", call. = FALSE)
}

.need <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required ", why, ". Please install it.",
         call. = FALSE)
  }
}

.default_reduction <- function(type, what = c("pca", "umap")) {
  what <- match.arg(what)
  switch(type,
         sce = c(pca = "PCA", umap = "UMAP")[[what]],
         seurat = c(pca = "pca", umap = "umap")[[what]],
         matrix = NULL)
}

#' Extract a low-dimensional embedding (cells x dimensions)
#'
#' @param object A `SingleCellExperiment`, a `Seurat` object or a matrix with
#'   cells in rows.
#' @param reduction Name of the reduction. Defaults to `"PCA"` for
#'   SingleCellExperiment and `"pca"` for Seurat. Ignored for matrices.
#' @return A numeric matrix with cells in rows and cell names as rownames.
#' @export
get_embedding <- function(object, reduction = NULL) {
  type <- .object_type(object)
  if (is.null(reduction)) reduction <- .default_reduction(type, "pca")
  emb <- switch(
    type,
    sce = {
      .need("SingleCellExperiment", "for SingleCellExperiment input")
      avail <- SingleCellExperiment::reducedDimNames(object)
      if (!reduction %in% avail)
        stop("Reduction '", reduction, "' not found. Available: ",
             paste(avail, collapse = ", "), call. = FALSE)
      SingleCellExperiment::reducedDim(object, reduction)
    },
    seurat = {
      .need("SeuratObject", "for Seurat input")
      avail <- SeuratObject::Reductions(object)
      if (!reduction %in% avail)
        stop("Reduction '", reduction, "' not found. Available: ",
             paste(avail, collapse = ", "), call. = FALSE)
      SeuratObject::Embeddings(object, reduction = reduction)
    },
    matrix = as.matrix(object)
  )
  emb <- as.matrix(emb)
  if (type != "matrix") rownames(emb) <- .cell_names(object)
  if (is.null(rownames(emb)))
    rownames(emb) <- paste0("cell", seq_len(nrow(emb)))
  if (anyDuplicated(rownames(emb)))
    stop("Cell names must be unique.", call. = FALSE)
  emb
}

.cell_names <- function(object) {
  type <- .object_type(object)
  nm <- switch(type,
               sce = colnames(object),
               seurat = colnames(object),
               matrix = rownames(object))
  if (is.null(nm))
    stop("The object has no cell names (colnames). Please set unique ",
         "cell names before running the co-clustering.", call. = FALSE)
  nm
}

#' Get cell-level metadata as a data.frame
#'
#' @inheritParams get_embedding
#' @return A data.frame with one row per cell (rownames = cell names).
#' @export
get_cell_metadata <- function(object) {
  switch(.object_type(object),
         sce = {
           .need("SummarizedExperiment", "for SingleCellExperiment input")
           as.data.frame(SummarizedExperiment::colData(object))
         },
         seurat = object[[]],
         matrix = data.frame(row.names = rownames(object)))
}

# Add the columns of `df` (rownames = cells) to the cell metadata.
.add_cell_metadata <- function(object, df) {
  df <- df[.cell_names(object), , drop = FALSE]
  switch(.object_type(object),
         sce = {
           for (cl in colnames(df))
             SummarizedExperiment::colData(object)[[cl]] <- df[[cl]]
           object
         },
         seurat = SeuratObject::AddMetaData(object, metadata = df),
         matrix = object)
}

.store_result <- function(object, result, slot = "coclustR") {
  switch(.object_type(object),
         sce = {
           .need("S4Vectors", "to store the result in a SingleCellExperiment")
           S4Vectors::metadata(object)[[slot]] <- result
           object
         },
         seurat = {
           # Misc<- drops the S3 class, so write to the slot directly
           misc <- methods::slot(object, "misc")
           misc[[slot]] <- result
           methods::slot(object, "misc") <- misc
           object
         },
         matrix = object)
}

#' Retrieve a co-clustering result stored inside an object
#'
#' @inheritParams get_embedding
#' @param slot Name under which [run_coclustering()] stored the result.
#' @return A `coclust_result` object.
#' @export
get_coclust_result <- function(object, slot = "coclustR") {
  if (inherits(object, "coclust_result")) return(object)
  res <- switch(.object_type(object),
                sce = S4Vectors::metadata(object)[[slot]],
                seurat = methods::slot(object, "misc")[[slot]],
                matrix = NULL)
  if (is.null(res))
    stop("No co-clustering result found in slot '", slot, "'.", call. = FALSE)
  res
}
