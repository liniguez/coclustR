#' Consensus co-clustering of a single-cell object (full pipeline)
#'
#' 1. [run_clustering_grid()]: Leiden/Louvain over `n_dims` x `k` x `resolutions`.
#' 2. [collapse_metabarcodes()]: cells with identical label vectors are merged.
#' 3. [coclustering_matrix()]: fraction of runs in which two meta-barcodes
#'    share a cluster.
#' 4. [core_stability()] + [flag_unstable()]: small, weakly co-clustered
#'    meta-barcodes are excluded from the tree.
#' 5. [consensus_clusters()]: hclust on `1 - co`, cut with dynamicTreeCut.
#' 6. Labels are mapped back to cells; unassigned/unstable cells are labelled
#'    `"Unstable"`.
#'
#' @inheritParams run_clustering_grid
#' @param stability_threshold Co-clustering threshold that defines the core
#'   of a meta-barcode (original `cluster_stability_thershold`).
#' @param min_core_cells,min_core_stability See [flag_unstable()].
#' @param linkage,methods,min_cluster_size,min_cluster_cells,deep_split,rescue_min_cells
#'   See [consensus_clusters()]. The rescue threshold is `stability_threshold`.
#' @param grid Optional pre-computed [run_clustering_grid()] result, to
#'   re-run only the consensus with different settings.
#' @param prefix Prefix of the metadata columns added to the object.
#' @param store Store the full result in the object
#'   (`metadata(sce)$coclustR` or `Misc(seurat, "coclustR")`).
#' @param store_matrix Keep the co-clustering matrix in the stored result
#'   (needed for [plot_coclustering_heatmap()]; can be large).
#' @return The input object with new metadata columns
#'   `<prefix>Hybrid` / `<prefix>DynamicTree` (factors) and, if `store`, the
#'   full `coclust_result`. If `object` is a matrix, the `coclust_result`
#'   itself is returned.
#' @examples
#' \dontrun{
#' future::plan(future::multisession, workers = 4)
#' sce <- run_coclustering(sce, reduction = "PCA")
#' seu <- run_coclustering(seu, reduction = "pca")
#' res <- get_coclust_result(sce)
#' plot_grid_metrics(res)
#' plot_clusters_embedding(sce, color_by = "coclust_Hybrid")
#' }
#' @export
run_coclustering <- function(object, reduction = NULL,
                             n_dims = c(15, 30, 50, 75, 100),
                             k = seq(10, 40, 5),
                             resolutions = seq(0.1, 1.4, 0.1),
                             algorithm = c("leiden", "louvain"),
                             prune = 1 / 15, n_iter = 2,
                             silhouette = FALSE, seed = 11,
                             stability_threshold = 0.5,
                             min_core_cells = 50, min_core_stability = 0.8,
                             linkage = "complete",
                             methods = c("hybrid", "dynamic"),
                             min_cluster_size = 50, min_cluster_cells = NULL,
                             deep_split = FALSE, rescue_min_cells = 50,
                             grid = NULL, prefix = "coclust_",
                             store = TRUE, store_matrix = TRUE,
                             BNPARAM = NULL, verbose = TRUE) {
  if (is.null(grid)) {
    grid <- run_clustering_grid(
      object, reduction = reduction, n_dims = n_dims, k = k,
      resolutions = resolutions, algorithm = algorithm, prune = prune,
      n_iter = n_iter, silhouette = silhouette, seed = seed,
      BNPARAM = BNPARAM, verbose = verbose)
  }
  res <- coclust_from_grid(
    grid, stability_threshold = stability_threshold,
    min_core_cells = min_core_cells, min_core_stability = min_core_stability,
    linkage = linkage, methods = methods, min_cluster_size = min_cluster_size,
    min_cluster_cells = min_cluster_cells, deep_split = deep_split,
    rescue_min_cells = rescue_min_cells, store_matrix = store_matrix, verbose = verbose)

  if (.object_type(object) == "matrix") return(res)

  cl <- res$cell_clusters
  if (!setequal(rownames(cl), .cell_names(object)))
    stop("Cell names in `grid` do not match the object.", call. = FALSE)
  colnames(cl) <- paste0(prefix, colnames(cl))
  object <- .add_cell_metadata(object, cl)
  if (store) object <- .store_result(object, res)
  object
}

#' Consensus step only, starting from a clustering grid
#'
#' Useful to try different consensus settings without re-running the grid.
#'
#' @param grid Output of [run_clustering_grid()].
#' @inheritParams run_coclustering
#' @return A `coclust_result` list with elements `grid`, `metabarcodes`,
#'   `coclust` (matrix or `NULL`), `stability`, `unstable`, `tree`,
#'   `meta_clusters`, `cell_clusters` and `params`.
#' @export
coclust_from_grid <- function(grid, stability_threshold = 0.5,
                              min_core_cells = 50, min_core_stability = 0.8,
                              linkage = "complete",
                              methods = c("hybrid", "dynamic"),
                              min_cluster_size = 50, min_cluster_cells = NULL,
                              deep_split = FALSE, rescue_min_cells = 50,
                              store_matrix = TRUE,
                              verbose = TRUE) {
  stopifnot(inherits(grid, "coclust_grid"))
  mb <- collapse_metabarcodes(grid$labels)
  n_mb <- nrow(mb$labels)
  if (verbose) message(nrow(grid$labels), " cells collapsed into ", n_mb,
                       " meta-barcodes")
  if (n_mb > 30000)
    warning(n_mb, " meta-barcodes: the co-clustering matrix needs ~",
            round(n_mb^2 * 8 * 3 / 1024^3, 1), " GB of RAM.", call. = FALSE)

  co <- coclustering_matrix(mb$labels)
  stab <- core_stability(co, mb$sizes, threshold = stability_threshold)
  unstable <- flag_unstable(stab, min_core_cells = min_core_cells,
                            min_core_stability = min_core_stability)
  if (verbose) message(sum(unstable), " meta-barcodes (",
                       sum(mb$sizes[unstable]), " cells) flagged as unstable")
  keep <- !unstable
  if (sum(keep) < 2) stop("Fewer than 2 stable meta-barcodes; relax ",
                          "`min_core_cells` / `min_core_stability`.",
                          call. = FALSE)

  cons <- consensus_clusters(co[keep, keep, drop = FALSE], mb$sizes[keep],
                             linkage = linkage, methods = methods,
                             min_cluster_size = min_cluster_size,
                             min_cluster_cells = min_cluster_cells,
                             deep_split = deep_split,
                             rescue_min_cells = rescue_min_cells,
                             rescue_threshold = stability_threshold)

  # Meta-barcode level labels (0 for the unstable ones)
  meta_lab <- matrix(0L, n_mb, ncol(cons$labels),
                     dimnames = list(rownames(mb$labels), colnames(cons$labels)))
  meta_lab[keep, ] <- cons$labels
  meta_clusters <- as.data.frame(lapply(as.data.frame(meta_lab),
                                        .labels_to_factor),
                                 row.names = rownames(meta_lab))

  cell_clusters <- meta_clusters[mb$cell2meta, , drop = FALSE]
  rownames(cell_clusters) <- names(mb$cell2meta)

  structure(
    list(grid = grid, metabarcodes = mb,
         coclust = if (store_matrix) co else NULL,
         stability = stab, unstable = unstable, tree = cons$tree,
         meta_clusters = meta_clusters, cell_clusters = cell_clusters,
         params = list(stability_threshold = stability_threshold,
                       min_core_cells = min_core_cells,
                       min_core_stability = min_core_stability,
                       linkage = linkage, methods = methods,
                       min_cluster_size = min_cluster_size,
                       min_cluster_cells = min_cluster_cells,
                       deep_split = deep_split,
                       rescue_min_cells = rescue_min_cells)),
    class = "coclust_result")
}

#' @export
print.coclust_result <- function(x, ...) {
  cat("<coclust_result>\n")
  cat("  cells:", nrow(x$cell_clusters), " meta-barcodes:",
      nrow(x$metabarcodes$labels), " runs:", ncol(x$grid$labels), "\n")
  cat("  unstable meta-barcodes:", sum(x$unstable), "(",
      sum(x$metabarcodes$sizes[x$unstable]), "cells )\n")
  for (m in colnames(x$cell_clusters)) {
    tab <- table(x$cell_clusters[[m]])
    cat("  ", m, ": ", length(tab) - 1, " clusters + Unstable (",
        tab[["Unstable"]], " cells)\n", sep = "")
  }
  invisible(x)
}
