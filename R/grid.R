#' Build pruned SNN graphs for every combination of PCs and neighbours
#'
#' Mirrors Seurat's `FindNeighbors()`: a shared-nearest-neighbour graph with
#' Jaccard weights, where edges with weight below `prune` are removed
#' (Seurat's `prune.SNN = 1/15`).
#'
#' Note: Seurat's `k.param` counts the cell itself as one of its neighbours,
#' whereas `bluster::makeSNNGraph(k = )` does not. `k = 20` here is therefore
#' roughly equivalent to Seurat's `k.param = 21`.
#'
#' @param emb Cells x dimensions matrix (e.g. from [get_embedding()]).
#' @param n_dims Integer vector with the number of leading dimensions to use.
#'   Values larger than `ncol(emb)` are truncated and duplicates removed.
#' @param k Integer vector with numbers of nearest neighbours.
#' @param prune Minimum Jaccard weight to keep an edge.
#' @param type Edge weighting passed to [bluster::makeSNNGraph()].
#' @param BNPARAM A `BiocNeighborParam` object controlling the KNN search.
#' @param verbose Print progress.
#' @return A named nested list: `graphs[[as.character(n_dims)]][[as.character(k)]]`
#'   holds an `igraph` object.
#' @export
build_snn_graphs <- function(emb, n_dims = c(15, 30, 50, 75, 100),
                             k = seq(10, 40, 5), prune = 1 / 15,
                             type = "jaccard", BNPARAM = NULL, verbose = TRUE) {
  n_dims <- .clean_n_dims(n_dims, ncol(emb))
  graphs <- lapply(n_dims, function(nd) {
    x <- emb[, seq_len(nd), drop = FALSE]
    gk <- lapply(k, function(kk) {
      if (verbose) message("SNN graph: n_dims = ", nd, ", k = ", kk)
      args <- list(x = x, k = kk, type = type)
      if (!is.null(BNPARAM)) args$BNPARAM <- BNPARAM
      g <- do.call(bluster::makeSNNGraph, args)
      if (prune > 0)
        g <- igraph::delete_edges(g, igraph::E(g)[igraph::E(g)$weight < prune])
      g
    })
    names(gk) <- as.character(k)
    gk
  })
  names(graphs) <- as.character(n_dims)
  graphs
}

.clean_n_dims <- function(n_dims, max_dims) {
  n_dims <- as.integer(n_dims)
  if (any(n_dims > max_dims)) {
    warning("Some `n_dims` exceed the ", max_dims,
            " available dimensions; they were truncated to ", max_dims,
            " and duplicated values removed.", call. = FALSE)
    n_dims <- pmin(n_dims, max_dims)
  }
  sort(unique(n_dims))
}

# Cluster one graph at one resolution. Returns integer membership (1-based).
.cluster_graph <- function(g, resolution, algorithm, seed, n_iter) {
  w <- igraph::E(g)$weight
  if (algorithm == "leiden") {
    res <- leidenbase::leiden_find_partition(
      g,
      partition_type = "RBConfigurationVertexPartition",
      edge_weights = w,
      resolution_parameter = resolution,
      num_iter = n_iter,
      seed = seed
    )
    memb <- as.integer(res$membership)
  } else {
    set.seed(seed)
    memb <- as.integer(igraph::membership(
      igraph::cluster_louvain(g, weights = w, resolution = resolution)))
  }
  # Make labels 1..K regardless of the backend's convention (0- or 1-based)
  match(memb, sort(unique(memb)))
}

# Clusters one graph at one resolution and computes the run metrics.
# Defined at top level (not inside run_clustering_grid()) on purpose: a
# closure would be serialised together with its enclosing frame and every
# worker would receive the whole single-cell object and all the graphs.
.grid_worker <- function(resolution, seed, g, x_sil, algorithm, n_iter) {
  memb <- .cluster_graph(g, resolution, algorithm, seed, n_iter)
  mod <- igraph::modularity(g, membership = memb,
                            weights = igraph::E(g)$weight)
  sil <- NA_real_
  if (!is.null(x_sil) && length(unique(memb)) > 1)
    sil <- mean(bluster::approxSilhouette(x_sil, memb)$width)
  list(membership = memb, modularity = mod, silhouette = sil)
}

#' Run graph clustering over a grid of parameters
#'
#' For every combination of number of PCs (`n_dims`), neighbours (`k`) and
#' `resolutions` the cells are clustered with Leiden (Seurat's
#' `RBConfigurationVertexPartition`, i.e. modularity with a resolution
#' parameter) or Louvain.
#'
#' Parallelisation uses the **future** framework: set a plan *before*
#' calling, e.g. `future::plan(future::multisession, workers = 4)`. The
#' package never changes the plan itself. Results are reproducible for a
#' given `seed` independently of the plan.
#'
#' @param object A `SingleCellExperiment`, `Seurat` object or a cells x dims
#'   matrix.
#' @param reduction Reduction to use (default `"PCA"` / `"pca"`).
#' @param n_dims,k See [build_snn_graphs()].
#' @param resolutions Numeric vector of resolutions.
#' @param algorithm `"leiden"` or `"louvain"`.
#' @param prune SNN pruning threshold (Seurat default 1/15).
#' @param n_iter Leiden iterations (leidenbase default 2; Seurat uses 10).
#' @param silhouette Compute the mean approximate silhouette
#'   ([bluster::approxSilhouette()]) of each run in the PCs used. Adds time.
#' @param seed Integer seed.
#' @param graphs Optional pre-computed output of [build_snn_graphs()].
#' @param BNPARAM Passed to [build_snn_graphs()].
#' @param verbose Print progress.
#' @return A `coclust_grid` object: a list with
#'   * `labels`: integer matrix (cells x runs) of cluster labels,
#'   * `metrics`: data.frame with one row per run (`run`, `k`, `resolution`,
#'     `n_dims`, `modularity`, `silhouette`, `n_clusters`),
#'   * `params`: list of the parameters used.
#' @export
run_clustering_grid <- function(object, reduction = NULL,
                                n_dims = c(15, 30, 50, 75, 100),
                                k = seq(10, 40, 5),
                                resolutions = seq(0.1, 1.4, 0.1),
                                algorithm = c("leiden", "louvain"),
                                prune = 1 / 15, n_iter = 2,
                                silhouette = FALSE, seed = 11,
                                graphs = NULL, BNPARAM = NULL,
                                verbose = TRUE) {
  algorithm <- match.arg(tolower(algorithm), c("leiden", "louvain"))
  emb <- get_embedding(object, reduction)
  n_dims <- .clean_n_dims(n_dims, ncol(emb))
  # Avoid floating point names such as 0.30000000000000004
  resolutions <- round(resolutions, 10)

  if (is.null(graphs)) {
    graphs <- build_snn_graphs(emb, n_dims = n_dims, k = k, prune = prune,
                               BNPARAM = BNPARAM, verbose = verbose)
  }

  grid <- expand.grid(resolution = resolutions, k = k, n_dims = n_dims,
                      KEEP.OUT.ATTRS = FALSE)
  grid$run <- sprintf("K%s_Res%s_Dims%s", grid$k, grid$resolution, grid$n_dims)
  grid$seed <- as.integer(seed) + seq_len(nrow(grid)) - 1L

  # Parallelise over resolutions within each graph, so that each worker only
  # receives the one graph it needs (instead of all graphs at once).
  out <- vector("list", nrow(grid))
  for (nd in n_dims) {
    x_sil <- if (silhouette) emb[, seq_len(nd), drop = FALSE] else NULL
    for (kk in k) {
      g <- graphs[[as.character(nd)]][[as.character(kk)]]
      if (is.null(g)) stop("Graph for n_dims = ", nd, ", k = ", kk,
                           " not found in `graphs`.", call. = FALSE)
      idx <- which(grid$n_dims == nd & grid$k == kk)
      if (verbose) message("Clustering n_dims = ", nd, ", k = ", kk,
                           " (", length(idx), " resolutions)")
      # Only the graph and the per-run parameters are sent to the workers.
      # .grid_worker() lives in the package namespace, so it does not carry
      # the frame of this function (object, embedding, all graphs) with it.
      out[idx] <- future.apply::future_mapply(
        .grid_worker,
        resolution = grid$resolution[idx],
        seed = grid$seed[idx],
        MoreArgs = list(g = g, x_sil = x_sil, algorithm = algorithm,
                        n_iter = n_iter),
        SIMPLIFY = FALSE, USE.NAMES = FALSE,
        future.seed = TRUE,  # avoids RNG warnings; seeds are set explicitly
        future.packages = c("igraph", "leidenbase", "bluster")
      )
    }
  }

  labels <- vapply(out, `[[`, integer(nrow(emb)), "membership")
  if (is.null(dim(labels))) labels <- matrix(labels, nrow = nrow(emb))
  dimnames(labels) <- list(rownames(emb), grid$run)

  metrics <- data.frame(
    run = grid$run, k = grid$k, resolution = grid$resolution,
    n_dims = grid$n_dims,
    modularity = vapply(out, `[[`, numeric(1), "modularity"),
    silhouette = vapply(out, `[[`, numeric(1), "silhouette"),
    n_clusters = apply(labels, 2, function(z) length(unique(z))),
    stringsAsFactors = FALSE
  )

  structure(
    list(labels = labels, metrics = metrics,
         params = list(reduction = reduction, n_dims = n_dims, k = k,
                       resolutions = resolutions, algorithm = algorithm,
                       prune = prune, n_iter = n_iter, seed = seed)),
    class = "coclust_grid"
  )
}

#' @export
print.coclust_grid <- function(x, ...) {
  cat("<coclust_grid>", nrow(x$labels), "cells x", ncol(x$labels), "runs\n")
  cat("  algorithm:", x$params$algorithm, "\n")
  cat("  n_dims:", paste(x$params$n_dims, collapse = ", "), "\n")
  cat("  k:", paste(x$params$k, collapse = ", "), "\n")
  cat("  resolutions:", paste(x$params$resolutions, collapse = ", "), "\n")
  cat("  clusters per run: ", paste(range(x$metrics$n_clusters), collapse = "-"),
      "\n", sep = "")
  invisible(x)
}
