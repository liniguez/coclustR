#' @importFrom ggplot2 ggplot aes .data
NULL

.as_result <- function(x) {
  if (inherits(x, c("coclust_result", "coclust_grid"))) return(x)
  get_coclust_result(x)
}

.metrics <- function(x) {
  x <- .as_result(x)
  m <- if (inherits(x, "coclust_grid")) x$metrics else x$grid$metrics
  m$n_dims_f <- factor(paste0(m$n_dims, " PCs"),
                       levels = paste0(sort(unique(m$n_dims)), " PCs"))
  m
}

#' Colour palette for cluster labels ("Unstable" in grey)
#'
#' @param levels Character vector of cluster levels.
#' @param palette Any palette name accepted by [grDevices::hcl.colors()].
#' @param unstable Label of the unassigned group.
#' @return Named character vector of colours.
#' @export
cluster_palette <- function(levels, palette = "viridis", unstable = "Unstable") {
  lv <- setdiff(levels, unstable)
  cols <- stats::setNames(grDevices::hcl.colors(max(length(lv), 1), palette)[seq_along(lv)], lv)
  if (unstable %in% levels) cols[unstable] <- "grey75"
  cols
}

#' Overview of the parameter grid
#'
#' One panel per number of PCs; x = k, y = resolution, colour = number of
#' clusters (also printed), size = modularity.
#'
#' @param x A `coclust_grid`, `coclust_result` or an object processed with
#'   [run_coclustering()].
#' @return A ggplot.
#' @export
plot_grid_overview <- function(x) {
  m <- .metrics(x)
  ggplot(m, aes(x = .data$k, y = .data$resolution)) +
    ggplot2::geom_point(aes(colour = .data$n_clusters, size = .data$modularity)) +
    ggplot2::geom_text(aes(label = .data$n_clusters), colour = "white",
                       size = 2.5) +
    ggplot2::scale_colour_viridis_c(name = "Clusters") +
    ggplot2::scale_size_continuous(name = "Modularity", range = c(2, 7)) +
    ggplot2::facet_wrap(~ n_dims_f) +
    ggplot2::labs(x = "k (neighbours)", y = "Resolution") +
    ggplot2::theme_bw()
}

#' Clustering metrics along resolution
#'
#' Rows = metric, columns = number of PCs, one line per `k`.
#'
#' @inheritParams plot_grid_overview
#' @param metrics Any of `"modularity"`, `"silhouette"`, `"n_clusters"`.
#'   Metrics that are all `NA` (e.g. silhouette not computed) are dropped.
#' @return A ggplot.
#' @export
plot_grid_metrics <- function(x, metrics = c("modularity", "silhouette",
                                             "n_clusters")) {
  m <- .metrics(x)
  metrics <- metrics[vapply(metrics, function(v) !all(is.na(m[[v]])), TRUE)]
  long <- do.call(rbind, lapply(metrics, function(v)
    data.frame(m[, c("k", "resolution", "n_dims_f")], metric = v,
               value = m[[v]])))
  long$metric <- factor(long$metric, levels = metrics)
  ggplot(long, aes(x = .data$resolution, y = .data$value,
                   colour = factor(.data$k), group = .data$k)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_grid(metric ~ n_dims_f, scales = "free_y", switch = "y") +
    ggplot2::scale_colour_viridis_d(name = "k") +
    ggplot2::labs(x = "Resolution", y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
                   strip.placement = "outside",
                   axis.title = ggplot2::element_text(face = "bold"))
}

#' Core stability diagnostic
#'
#' Each point is a meta-barcode: cells in its core vs core stability; the
#' dashed lines are the thresholds used by [flag_unstable()].
#'
#' @inheritParams plot_grid_overview
#' @return A ggplot.
#' @export
plot_stability <- function(x) {
  res <- .as_result(x)
  if (!inherits(res, "coclust_result")) stop("Needs a coclust_result.")
  s <- res$stability
  s$status <- ifelse(res$unstable[s$metabarcode], "Unstable", "Kept")
  ggplot(s, aes(x = .data$core_cells, y = .data$core_stability)) +
    ggplot2::geom_point(aes(colour = .data$status, size = .data$cells),
                        alpha = 0.6) +
    ggplot2::geom_vline(xintercept = res$params$min_core_cells, linetype = 2) +
    ggplot2::geom_hline(yintercept = res$params$min_core_stability, linetype = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c(Kept = "#2B2D42", Unstable = "#E63946")) +
    ggplot2::scale_size_area(max_size = 5, name = "Cells") +
    ggplot2::labs(x = "Cells in core (log10)", y = "Core stability",
                  colour = NULL) +
    ggplot2::theme_bw()
}

#' Plot cells on an embedding coloured by (consensus) clusters
#'
#' Works with SingleCellExperiment and Seurat objects.
#'
#' @param object A `SingleCellExperiment` or `Seurat` object.
#' @param color_by One or more metadata columns; several columns are shown as
#'   facets.
#' @param reduction Embedding (default `"UMAP"` / `"umap"`).
#' @param dims Which two dimensions to plot.
#' @param label Print the label at the median position of each group.
#' @param point_size Point size.
#' @param palette Passed to [cluster_palette()].
#' @return A ggplot.
#' @export
plot_clusters_embedding <- function(object, color_by = "coclust_Hybrid",
                                    reduction = NULL, dims = c(1, 2),
                                    label = TRUE, point_size = 0.5,
                                    palette = "viridis") {
  type <- .object_type(object)
  if (is.null(reduction)) reduction <- .default_reduction(type, "umap")
  emb <- get_embedding(object, reduction)[, dims, drop = FALSE]
  md <- get_cell_metadata(object)
  miss <- setdiff(color_by, colnames(md))
  if (length(miss)) stop("Columns not found in metadata: ",
                         paste(miss, collapse = ", "), call. = FALSE)
  long <- do.call(rbind, lapply(color_by, function(v)
    data.frame(x = emb[, 1], y = emb[, 2], panel = v,
               group = as.character(md[rownames(emb), v]))))
  lv <- unique(unlist(lapply(color_by, function(v) {
    z <- md[[v]]
    if (is.factor(z)) levels(z) else sort(unique(as.character(z)))
  })))
  lv <- c(setdiff(lv, "Unstable"), intersect(lv, "Unstable"))
  long$group <- factor(long$group, levels = lv)
  long$panel <- factor(long$panel, levels = color_by)
  long <- long[order(long$group == "Unstable", decreasing = TRUE), ]  # grey below

  p <- ggplot(long, aes(x = .data$x, y = .data$y, colour = .data$group)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_colour_manual(values = cluster_palette(lv, palette),
                                 name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 3))) +
    ggplot2::labs(x = paste0(reduction, "_", dims[1]),
                  y = paste0(reduction, "_", dims[2])) +
    ggplot2::theme_bw()
  if (length(color_by) > 1) p <- p + ggplot2::facet_wrap(~ panel)
  if (label && any(long$group != "Unstable")) {
    cent <- stats::aggregate(cbind(x, y) ~ group + panel,
                             data = long[long$group != "Unstable", ],
                             FUN = stats::median)
    p <- p + ggplot2::geom_text(data = cent, aes(label = .data$group),
                                colour = "black", fontface = "bold")
  }
  p
}

#' Heatmap of the co-clustering matrix on a sample of cells
#'
#' Cells are sampled (optionally stratified, e.g. by sample and cluster),
#' each cell takes the row of its meta-barcode, and rows/columns are ordered
#' by a hierarchical tree on `1 - co` using the same linkage as the
#' consensus.
#'
#' @param x A `coclust_result` or an object processed with
#'   [run_coclustering()] (with `store_matrix = TRUE`).
#' @param object Optional SingleCellExperiment / Seurat object to take extra
#'   annotation columns from (defaults to `x` when `x` is such an object).
#' @param n_cells Number of cells to sample (total).
#' @param stratify_by Metadata columns to stratify the sampling by
#'   (each stratum contributes proportionally to its size, at least 1 cell).
#'   Consensus columns (`"Hybrid"`, `"DynamicTree"`) can be used directly.
#' @param annotate_by Extra metadata columns to show as annotations.
#' @param annotation_colors Named list of named colour vectors, as in
#'   [ComplexHeatmap::HeatmapAnnotation()].
#' @param engine `"ComplexHeatmap"` (default) or `"ggplot"`.
#' @param seed Random seed for the sampling.
#' @return A `Heatmap` object or a ggplot.
#' @export
plot_coclustering_heatmap <- function(x, object = NULL, n_cells = 1000,
                                      stratify_by = NULL, annotate_by = NULL,
                                      annotation_colors = list(),
                                      engine = c("ComplexHeatmap", "ggplot"),
                                      seed = 1) {
  engine <- match.arg(engine)
  if (is.null(object) && !inherits(x, "coclust_result")) object <- x
  res <- .as_result(x)
  if (is.null(res$coclust))
    stop("The result was stored without the co-clustering matrix ",
         "(store_matrix = FALSE).", call. = FALSE)

  ann <- res$cell_clusters
  if (!is.null(object) && length(c(stratify_by, annotate_by))) {
    md <- get_cell_metadata(object)
    extra <- setdiff(c(stratify_by, annotate_by), colnames(ann))
    miss <- setdiff(extra, colnames(md))
    if (length(miss)) stop("Columns not found: ", paste(miss, collapse = ", "),
                           call. = FALSE)
    ann <- cbind(ann, md[rownames(ann), extra, drop = FALSE])
  }

  # --- stratified sampling of cells
  set.seed(seed)
  cells <- rownames(ann)
  n_cells <- min(n_cells, length(cells))
  if (is.null(stratify_by)) {
    picked <- sample(cells, n_cells)
  } else {
    strata <- interaction(ann[, stratify_by, drop = FALSE], drop = TRUE)
    by_s <- split(cells, strata)
    n_s <- pmax(1, round(lengths(by_s) / length(cells) * n_cells))
    picked <- unlist(mapply(function(cc, nn) cc[sample.int(length(cc), min(nn, length(cc)))],
                            by_s, n_s, SIMPLIFY = FALSE), use.names = FALSE)
  }

  mb <- res$metabarcodes$cell2meta[picked]
  M <- res$coclust[mb, mb, drop = FALSE]
  M[is.na(M)] <- 0
  dimnames(M) <- list(picked, picked)
  tree <- stats::hclust(stats::as.dist(1 - M), method = res$params$linkage)

  ann_cols <- c(colnames(res$cell_clusters), setdiff(annotate_by, colnames(res$cell_clusters)))
  A <- ann[picked, ann_cols, drop = FALSE]
  for (m in colnames(res$cell_clusters)) {
    if (is.null(annotation_colors[[m]]))
      annotation_colors[[m]] <- cluster_palette(levels(A[[m]]))
  }

  if (engine == "ComplexHeatmap") {
    .need("ComplexHeatmap", "for engine = 'ComplexHeatmap'")
    .need("circlize", "for engine = 'ComplexHeatmap'")
    col_fun <- circlize::colorRamp2(seq(0, 1, 0.1),
                                    grDevices::hcl.colors(11, "viridis"))
    annotation_colors <- annotation_colors[intersect(names(annotation_colors),
                                                     colnames(A))]
    top <- ComplexHeatmap::HeatmapAnnotation(df = A, col = annotation_colors)
    return(ComplexHeatmap::Heatmap(
      M, name = "Co-clustering", col = col_fun,
      cluster_rows = tree, cluster_columns = tree,
      show_row_names = FALSE, show_column_names = FALSE,
      top_annotation = top))
  }

  ord <- tree$labels[tree$order]
  long <- data.frame(c1 = factor(rep(rownames(M), ncol(M)), levels = ord),
                     c2 = factor(rep(colnames(M), each = nrow(M)),
                                 levels = rev(ord)),
                     value = as.vector(M))
  ggplot(long, aes(x = .data$c1, y = .data$c2, fill = .data$value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(name = "Co-clustering", limits = c(0, 1)) +
    ggplot2::labs(x = "Cells", y = "Cells") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
}
