#' Consensus clusters from a co-clustering matrix
#'
#' Builds a hierarchical tree on `1 - co` and cuts it with
#' [dynamicTreeCut::cutreeDynamicTree()] and/or
#' [dynamicTreeCut::cutreeHybrid()].
#'
#' The tree is built on meta-barcodes. `sizes` is passed to
#' [stats::hclust()] as `members`, so for size-aware linkages (`"average"`,
#' `"ward.D2"`, ...) the tree is the one you would get on the individual
#' cells. (For `"complete"` and `"single"` the sizes do not matter.)
#' Note that `min_cluster_size` in dynamicTreeCut counts **meta-barcodes**,
#' not cells; use `min_cluster_cells` to additionally require a minimum
#' number of cells per consensus cluster.
#'
#' @param co Co-clustering matrix (meta-barcodes x meta-barcodes).
#' @param sizes Cells per meta-barcode.
#' @param linkage Linkage for [stats::hclust()].
#' @param methods Any of `"hybrid"`, `"dynamic"`.
#' @param min_cluster_size Minimum module size (in meta-barcodes) for
#'   dynamicTreeCut.
#' @param min_cluster_cells Optional minimum number of cells; smaller
#'   consensus clusters are set to 0 (unstable). `NULL` to skip.
#' @param deep_split Sensitivity of the cut. Logical for `"dynamic"`;
#'   for `"hybrid"` it can also be an integer 0-4.
#' @param pam_respects_dendro Passed to `cutreeHybrid()`.
#' @param max_tree_height `maxTreeHeight` for `cutreeDynamicTree()`. The
#'   default of dynamicTreeCut (1) lets it merge groups that *never*
#'   co-clustered (distance exactly 1); 0.99 requires every pair in a
#'   complete-linkage cluster to have co-clustered in at least 1% of runs.
#' @param rescue_min_cells Rescue of very stable populations (see Details).
#'   `NULL` disables it.
#' @param rescue_threshold Co-clustering threshold that defines a rescue
#'   group (complete-linkage groups in which every pair co-clusters in at
#'   least this fraction of runs).
#'
#' @details
#' **Rescue of stable populations.** A population whose cells get the same
#' label in every run collapses into a *single* meta-barcode. Since
#' dynamicTreeCut counts leaves (meta-barcodes), such a population can never
#' reach `min_cluster_size`; DynamicTree leaves it unassigned and the PAM
#' stage of Hybrid merges it into the nearest cluster even if they never
#' co-cluster. These are usually the most distinct (often rare) cell types.
#' With `rescue_min_cells`, groups of meta-barcodes that co-cluster in at
#' least `rescue_threshold` of the runs, that were left unassigned by the tree
#' cut (before PAM for Hybrid) and that contain at least `rescue_min_cells`
#' cells become clusters of their own.
#' @return A list with `tree` (hclust) and `labels` (integer matrix
#'   meta-barcodes x methods; 0 = unassigned).
#' @export
consensus_clusters <- function(co, sizes = rep(1L, nrow(co)),
                               linkage = "complete",
                               methods = c("hybrid", "dynamic"),
                               min_cluster_size = 50,
                               min_cluster_cells = NULL,
                               deep_split = FALSE,
                               pam_respects_dendro = FALSE,
                               max_tree_height = 0.99,
                               rescue_min_cells = 50,
                               rescue_threshold = 0.5) {
  methods <- match.arg(methods, c("hybrid", "dynamic"), several.ok = TRUE)
  D <- 1 - co
  D[is.na(D)] <- 1
  diag(D) <- 0
  d <- stats::as.dist(D)
  tree <- stats::hclust(d, method = linkage, members = as.numeric(sizes))

  # dynamicTreeCut needs at least 2 leaves and a minimum size <= n leaves
  mcs <- min(min_cluster_size, max(1L, nrow(co) - 1L))
  ds_h <- if (is.logical(deep_split)) as.integer(deep_split) else deep_split
  # dynamicTreeCut prints "No module detected" with print(): silence it
  quiet <- function(expr) {
    utils::capture.output(out <- suppressMessages(expr))
    out
  }
  hybrid <- function(pam) as.integer(quiet(dynamicTreeCut::cutreeHybrid(
    dendro = tree, distM = D, deepSplit = ds_h, minClusterSize = mcs,
    pamStage = pam, pamRespectsDendro = pam_respects_dendro, verbose = 0
  )$labels))

  labs <- list(); pre <- list()
  if ("hybrid" %in% methods) {
    labs$Hybrid <- hybrid(TRUE)
    pre$Hybrid <- if (is.null(rescue_min_cells)) labs$Hybrid else hybrid(FALSE)
  }
  if ("dynamic" %in% methods) {
    labs$DynamicTree <- as.integer(quiet(
      dynamicTreeCut::cutreeDynamicTree(
        dendro = tree, maxTreeHeight = max_tree_height,
        deepSplit = as.logical(deep_split), minModuleSize = mcs)))
    pre$DynamicTree <- labs$DynamicTree
  }

  if (!is.null(rescue_min_cells)) {
    groups <- stats::cutree(stats::hclust(d, method = "complete"),
                            h = 1 - rescue_threshold)
    for (m in names(labs)) {
      labs[[m]] <- .rescue_groups(labs[[m]], pre[[m]], groups, sizes,
                                  rescue_min_cells)
    }
  }

  if (!is.null(min_cluster_cells)) {
    labs <- lapply(labs, function(l) {
      tot <- tapply(sizes, l, sum)
      small <- as.integer(names(tot)[tot < min_cluster_cells])
      l[l %in% small] <- 0L
      l
    })
  }
  labels <- do.call(cbind, labs)
  rownames(labels) <- rownames(co)
  list(tree = tree, labels = labels)
}

# Groups entirely unassigned before PAM and large enough (in cells) become
# new clusters, numbered after the existing ones.
.rescue_groups <- function(final, pre, groups, sizes, min_cells) {
  unassigned <- tapply(pre == 0L, groups, all)
  cells <- tapply(sizes, groups, sum)
  ok <- names(which(unassigned & cells >= min_cells))
  ok <- ok[order(cells[ok], decreasing = TRUE)]
  nxt <- max(c(0L, final))
  for (gid in ok) {
    nxt <- nxt + 1L
    final[groups == as.integer(gid)] <- nxt
  }
  final
}

# Turn integer labels (0/NA = unassigned) into a factor with ordered levels
# 1..K followed by "Unstable". Fixes the original script, which built levels
# as 1:(n_unique - 1) and silently produced NAs when no cell was unstable or
# when labels were not consecutive.
.labels_to_factor <- function(x, unstable = "Unstable") {
  x <- as.integer(x)
  bad <- is.na(x) | x == 0L
  lv <- sort(unique(x[!bad]))
  out <- as.character(x)
  out[bad] <- unstable
  factor(out, levels = c(as.character(lv), unstable))
}
