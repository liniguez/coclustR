#' Collapse cells with identical label vectors into meta-barcodes
#'
#' Two cells that received the same label in every run are indistinguishable
#' for the consensus; collapsing them makes the co-clustering matrix much
#' smaller.
#'
#' @param labels Cells x runs matrix of labels (e.g. `grid$labels`). `NA`
#'   means "unassigned in this run".
#' @return A list with
#'   * `labels`: meta-barcodes x runs matrix,
#'   * `sizes`: number of cells in each meta-barcode,
#'   * `representative`: one cell name per meta-barcode,
#'   * `cell2meta`: named integer vector mapping every cell to its
#'     meta-barcode (row index of `labels`).
#' @export
collapse_metabarcodes <- function(labels) {
  if (is.null(rownames(labels))) stop("`labels` needs cell names as rownames.")
  # A separator is essential: without it, labels (1, 12) and (11, 2) would
  # produce the same key "112" (bug in the original script).
  key <- do.call(paste, c(as.data.frame(labels, stringsAsFactors = FALSE),
                          sep = "\r"))
  first <- !duplicated(key)
  cell2meta <- match(key, key[first])
  meta_ids <- sprintf("MB%0*d", nchar(sum(first)), seq_len(sum(first)))
  meta_labels <- labels[first, , drop = FALSE]
  representative <- rownames(labels)[first]
  rownames(meta_labels) <- meta_ids
  list(labels = meta_labels,
       sizes = stats::setNames(tabulate(cell2meta, nbins = sum(first)), meta_ids),
       representative = stats::setNames(representative, meta_ids),
       cell2meta = stats::setNames(cell2meta, rownames(labels)))
}

#' Co-clustering (consensus) matrix
#'
#' Entry (i, j) is the fraction of runs in which i and j were assigned to the
#' same cluster, among the runs in which both were assigned (non-`NA`).
#'
#' Uses block updates (`co[idx, idx] + 1` for each cluster) which costs
#' `sum(cluster_size^2)` per run instead of `n^2 * n_clusters` for the
#' indicator-matrix product in the original script.
#'
#' @param labels Rows = items (cells or meta-barcodes), columns = runs.
#' @return Symmetric numeric matrix with values in \[0, 1\] (`NA` where two
#'   items were never assigned together).
#' @export
coclustering_matrix <- function(labels) {
  labels <- as.matrix(labels)
  n <- nrow(labels)
  co <- matrix(0, n, n)
  has_na <- anyNA(labels)
  valid <- if (has_na) matrix(0, n, n) else NULL

  for (j in seq_len(ncol(labels))) {
    lab <- labels[, j]
    ok <- which(!is.na(lab))
    groups <- split(ok, lab[ok])
    for (idx in groups) co[idx, idx] <- co[idx, idx] + 1
    if (has_na) valid[ok, ok] <- valid[ok, ok] + 1
  }

  if (has_na) {
    co <- co / valid
    co[valid == 0] <- NA_real_
  } else {
    co <- co / ncol(labels)
  }
  dimnames(co) <- list(rownames(labels), rownames(labels))
  co
}

#' Core stability of each meta-barcode
#'
#' The "core" of a meta-barcode is the set of meta-barcodes it co-clusters
#' with in at least `threshold` of the runs (including itself).
#'
#' * `core_cells`: number of *cells* in the core (itself included).
#' * `core_stability`: mean co-clustering with the *other* members of the
#'   core (the diagonal is excluded; 0 if the core has no other member).
#' * `core_stability_w`: same, weighted by the number of cells of each member.
#'
#' Excluding the diagonal fixes a problem of the original script, where an
#' isolated meta-barcode (nobody above the threshold) got a stability of 1
#' (its self-similarity) and therefore was never flagged as unstable.
#'
#' @param co Co-clustering matrix from [coclustering_matrix()].
#' @param sizes Cells per row of `co`.
#' @param threshold Minimum co-clustering to be part of the core.
#' @return data.frame with one row per meta-barcode.
#' @export
core_stability <- function(co, sizes, threshold = 0.5) {
  stopifnot(nrow(co) == length(sizes))
  A <- co >= threshold
  A[is.na(A)] <- FALSE
  diag(A) <- TRUE
  core_cells <- as.vector(A %*% sizes)

  A_off <- A
  diag(A_off) <- FALSE
  co0 <- co
  co0[is.na(co0)] <- 0
  M <- co0 * A_off
  n_other <- rowSums(A_off)
  stab <- ifelse(n_other > 0, rowSums(M) / pmax(n_other, 1), 0)
  w_other <- as.vector(A_off %*% sizes)
  stab_w <- ifelse(w_other > 0, as.vector(M %*% sizes) / pmax(w_other, 1), 0)

  data.frame(metabarcode = rownames(co), cells = as.integer(sizes),
             core_members = n_other + 1L, core_cells = core_cells,
             core_stability = stab, core_stability_w = stab_w,
             row.names = rownames(co), stringsAsFactors = FALSE)
}

#' Flag unstable meta-barcodes
#'
#' A meta-barcode is flagged when its core is small **and** weakly
#' co-clustered: `core_cells < min_core_cells & core_stability < min_core_stability`
#' (same rule as the original script, now with explicit parameters).
#'
#' @param stability Output of [core_stability()].
#' @param min_core_cells Minimum number of cells in the core.
#' @param min_core_stability Minimum core stability.
#' @param weighted Use `core_stability_w` instead of `core_stability`.
#' @return Logical vector (`TRUE` = unstable), named by meta-barcode.
#' @export
flag_unstable <- function(stability, min_core_cells = 50,
                          min_core_stability = 0.8, weighted = FALSE) {
  s <- if (weighted) stability$core_stability_w else stability$core_stability
  stats::setNames(stability$core_cells < min_core_cells &
                    s < min_core_stability,
                  stability$metabarcode)
}
