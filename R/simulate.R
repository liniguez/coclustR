#' Simulate a single-cell count matrix with known populations
#'
#' Negative-binomial counts for a set of populations defined by marker genes.
#' Populations 1 and 2 share part of their markers (so they are hard to
#' separate and split differently across resolutions), and the last
#' population is small and very distinct (a "rare cell type", useful to see
#' the rescue step of [consensus_clusters()]).
#'
#' @param n_cells Cells per population.
#' @param n_genes Number of genes.
#' @param n_markers Marker genes per population.
#' @param fold_change Fold change of the markers.
#' @param rare_fold_change Fold change of the markers of the last population.
#' @param shared Fraction of the markers of population 2 taken from population 1.
#' @param n_samples Number of samples (cells are assigned at random).
#' @param dispersion Negative-binomial size parameter.
#' @param seed Random seed.
#' @return A list with `counts` (genes x cells integer matrix) and `cell_data`
#'   (data.frame with `population` and `Sample`, rownames = cell names).
#' @examples
#' sim <- simulate_sc_counts()
#' dim(sim$counts)
#' table(sim$cell_data$population)
#' @export
simulate_sc_counts <- function(n_cells = c(400, 300, 250, 150, 60),
                               n_genes = 1000, n_markers = 40,
                               fold_change = 2.5, rare_fold_change = 5,
                               shared = 0.6, n_samples = 2,
                               dispersion = 2, seed = 1) {
  set.seed(seed)
  K <- length(n_cells)
  base <- stats::rgamma(n_genes, shape = 0.6, rate = 0.3)
  mu <- matrix(base, n_genes, K)
  markers <- vector("list", K)
  for (k in seq_len(K)) {
    m <- sample(n_genes, n_markers)
    if (k == 2 && K > 2) {
      n_sh <- round(shared * n_markers)
      m <- c(markers[[1]][seq_len(n_sh)], m[seq_len(n_markers - n_sh)])
    }
    markers[[k]] <- m
    fc <- if (k == K && K > 1) rare_fold_change else fold_change
    mu[m, k] <- mu[m, k] * fc
  }
  pop <- rep(seq_len(K), n_cells)
  N <- length(pop)
  lib <- stats::rlnorm(N, 0, 0.3)
  counts <- matrix(stats::rnbinom(n_genes * N, size = dispersion,
                                  mu = mu[, pop] * rep(lib, each = n_genes)),
                   n_genes, N)
  cells <- sprintf("Cell%04d", seq_len(N))
  dimnames(counts) <- list(sprintf("Gene%04d", seq_len(n_genes)), cells)
  cell_data <- data.frame(
    population = factor(paste0("Pop", pop)),
    Sample = factor(sample(paste0("S", seq_len(n_samples)), N, replace = TRUE)),
    row.names = cells)
  list(counts = counts, cell_data = cell_data)
}
