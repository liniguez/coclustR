sim_embedding <- function(n_per = c(150, 120, 100, 80), n_dims = 20, sep = 6,
                          seed = 1) {
  set.seed(seed)
  K <- length(n_per)
  centers <- matrix(rnorm(K * n_dims, sd = sep), K, n_dims)
  truth <- rep(seq_len(K), n_per)
  x <- centers[truth, ] + matrix(rnorm(sum(n_per) * n_dims), ncol = n_dims)
  rownames(x) <- paste0("cell", seq_len(nrow(x)))
  colnames(x) <- paste0("PC_", seq_len(n_dims))
  list(x = x, truth = truth)
}

sim_sce <- function(...) {
  s <- sim_embedding(...)
  counts <- matrix(rpois(50 * nrow(s$x), 1), 50,
                   dimnames = list(paste0("g", 1:50), rownames(s$x)))
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  SingleCellExperiment::reducedDim(sce, "PCA") <- s$x
  SingleCellExperiment::reducedDim(sce, "UMAP") <- s$x[, 1:2]
  sce$truth <- factor(s$truth)
  sce$Sample <- rep(c("A", "B"), length.out = ncol(sce))
  sce
}

sim_seurat <- function(...) {
  s <- sim_embedding(...)
  counts <- Matrix::Matrix(matrix(rpois(50 * nrow(s$x), 1), 50,
                                  dimnames = list(paste0("g", 1:50), rownames(s$x))),
                           sparse = TRUE)
  seu <- SeuratObject::CreateSeuratObject(counts)
  seu[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = s$x, key = "PC_", assay = "RNA")
  u <- s$x[, 1:2]; colnames(u) <- c("umap_1", "umap_2")
  seu[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = u, key = "umap_", assay = "RNA")
  seu$truth <- factor(s$truth)
  seu
}

small_grid_args <- list(n_dims = c(5, 10, 20), k = c(10, 20),
                        resolutions = c(0.2, 0.5, 0.8), verbose = FALSE,
                        min_cluster_size = 5, min_core_cells = 10)
