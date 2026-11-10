test_that("SingleCellExperiment workflow", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- sim_sce()
  sce <- do.call(run_coclustering, c(list(object = sce), small_grid_args))
  expect_true(all(c("coclust_Hybrid", "coclust_DynamicTree") %in%
                    colnames(SummarizedExperiment::colData(sce))))
  expect_s3_class(get_coclust_result(sce), "coclust_result")
  expect_s3_class(plot_grid_overview(sce), "ggplot")
  expect_s3_class(plot_grid_metrics(sce), "ggplot")
  expect_s3_class(plot_stability(sce), "ggplot")
  expect_s3_class(plot_clusters_embedding(
    sce, color_by = c("coclust_Hybrid", "coclust_DynamicTree")), "ggplot")
  expect_s3_class(plot_coclustering_heatmap(sce, n_cells = 200,
                                            stratify_by = c("Sample", "Hybrid"),
                                            annotate_by = "Sample",
                                            engine = "ggplot"), "ggplot")
  skip_if_not_installed("ComplexHeatmap")
  expect_s4_class(plot_coclustering_heatmap(sce, n_cells = 200,
                                            annotate_by = "Sample"), "Heatmap")
})

test_that("Seurat workflow gives the same clusters as SCE", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("SingleCellExperiment")
  seu <- sim_seurat()
  sce <- sim_sce()
  seu <- do.call(run_coclustering, c(list(object = seu), small_grid_args))
  sce <- do.call(run_coclustering, c(list(object = sce), small_grid_args))
  expect_true("coclust_Hybrid" %in% colnames(seu[[]]))
  expect_equal(as.character(seu$coclust_Hybrid),
               as.character(sce$coclust_Hybrid))
  expect_s3_class(get_coclust_result(seu), "coclust_result")
  expect_s3_class(plot_clusters_embedding(seu), "ggplot")
  expect_s3_class(plot_grid_metrics(seu), "ggplot")
})

test_that("re-running only the consensus from a stored grid", {
  s <- sim_embedding()
  g <- run_clustering_grid(s$x, n_dims = c(5, 10), k = c(10, 20),
                           resolutions = c(.3, .8), verbose = FALSE)
  r1 <- coclust_from_grid(g, min_cluster_size = 5, min_core_cells = 10,
                          verbose = FALSE)
  r2 <- coclust_from_grid(g, min_cluster_size = 5, min_core_cells = 10,
                          linkage = "average", verbose = FALSE)
  expect_equal(nrow(r1$cell_clusters), nrow(s$x))
  expect_equal(nrow(r2$cell_clusters), nrow(s$x))
})

test_that("parallel future plan gives identical results", {
  skip_if_not_installed("future")
  skip_on_cran()
  s <- sim_embedding()
  args <- list(object = s$x, n_dims = c(5, 10), k = 15,
               resolutions = c(.3, .8), verbose = FALSE, seed = 5)
  future::plan(future::sequential)
  a <- do.call(run_clustering_grid, args)
  future::plan(future::multisession, workers = 2)
  on.exit(future::plan(future::sequential), add = TRUE)
  b <- do.call(run_clustering_grid, args)
  expect_identical(a$labels, b$labels)
})
