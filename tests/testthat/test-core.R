test_that("metabarcode keys do not collide (paste0 bug)", {
  # (1, 12) and (11, 2) would both give "112" without a separator
  lab <- matrix(c(1L, 11L, 12L, 2L), 2, dimnames = list(c("a", "b"), NULL))
  mb <- collapse_metabarcodes(lab)
  expect_equal(nrow(mb$labels), 2)
  expect_equal(unname(mb$sizes), c(1, 1))
})

test_that("collapse_metabarcodes maps every cell back", {
  lab <- matrix(c(1, 1, 2, 2, 1, 1, 3, 3), 4,
                dimnames = list(paste0("c", 1:4), NULL))
  mb <- collapse_metabarcodes(lab)
  expect_equal(nrow(mb$labels), 2)
  expect_equal(mb$labels[mb$cell2meta, ], lab, ignore_attr = TRUE)
})

test_that("coclustering_matrix matches brute force, with and without NA", {
  set.seed(2)
  lab <- matrix(sample(1:3, 40 * 6, TRUE), 40,
                dimnames = list(paste0("c", 1:40), NULL))
  brute <- function(L) {
    n <- nrow(L); out <- matrix(NA_real_, n, n)
    for (i in 1:n) for (j in 1:n) {
      ok <- !is.na(L[i, ]) & !is.na(L[j, ])
      if (any(ok)) out[i, j] <- mean(L[i, ok] == L[j, ok])
    }
    out
  }
  expect_equal(coclustering_matrix(lab), brute(lab), ignore_attr = TRUE)
  lab[sample(length(lab), 30)] <- NA
  expect_equal(coclustering_matrix(lab), brute(lab), ignore_attr = TRUE)
})

test_that("isolated meta-barcodes get stability 0 and are flagged", {
  co <- matrix(c(1, .9, .1,
                 .9, 1, .1,
                 .1, .1, 1), 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  st <- core_stability(co, sizes = c(30, 30, 5), threshold = 0.5)
  expect_equal(st$core_cells, c(60, 60, 5))
  expect_equal(st$core_stability, c(.9, .9, 0))
  fl <- flag_unstable(st, min_core_cells = 50, min_core_stability = 0.8)
  expect_equal(unname(fl), c(FALSE, FALSE, TRUE))
})

test_that(".labels_to_factor handles no-unstable and gaps", {
  f <- coclustR:::.labels_to_factor(c(1, 2, 2, 3))
  expect_false(anyNA(f))
  expect_equal(levels(f), c("1", "2", "3", "Unstable"))
  f2 <- coclustR:::.labels_to_factor(c(1, 0, 4, NA))
  expect_equal(as.character(f2), c("1", "Unstable", "4", "Unstable"))
})

test_that("leiden uses edge weights and is reproducible", {
  s <- sim_embedding()
  g1 <- run_clustering_grid(s$x, n_dims = 10, k = 15, resolutions = 0.5,
                            verbose = FALSE, seed = 3)
  g2 <- run_clustering_grid(s$x, n_dims = 10, k = 15, resolutions = 0.5,
                            verbose = FALSE, seed = 3)
  expect_identical(g1$labels, g2$labels)
  expect_s3_class(g1, "coclust_grid")
})

test_that("full pipeline recovers simulated clusters (matrix input)", {
  s <- sim_embedding()
  res <- do.call(run_coclustering, c(list(object = s$x), small_grid_args))
  expect_s3_class(res, "coclust_result")
  h <- res$cell_clusters$Hybrid
  ok <- h != "Unstable"
  expect_gt(mean(ok), 0.9)
  # each true cluster maps to a single consensus cluster
  tab <- table(s$truth[ok], droplevels(h[ok]))
  expect_equal(nrow(tab), ncol(tab))
  expect_gt(sum(apply(tab, 1, max)) / sum(tab), 0.98)
})

test_that("louvain also runs", {
  s <- sim_embedding()
  g <- run_clustering_grid(s$x, n_dims = 10, k = 15, resolutions = c(.5, 1),
                           algorithm = "louvain", verbose = FALSE)
  expect_equal(ncol(g$labels), 2)
})
