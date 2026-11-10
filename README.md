# coclustR

Consensus co-clustering for single-cell data. Works with **SingleCellExperiment** and **Seurat** objects, or with a plain embedding matrix (cells × PCs).

Choosing a single combination of number of PCs, number of neighbours and resolution is arbitrary: a different combination gives different clusters. `coclustR` clusters the cells with **every** combination of a parameter grid and keeps what is reproducible:

1. Leiden/Louvain clustering on Seurat-style SNN graphs over a grid of `n_dims × k × resolutions`.
2. Cells with the same label vector across all runs are collapsed into **meta-barcodes**.
3. A **co-clustering matrix** is computed (fraction of runs in which two meta-barcodes fall in the same cluster).
4. Unstable meta-barcodes (small, weakly co-clustered core) are flagged and excluded.
5. An `hclust` tree is built on `1 - co-clustering` and cut with **dynamicTreeCut** (`Hybrid` and `DynamicTree` methods).
6. Labels are mapped back to cells; unassigned cells are labelled `"Unstable"`.

## Installation

```r
install.packages("remotes")
remotes::install_github("liniguez/coclustR",
                        dependencies = TRUE,     # also installs Seurat, ComplexHeatmap, ...
                        build_vignettes = TRUE)
```

Bioconductor dependencies (`bluster`, `SingleCellExperiment`, `ComplexHeatmap`) are resolved automatically. If you do not have Seurat and scater installed, drop `build_vignettes = TRUE`.

## Vignettes

Two vignettes with simulated data (`simulate_sc_counts()`):

```r
vignette("coclustR-sce")     # SingleCellExperiment workflow (detailed explanation)
vignette("coclustR-seurat")  # Seurat workflow
browseVignettes("coclustR")
```

## Quick start

```r
library(coclustR)
future::plan(future::multisession, workers = 4)   # the package never changes the plan
options(future.globals.maxSize = 6 * 1024^3)

# SingleCellExperiment
sce <- run_coclustering(sce, reduction = "PCA")
# Seurat
seu <- run_coclustering(seu, reduction = "pca")

table(sce$coclust_Hybrid)
res <- get_coclust_result(sce)          # metadata(sce)$coclustR / seu@misc$coclustR
res
```

### Main parameters

| Argument | Default | Meaning |
|---|---|---|
| `reduction` | `"PCA"` / `"pca"` | Embedding used to build the graphs |
| `n_dims`, `k`, `resolutions` | `c(15,30,50,75,100)`, `seq(10,40,5)`, `seq(0.1,1.4,0.1)` | Parameter grid (490 runs) |
| `algorithm` | `"leiden"` | `"leiden"` or `"louvain"` |
| `prune` | `1/15` | SNN pruning threshold (as Seurat's `prune.SNN`) |
| `seed` | `11` | Results are reproducible for any `future` plan |
| `stability_threshold` | `0.5` | Co-clustering needed to belong to a meta-barcode's core |
| `min_core_cells`, `min_core_stability` | `50`, `0.8` | Meta-barcodes below both are flagged as unstable |
| `linkage` | `"complete"` | Linkage of the consensus tree |
| `min_cluster_size` | `50` | Minimum module size for dynamicTreeCut, **in meta-barcodes** |
| `deep_split` | `FALSE` | Sensitivity of the tree cut |
| `rescue_min_cells` | `50` | Very stable populations collapsed into few meta-barcodes become clusters if they have at least this many cells |

### Step by step

```r
grid <- run_clustering_grid(sce, reduction = "PCA", silhouette = TRUE)   # slow part
res  <- coclust_from_grid(grid, linkage = "average", min_cluster_size = 30)
# reuse the grid with another consensus, without re-running the clusterings:
sce  <- run_coclustering(sce, grid = grid, linkage = "average")
```

Low-level building blocks: `build_snn_graphs()`, `collapse_metabarcodes()`,
`coclustering_matrix()`, `core_stability()`, `flag_unstable()`,
`consensus_clusters()`.

### Visualisation

```r
plot_grid_overview(sce)                      # k × resolution, number of clusters and modularity
plot_grid_metrics(sce)                       # modularity / silhouette / n clusters vs resolution
plot_stability(sce)                          # diagnostic of the unstable filter
plot_clusters_embedding(sce, color_by = c("coclust_Hybrid", "coclust_DynamicTree"))
plot_coclustering_heatmap(sce, n_cells = 1500,
                          stratify_by = c("Sample", "Hybrid"),
                          annotate_by = c("Sample", "Sex"),
                          annotation_colors = list(Sex = c(Male = "green", Female = "purple")))
```

All return ggplot objects (or a ComplexHeatmap `Heatmap`) that can be combined with patchwork.
