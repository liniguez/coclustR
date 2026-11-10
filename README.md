# coclustR

Co-clustering de consenso para single-cell. Funciona con objetos **SingleCellExperiment**, **Seurat** o con una matriz de embedding (células × PCs).

1. Clustering Leiden/Louvain sobre grafos SNN (estilo Seurat) en una rejilla de `n_dims × k × resolutions`.
2. Las células con el mismo vector de etiquetas en todas las corridas se colapsan en **meta-barcodes**.
3. Se calcula la **matriz de co-clustering** (fracción de corridas en las que dos meta-barcodes caen en el mismo cluster).
4. Los meta-barcodes inestables (core pequeño y poco estable) se marcan y se excluyen.
5. Se construye un `hclust` sobre `1 - co` y se corta con **dynamicTreeCut** (`Hybrid` y `DynamicTree`).
6. Las etiquetas vuelven a las células; las no asignadas quedan como `"Unstable"`.

## Instalación

```r
# Dependencias de Bioconductor
BiocManager::install(c("bluster", "SingleCellExperiment"))
install.packages(c("leidenbase", "igraph", "future.apply", "dynamicTreeCut", "ggplot2"))
# Opcionales
install.packages(c("SeuratObject", "circlize")); BiocManager::install("ComplexHeatmap")

install.packages("coclustR_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Viñetas

El paquete trae dos viñetas con datos simulados (`simulate_sc_counts()`):

```r
vignette("coclustR-sce")     # flujo con SingleCellExperiment (explicación detallada)
vignette("coclustR-seurat")  # flujo con Seurat
browseVignettes("coclustR")
```

## Uso rápido

```r
library(coclustR)
future::plan(future::multisession, workers = 4)   # el paquete NO cambia el plan
options(future.globals.maxSize = 6 * 1024^3)

# SingleCellExperiment
sce <- run_coclustering(sce, reduction = "PCA")
# Seurat
seu <- run_coclustering(seu, reduction = "pca")

table(sce$coclust_Hybrid)
res <- get_coclust_result(sce)          # metadata(sce)$coclustR / seu@misc$coclustR
res
```

### Equivalencia con los parámetros del script original

| Script | Argumento |
|---|---|
| `resolutions`, `neighbours`, `n_dims` | `resolutions`, `k`, `n_dims` |
| `pca.name` | `reduction` |
| `random_seed` | `seed` |
| `clustering_algorithm` | `algorithm = "leiden"/"louvain"` |
| `prune_thresh` | `prune` |
| `cluster_stability_thershold` | `stability_threshold` |
| `Cells < 50 & CoreStab < 0.8` (hard-coded) | `min_core_cells`, `min_core_stability` |
| `hclustMethod` | `linkage` |
| `minModuleSize / minClusterSize = 50` | `min_cluster_size` |
| `deepSplit` | `deep_split` |

### Funciones por pasos

```r
grid <- run_clustering_grid(sce, reduction = "PCA", silhouette = TRUE)
res  <- coclust_from_grid(grid, linkage = "average", min_cluster_size = 30)
# reutilizar la rejilla con otro consenso, sin repetir los clusterings:
sce  <- run_coclustering(sce, grid = grid, linkage = "average")
```

Bloques de bajo nivel: `build_snn_graphs()`, `collapse_metabarcodes()`,
`coclustering_matrix()`, `core_stability()`, `flag_unstable()`,
`consensus_clusters()`.

### Visualización

```r
plot_grid_overview(sce)                      # k × resolución, nº de clusters y modularidad
plot_grid_metrics(sce)                       # modularidad / silhouette / nº clusters vs resolución
plot_stability(sce)                          # diagnóstico del filtro de inestables
plot_clusters_embedding(sce, color_by = c("coclust_Hybrid", "coclust_DynamicTree"))
plot_coclustering_heatmap(sce, n_cells = 1500,
                          stratify_by = c("Sample", "Hybrid"),
                          annotate_by = c("Sample", "Sex"),
                          annotation_colors = list(Sex = c(Male = "green", Female = "purple")))
```

Todas devuelven objetos ggplot (o `Heatmap` de ComplexHeatmap), combinables con patchwork.
