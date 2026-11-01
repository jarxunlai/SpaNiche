# SpaNiche 1.2.0

The analysis method is unchanged; results agree with 1.1.0 to floating-point
precision (see "Consistency" below).

## Compatibility

* Works with Seurat v5 (`Assay5`, `VisiumV2`) as well as Seurat v4. Expression
  is read through `LayerData()` / `GetAssayData()` instead of the `@data` /
  `@scale.data` slots; split layers (merged objects) and missing `data` layers
  give clear errors.
* New `get_spatial_coordinates()`: the `spatialdf` input (array `row`, `col / sqrt(3)`)
  from a Seurat object. Seurat v5 Visium images no longer store array
  coordinates, so they are read from the Space Ranger output (`data.dir`).
* New `read_tissue_positions()`: reads `tissue_positions_list.csv`
  (Space Ranger 1.x), `tissue_positions.csv` (Space Ranger >= 2.0) and
  `tissue_positions.parquet` (Visium HD, needs arrow).
* New `assay` argument (default `"Spatial"` if present, otherwise the default
  assay) in `spaniche_nmf_celltype_and_lr()`, `spaniche_nmf_celltype_and_gene()`,
  `spaniche_nmf_comm()` and `prepare_for_LRintegratedmatrix()`; new `image`
  argument in `moransi_for_LRintegratedmatrix()`.
* The package no longer needs `library(tidyverse)` to be attached: dplyr,
  stringr, magrittr and ggplot2 functions are imported explicitly (previously
  an unqualified `filter()` could silently resolve to `stats::filter`).

## Threads

* New `n_threads` argument in the three `spaniche_nmf_*()` functions, default
  `getOption("SpaNiche.n_threads", 1L)`. It sets the number of BLAS/OpenMP
  threads for the call and restores the previous setting on exit. Previously
  OpenBLAS/MKL used every core (up to 32 of 32 measured). `n_threads = NULL`
  keeps the old behaviour.
* In `spaniche_nmf_comm(num_core = N)`, each scSeqComm worker is limited to one
  BLAS thread, so N workers no longer start N x (all cores) threads.

## Speed

* Moran's I of all LR interactions is computed with matrix products (same
  definition and p-values as Seurat's `moransi` method with ape) instead of one
  feature at a time: ~50-70 min -> ~5 s on 4,007 spots x 6,312 interactions.
* Neighbour matrices are sparse and built from a windowed pair search; LR
  integration multiplies each distinct gene once.
* `nmf_modified()`: no explicit Laplacian, `D %*% W` as a vector product,
  `W %*% (H %*% t(H))`, and `A %*% W` reused between objective and update.
* The initial integrative NMF reproduces `IntNMF::nmf.mnnals()` bit for bit
  without building its unused n x n consensus matrices.

## Dependencies

* Removed: raster, reticulate, plyr, pheatmap, tidyverse (meta-package).
* Added: SeuratObject, Matrix, RhpcBLASctl, dplyr, stringr, magrittr, ggplot2.
* scSeqComm (GitLab) and Nebulosa (Bioconductor) moved to Suggests; they are
  only needed by `spaniche_nmf_comm()` and
  `distribution_2d(visualization_type = "density")`.

## Other changes

* `make_distance_matrix()` returns sparse `dgCMatrix` matrices (same values;
  use `as.matrix()` for a dense copy).
* Spots without any neighbour within the distance threshold now give an
  all-zero adjacency row and a warning, instead of a silent row of `NaN`.
* `moransi_for_LRintegratedmatrix()` p-values always follow the ape
  (normal approximation) definition; in 1.1.0 they depended on whether Rfast2
  was installed. Observed Moran's I values are unaffected.
* Added unit tests comparing every optimised step with the original
  implementation.
