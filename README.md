## SpaNiche: spatial niche analysis to explore colocalization patterns and cellular interactions in spatial transcriptomics data

We propose a computational framework for spatial niche analysis
(SpaNiche) in spatial transcriptomics data to uncover colocalization
patterns and infer potential ligand-receptor interactions. SpaNiche
leverages graph-regularized joint non-negative matrix factorization to
integrate information from cell abundance and ligand-receptor
expression, identifying spatial colocalization patterns among cell types
while providing insights into associated ligand-receptor interactions.
In addition, SpaNiche employs consensus clustering to define “ecotypes”,
enhancing its utility in multi-sample spatial transcriptomics datasets.

<figure>
<img src="./tutorial/example/figure_1_framework.png"
alt="Schematic overview of SpaNiche" />
<figcaption aria-hidden="true">Schematic overview of
SpaNiche</figcaption>
</figure>

> **About this repository.** This is SpaNiche **1.2.0**, a maintenance and
> performance release of the original SpaNiche
> ([SiyuanHuang1/SpaNiche](https://github.com/SiyuanHuang1/SpaNiche), see
> [Citation](#citation)). The method is unchanged and results agree with
> 1.1.0 to floating-point precision; what changed is compatibility (Seurat v5),
> resource use (explicit thread control) and speed. See [NEWS.md](./NEWS.md).

## What's new in 1.2.0

| | 1.1.0 | 1.2.0 |
|---|---|---|
| Seurat v5 objects (`Assay5`, `VisiumV2`) | fails (`no slot of name "data"`, `"coordinates"`, `"scale.data"`) | supported (and Seurat v4) |
| Space Ranger positions | `tissue_positions_list.csv` via the Seurat v4 image | also `tissue_positions.csv` (SR >= 2.0) and `.parquet` (Visium HD) via `read_tissue_positions()` |
| CPU threads | BLAS uses every core (measured peak: 32 of 32) | `n_threads` argument, default 1 (`options(SpaNiche.n_threads = N)`) |
| Tutorial 1a run time (4,007 spots) | 55.5 min | 3.4 min |
| Tutorial 1b run time | 83 min | 10.4 min |
| `library(tidyverse)` required | yes (NAMESPACE imported nothing usable) | no |

Main changes:

- **Seurat v5 support.** Expression is read with `LayerData()` /
  `GetAssayData()` instead of `@data` / `@scale.data`. Because Seurat v5 Visium
  images only keep pixel centroids, the new `get_spatial_coordinates(object,
  data.dir)` returns `spatialdf` (array `row`, `col / sqrt(3)`) from the Space
  Ranger output. All analysis functions get an `assay` argument.
- **Thread control.** The three `spaniche_nmf_*()` functions take
  `n_threads` (BLAS/OpenMP threads for the call, restored afterwards). With
  `spaniche_nmf_comm(num_core = N)` every scSeqComm worker uses one BLAS
  thread. `n_threads = NULL` keeps the previous unrestricted behaviour.
  More threads do not make SpaNiche faster (8 threads: 392 s vs 1 thread:
  349 s on tutorial 1a in an earlier measurement), hence the default of 1.
- **Speed, with identical maths.**
  - Moran's I for all LR interactions at once with matrix products (same
    weights, statistic and p-values as Seurat's `moransi` method with ape),
    instead of one feature at a time: 50 min -> 5 s.
  - Sparse neighbour matrices; each distinct gene multiplied once in the LR
    integration.
  - `nmf_modified()` without an explicit n x n Laplacian, with `D %*% W` as a
    vector product, `W %*% (H %*% t(H))`, and `A %*% W` reused.
  - The initial integrative NMF reproduces `IntNMF::nmf.mnnals()` bit for bit
    without its unused n x n consensus matrices.
- **Packaging.** Explicit imports of dplyr/stringr/magrittr/ggplot2; raster,
  reticulate, plyr, pheatmap removed; scSeqComm and Nebulosa are optional;
  unit tests compare each optimised step with the original code.

## Benchmark

**Data:** tutorial example ST-colon3 (10x Visium, 4,007 spots, 36,601 genes,
7,704 LR integration features, 15 topics, parameters as in the tutorials).
**Machine:** AMD Ryzen 9 9950X (16 cores / 32 threads), 30 GB RAM, Linux
(WSL2, kernel 6.6), OpenBLAS 0.3.34. **Software:** R 4.3.3 + Seurat 4.4.0 /
SeuratObject 4.1.4 ("Seurat v4"), R 4.4.3 + Seurat 5.5.1 / SeuratObject 5.4.0
("Seurat v5"), IntNMF 1.3.0, scSeqComm 2.0.0. Time is the wall time of the
analysis function call; cores are the CPU cores busy during the call, sampled
every second.

### Run time and CPU use

| Pipeline | 1.1.0, Seurat v4 (default threads) | 1.2.0, Seurat v4 | 1.2.0, Seurat v5 | Speed-up |
|---|---|---|---|---|
| Tutorial 1a `spaniche_nmf_celltype_and_lr()` | 3,327 s; peak 32.2 / mean 5.9 cores | 207 s; 1.0 core | 205 s; 1.0 core | 16x |
| `spaniche_nmf_celltype_and_gene()` | 564 s; peak 28.7 / mean 10.8 cores | 220 s; 1.0 core | 223 s; 1.0 core | 2.6x |
| Tutorial 1b `spaniche_nmf_comm()` | 4,987 s; peak 32.0 / mean 4.4 cores | 623 s; 1.0 core | 603 s; 1.1 cores | 8.0x |
| Any of the above with Seurat v5 | error | - | runs | - |

In CPU time (cores x wall time), tutorial 1a goes from ~19,600 to ~205
core-seconds (~95x less).

### Where the time went (tutorial 1a, Seurat v4)

| Stage | 1.1.0 | 1.2.0 |
|---|---|---|
| Spatial smoothing + LR integrated matrix | 10 s | 7 s |
| Moran's I (6,312 LR features) | 3,008 s | 5 s |
| Integrative NMF (initialisation + spatially regularised NMF) | 305 s | 194 s |
| Tutorial 1b only: scSeqComm + GO per topic | 429 s | 405 s |

The 1.2.0 runs were measured one at a time on an otherwise idle machine
(`dev_test/scripts/bench.sh`). The 1.1.0 runs took 55-83 min each and were
measured once; parts of them overlapped with other (single-core) jobs, which
mainly slowed the single-threaded Moran's I stage (50 min in tutorial 1a,
70 min in tutorial 1b), so their numbers are an approximation.

### Consistency of results

Same input and parameters, 1.2.0 vs 1.1.0 under Seurat v4:

| Output | Difference |
|---|---|
| LR integrated matrix (4,007 x 7,704) | max abs 3.6e-15 |
| Moran's I observed / p-value | max abs 1.1e-16 / 1.1e-12 |
| NMF `W`, `H` (tutorial 1a, 1b) | max relative 5e-12; same iterations (199), objective, clusters |
| NMF `W`, `H` (`celltype_and_gene`) | max relative 2e-13; same iterations, objective, clusters |
| Topic-LR table (1b) | same 3,739 rows; values within 3.4e-10 |
| GO enrichment (1b) | identical terms, order and p-values |

- 1.2.0 gives bit-identical NMF and topic-LR results under Seurat v4 and v5.
- Differences are floating-point rounding from a different order of
  summation (sparse products, matrix associativity).
- In the Moran's I table, interactions whose statistics are tied to 1e-19
  may swap rows; ranks only move among non-significant ones (p > 0.27) and
  the `moransi.spatially.variable` flags are identical.
- Against the authors' published ST-colon3 output: LR matrix within 5e-15,
  observed Moran's I within 5e-14.
- GO results also depend on the installed annotation (GO.db / org.Hs.eg.db)
  version, independent of SpaNiche.

### Reproducing

`dev_test/` contains a pixi project with both Seurat versions and the
scripts used above (see [dev_test/README.md](./dev_test/README.md)):

    cd dev_test
    pixi install --all
    pixi run -e v5 install-extra && pixi run -e v4 install-extra
    pixi run -e v5 test                 # unit tests (new vs original code)
    bash scripts/bench.sh               # benchmark of the current package

## Installation

SpaNiche has been tested on both Windows (R version 4.2.1) and
Linux-based systems (R version 4.2.3, 4.3.3 and 4.4.3), with Seurat v4
(4.4.0) and Seurat v5 (5.5.1), and should be compatible with other R
environments.

You can install SpaNiche 1.2.0 from GitHub using:

    remotes::install_github("jarxunlai/SpaNiche@v1.2.0")

or from the source package attached to the
[v1.2.0 release](https://github.com/jarxunlai/SpaNiche/releases/tag/v1.2.0):

    install.packages("SpaNiche_1.2.0.tar.gz", repos = NULL, type = "source")

`spaniche_nmf_comm()` (LR-induced pathway inference) additionally needs
scSeqComm, and `distribution_2d(visualization_type = "density")` needs
Nebulosa:

    remotes::install_gitlab("sysbiobig/scseqcomm")
    BiocManager::install("Nebulosa")

To control the number of CPU threads (default 1):

    options(SpaNiche.n_threads = 4)   # or n_threads = 4 in each call

## Application examples

-   [1a. Graph-regularized integrated NMF in single-sample spatial transcriptomics data using SpaNiche](./tutorial/1a.singlesample_iNMF.md)

This analysis typically reveals colocalization patterns of different cell type combinations within a spatial transcriptomics slide, along with their corresponding ligand-receptor interaction features. Because it jointly factorizes cell type abundance and ligand-receptor matrices, it often yields better correspondence between patterns and interactions. Moreover, the workflow is more streamlined compared to two-step approaches (e.g., spatial clustering followed by CellphoneDB/cellchat analysis).

-   [1b. Graph-regularized integrated NMF in single-sample ST data, and LR-induced pathway inference](./tutorial/1b.singlesample_iNMF_LRpathway.md)

We add inference of pathway activation induced by high expression of specific ligand-receptor interactions based on **1a**.

-   [2a. Ecotype construction using multi-sample ST data](./tutorial/2a.multisample_ecotype.md)

This section provides a detailed tutorial on ecotype construction using multi-sample spatial transcriptomics datasets, which represents one of the core components of the SpaNiche framework.

-   [2b. Representative cell types and LR interactions of each ecotype](./tutorial/2b.ecotype_related_celltype_and_LR.md)

After defining ecotypes, the next important step is to characterize the biological features associated with each ecotype. In SpaNiche, two major aspects are typically investigated: 1. Representative cell types enriched within each ecotype; 2. Representative ligand-receptor (LR) interactions.

## Citation

If you have used SpaNiche or referenced some concepts from our paper (such as topic, niche, ecotype, colocalization, etc.), please cite our paper:

Huang, S., Ran, Q., Tang, J., Wang, X., Xi, J., Ma, S., & Xi, R. (2026). SpaNiche: spatial niche analysis to explore colocalization patterns and cellular interactions in spatial transcriptomics data. Genome Biology.
