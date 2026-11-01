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

## Installation

SpaNiche has been tested on both Windows (R version 4.2.1) and
Linux-based systems (R version 4.2.3, 4.3.3 and 4.4.3), with Seurat v4
(4.4.0) and Seurat v5 (5.5.1), and should be compatible with other R
environments.

You can install SpaNiche from GitHub using:

    devtools::install_github("SiyuanHuang1/SpaNiche")

`spaniche_nmf_comm()` (LR-induced pathway inference) additionally needs
scSeqComm, and `distribution_2d(visualization_type = "density")` needs
Nebulosa:

    remotes::install_gitlab("sysbiobig/scseqcomm")
    BiocManager::install("Nebulosa")

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
