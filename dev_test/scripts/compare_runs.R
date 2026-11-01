# Compare intermediate outputs of a run against a reference.
# Usage: Rscript scripts/compare_runs.R <run_dir> <ref_dir|author>
#   author: the author's published ST-colon3 output (data/baidu/ST-colon3),
#           whose interaction names use the old "loop" wording.
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
run <- args[1]
ref <- args[2]
author <- identical(ref, "author")
ref_dir <- if (author) "data/baidu/ST-colon3" else ref

rename_loops <- function(x) {
  x <- sub("\\(loop0\\+loop1\\)", "(view1)", x)
  sub("\\(loop0\\+loop2\\)", "(view2)", x)
}
report <- function(name, a, b) {
  d <- max(abs(a - b))
  rel <- d / max(abs(b))
  cat(sprintf("%-42s max|diff| = %.3g  (relative %.3g)\n", name, d, rel))
}

# 1. LR integrated matrix (saved before variance filtering)
a <- readRDS(file.path(run, "LRintegratedmatrix.rds"))
b <- readRDS(file.path(ref_dir, "LRintegratedmatrix.rds"))
if (author) colnames(b) <- rename_loops(colnames(b))
cat("LRintegratedmatrix dims:", dim(a), "vs", dim(b),
    "| same spots:", setequal(rownames(a), rownames(b)),
    "| same interactions:", setequal(colnames(a), colnames(b)), "\n")
common <- intersect(colnames(a), colnames(b))
report("LRintegratedmatrix", a[rownames(b), common], b[, common])

# 2. Moran's I
ma <- readRDS(file.path(run, "LRintegratedmatrix_with_moransi.rds"))
mb <- readRDS(file.path(ref_dir, "LRintegratedmatrix_with_moransi.rds"))
if (author) mb$interaction_name_detailed <- rename_loops(mb$interaction_name_detailed)
key <- function(m) paste(m$interaction_name, m$interaction_name_detailed)
cat("Moran rows:", nrow(ma), "vs", nrow(mb), "| same set:", setequal(key(ma), key(mb)),
    "| same columns:", identical(colnames(ma), colnames(mb)), "\n")
mb <- mb[match(key(ma), key(mb)), ]
report("MoransI_observed", ma$MoransI_observed, mb$MoransI_observed)
if (!author) {
  report("MoransI_p.value", ma$MoransI_p.value, mb$MoransI_p.value)
  cat("identical row order:", identical(key(ma), key(readRDS(file.path(ref_dir, "LRintegratedmatrix_with_moransi.rds")))),
      "| same rank:", identical(ma$moransi.spatially.variable.rank, mb$moransi.spatially.variable.rank), "\n")
}

# 3. NMF result (only meaningful against a run with identical parameters)
if (!author) {
  fa <- readRDS(file.path(run, "nmf_modified_fit.rds"))
  fb <- readRDS(file.path(ref_dir, "nmf_modified_fit.rds"))
  report("NMF W", fa$W, fb$W)
  report("NMF H1 (cell types)", fa$H$H1, fb$H$H1)
  report("NMF H2 (LR)", fa$H$H2, fb$H$H2)
  cat("iterations:", nrow(fa$convergence), "vs", nrow(fb$convergence),
      "| final f.WH:", tail(fa$convergence$f.WH, 1), "vs", tail(fb$convergence$f.WH, 1),
      "| identical clusters:", identical(fa$clusters, fb$clusters), "\n")
  ra <- readRDS(file.path(run, "res.list.rds"))
  rb <- readRDS(file.path(ref_dir, "res.list.rds"))
  cat("filtered LR matrix identical columns:", identical(colnames(ra$LRintegratedmatrix), colnames(rb$LRintegratedmatrix)), "\n")
}
