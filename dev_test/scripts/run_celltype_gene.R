# spaniche_nmf_celltype_and_gene on ST-colon3 (cell types + 2000 HVGs).
# Usage: Rscript scripts/run_celltype_gene.R <code_dir|pkg> <out_dir> [key=value ...]
source("scripts/common.R")
args <- commandArgs(trailingOnly = TRUE)
code_dir <- args[1]
out_dir <- args[2]
extra <- if (length(args) > 2) {
  kv <- strsplit(args[-(1:2)], "=", fixed = TRUE)
  setNames(lapply(kv, function(x) type.convert(x[2], as.is = TRUE)), vapply(kv, `[`, "", 1))
} else list()
load_spaniche(code_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
ex <- example_dir()
setwd(out_dir)

spatial.seu <- Load10X_Spatial(ex)
spatialdf <- if (exists("get_spatial_coordinates")) {
  get_spatial_coordinates(spatial.seu, data.dir = ex)
} else {
  d <- spatial.seu@images$slice1@coordinates; d$col <- d$col / sqrt(3); d
}
spot_by_celltype <- read.csv(file.path(ex, "cell2location_map", "res_spotXcelltype.csv"),
                             header = TRUE, row.names = 1, check.names = FALSE)
t0 <- Sys.time()
res <- do.call(spaniche_nmf_celltype_and_gene, c(list(
  spatial.seu, spatialdf, spot_by_celltype,
  topic_num = 15, defined_weight = c(0.6, 0.4), lambda_v = 0.5, sigma_v = 0.5
), extra))
cat(sprintf("TOTAL spaniche_nmf_celltype_and_gene: %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
