# Tutorial 1a end to end (spaniche_nmf_celltype_and_lr on ST-colon3).
# Usage: Rscript scripts/run_tutorial_1a.R <code_dir|pkg> <out_dir> [extra args as key=value]
#   e.g. n_threads=4
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
setwd(out_dir)  # the package writes intermediate files into the working directory

cat("Seurat", as.character(packageVersion("Seurat")), "| code:", code_dir,
    "| BLAS threads at start:", RhpcBLASctl::blas_get_num_procs(), "\n")

spatial.seu <- Load10X_Spatial(ex)
spatial.seu <- NormalizeData(spatial.seu, verbose = FALSE)
spatial.seu@meta.data$SB <- rownames(spatial.seu@meta.data)
if (exists("get_spatial_coordinates")) {
  spatialdf <- get_spatial_coordinates(spatial.seu, data.dir = ex)
} else {
  spatialdf <- spatial.seu@images$slice1@coordinates
  spatialdf$col <- spatialdf$col / (3^0.5)
}
spot_by_celltype <- read.csv(file.path(ex, "cell2location_map", "res_spotXcelltype.csv"),
                             header = TRUE, row.names = 1, check.names = FALSE)

profile <- nzchar(Sys.getenv("SPANICHE_PROFILE"))
if (profile) Rprof(file.path(out_dir, "Rprof.out"), interval = 0.05)
t0 <- Sys.time()
res.list <- do.call(spaniche_nmf_celltype_and_lr, c(list(
  spatial.seu,
  spatialdf,
  spot_by_celltype,
  LRDB = CellChatDB.human,
  smoothing_type = c("view0", "view1", "view2"),
  distance_thre = c(NA, 2 / (3^0.5), 4 / (3^0.5)),
  LR_dm_type = c("view1", "view2"),
  LR_distance_thre = c(2 / (3^0.5), 4 / (3^0.5)),
  topic_num = 15,
  defined_weight = c(0.6, 0.4),
  lambda_v = 0.5,
  sigma_v = 0.5
), extra))
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (profile) {
  Rprof(NULL)
  pr <- summaryRprof(file.path(out_dir, "Rprof.out"))$by.total
  fns <- c("spaniche_nmf_celltype_and_lr", "get_spot_by_celltype_extended", "prepare_for_LRintegratedmatrix",
           "make_distance_matrix_LR", "generate_LRintegratedmatrix", "moransi_for_LRintegratedmatrix",
           ".moransi_svf", "saveRDS", "write.csv", "nmf.mnnals", "define_good_weight", "nmf_modified", "apply")
  print(pr[intersect(paste0('"', fns, '"'), rownames(pr)), c("total.time", "total.pct")])
}
saveRDS(res.list, "res.list.rds")
cat(sprintf("TOTAL spaniche_nmf_celltype_and_lr: %.1f s\n", elapsed))
cat("BLAS threads at end:", RhpcBLASctl::blas_get_num_procs(), "\n")
