# Probe every Seurat-object access in SpaNiche independently, so one failure
# does not hide the next. Usage: Rscript scripts/probe_seurat.R <code_dir|pkg>
source("scripts/common.R")
args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (length(args)) args[1] else "baseline"
load_spaniche(code_dir)

cat("Seurat", as.character(packageVersion("Seurat")),
    "/ SeuratObject", as.character(packageVersion("SeuratObject")),
    "/ code:", code_dir, "\n\n")

seu <- Load10X_Spatial(example_dir())
seu <- NormalizeData(seu, verbose = FALSE)
cat("assay class:", class(seu[["Spatial"]]), " image class:", class(seu@images[[1]]),
    " image names:", names(seu@images), "\n\n")

# Fallback coordinates read straight from Space Ranger output, so later probes
# can still run when the tutorial's own way of getting them fails.
pos <- read.csv(file.path(example_dir(), "spatial", "tissue_positions_list.csv"),
                header = FALSE, row.names = 1,
                col.names = c("barcode", "tissue", "row", "col", "imagerow", "imagecol"))
fallback_df <- pos[colnames(seu), c("row", "col")]
fallback_df$col <- fallback_df$col / sqrt(3)

spatialdf <- step("tutorial: spatial.seu@images$slice1@coordinates", {
  d <- seu@images$slice1@coordinates
  d$col <- d$col / sqrt(3)
  d
})
if (is.null(spatialdf)) spatialdf <- fallback_df

if (exists("get_spatial_coordinates")) {
  step("new: get_spatial_coordinates(seu)", {
    d <- get_spatial_coordinates(seu, data.dir = example_dir())
    stopifnot(isTRUE(all.equal(d[, c("row", "col")], fallback_df[, c("row", "col")])))
    d
  })
}

lr <- step("prepare_for_LRintegratedmatrix(seu, CellChatDB)",
           prepare_for_LRintegratedmatrix(seu, CellChatDB.human))

get_data <- if (exists(".get_layer_data")) {
  function(x) .get_layer_data(x, "data")
} else {
  function(x) x@assays$Spatial@data
}
step("spaniche_nmf_celltype_and_gene: normalized data of HVGs", {
  s <- FindVariableFeatures(seu, nfeatures = 200, verbose = FALSE)
  dim(get_data(s)[VariableFeatures(s), ])
})

step("spaniche_nmf_comm: normalized data", dim(get_data(seu)))

# Moran's I on a small synthetic LR matrix (features x spots).
set.seed(1)
small <- matrix(rexp(30 * ncol(seu)), nrow = 30,
                dimnames = list(paste0("L", 1:30, "->R", 1:30, " (CisL, NisR)"), colnames(seu)))
cc <- data.frame(ligand = paste0("L", 1:30), receptor = paste0("R", 1:30))
step("moransi_for_LRintegratedmatrix (30 features)", {
  m <- moransi_for_LRintegratedmatrix(small, seu, cc)
  stopifnot(nrow(m) == 30, all(c("MoransI_observed", "MoransI_p.value") %in% colnames(m)))
  m
})
