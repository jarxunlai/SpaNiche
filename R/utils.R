# Internal helpers shared across SpaNiche.

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------------
# Seurat v4 / v5 compatible access to an expression matrix.
# Seurat v5 stores matrices as layers of an Assay5 object (no @data slot);
# Seurat v4 stores them as slots of an Assay object.
# ------------------------------------------------------------------
# Assay used when none is given: "Spatial" (what SpaNiche has always used)
# if present, otherwise the default assay (e.g. "Spatial.008um" for Visium HD).
.resolve_assay <- function(object, assay = NULL) {
  if (!is.null(assay)) return(assay)
  if ("Spatial" %in% SeuratObject::Assays(object)) return("Spatial")
  SeuratObject::DefaultAssay(object)
}

.get_layer_data <- function(object, layer = "data", assay = NULL) {
  assay <- .resolve_assay(object, assay)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    assay_obj <- object[[assay]]
    if (inherits(assay_obj, "StdAssay")) {
      # Looked up at run time: these functions only exist in SeuratObject >= 5.
      found <- getExportedValue("SeuratObject", "Layers")(assay_obj, search = layer)
      if (length(found) == 0) {
        stop("Layer '", layer, "' not found in assay '", assay, "'.",
             if (layer == "data") " Run NormalizeData() first." else "", call. = FALSE)
      }
      if (length(found) > 1) {
        stop("Assay '", assay, "' has split layers (", paste(found, collapse = ", "),
             "). Run JoinLayers() first, or analyse one sample at a time.", call. = FALSE)
      }
    }
    return(getExportedValue("SeuratObject", "LayerData")(object, assay = assay, layer = layer))
  }
  SeuratObject::GetAssayData(object, assay = assay, slot = layer)
}

# ------------------------------------------------------------------
# Thread control for the BLAS/OpenMP libraries R is linked against.
# Matrix products in SpaNiche go through BLAS, which by default (OpenBLAS,
# MKL) starts one thread per core. Returns the previous settings so the
# caller can restore them with .restore_threads() via on.exit().
# ------------------------------------------------------------------
.set_threads <- function(n_threads) {
  if (is.null(n_threads)) return(NULL)
  n_threads <- as.integer(n_threads)
  if (length(n_threads) != 1 || is.na(n_threads) || n_threads < 1) {
    stop("n_threads must be NULL or a single positive integer.", call. = FALSE)
  }
  old <- list(blas = RhpcBLASctl::blas_get_num_procs(),
              omp = RhpcBLASctl::omp_get_max_threads())
  RhpcBLASctl::blas_set_num_threads(n_threads)
  RhpcBLASctl::omp_set_num_threads(n_threads)
  old
}

.restore_threads <- function(old) {
  if (is.null(old)) return(invisible(NULL))
  if (!is.na(old$blas)) RhpcBLASctl::blas_set_num_threads(old$blas)
  if (!is.na(old$omp)) RhpcBLASctl::omp_set_num_threads(old$omp)
  invisible(NULL)
}

# ------------------------------------------------------------------
# All pairs of points (including each point with itself) whose Euclidean
# distance, rounded to `digits`, is <= max_dist. Points are processed in
# blocks against an x-sorted window, so memory stays O(block * window)
# instead of O(n^2).
# Returns a data.frame with 1-based indices i, j and the rounded distance d.
# ------------------------------------------------------------------
.neighbor_pairs <- function(coords, max_dist, digits, block = 1024L) {
  x <- as.numeric(coords[, 1])
  y <- as.numeric(coords[, 2])
  n <- length(x)
  ord <- order(x)
  xs <- x[ord]
  # A pair whose true distance is slightly above max_dist can still round down to it.
  margin <- max_dist + 10^(-digits)

  res <- vector("list", ceiling(n / block))
  for (b in seq_along(res)) {
    rows <- ord[((b - 1L) * block + 1L):min(n, b * block)]
    lo <- findInterval(min(x[rows]) - margin, xs, left.open = TRUE) + 1L
    hi <- findInterval(max(x[rows]) + margin, xs)
    cand <- ord[lo:hi]
    d <- round(sqrt(outer(x[rows], x[cand], "-")^2 + outer(y[rows], y[cand], "-")^2), digits)
    hit <- which(d <= max_dist, arr.ind = TRUE)
    res[[b]] <- data.frame(i = rows[hit[, 1]], j = cand[hit[, 2]], d = d[hit])
  }
  do.call(rbind, res)
}

# Sparse, row-normalised (rows sum to 1) adjacency matrix from pair indices.
.row_normalized_adjacency <- function(i, j, n, names) {
  adj <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n, n),
                              dimnames = list(names, names))
  rs <- Matrix::rowSums(adj)
  if (any(rs == 0)) {
    # The former dense implementation silently produced rows of NaN here.
    warning(sum(rs == 0), " spot(s) have no neighbours within the distance threshold; ",
            "their rows of the adjacency matrix are all zero.", call. = FALSE)
  }
  adj <- Matrix::Diagonal(x = 1 / rs) %*% adj
  dimnames(adj) <- list(names, names)
  methods::as(adj, "CsparseMatrix")
}

# Evaluate `expr` with BLAS/OpenMP limited to one thread in the worker
# processes it starts (forked workers inherit the current BLAS setting,
# fresh R workers read the environment variables), so that `n_workers`
# workers do not each spawn one BLAS thread per core.
.with_single_thread_workers <- function(expr, n_workers) {
  if (n_workers <= 1) return(expr)
  vars <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
  old <- Sys.getenv(vars, unset = NA)
  do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(vars)), vars)))
  on.exit({
    for (v in vars) {
      if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, stats::setNames(list(old[[v]]), v))
    }
  }, add = TRUE)
  old_threads <- .set_threads(1L)
  on.exit(.restore_threads(old_threads), add = TRUE)
  expr
}
