# Shared helpers for the dev_test scripts.
#
# load_spaniche(code_dir): makes SpaNiche functions available.
#   - a directory containing R/*.R sourced into the global env (used for the
#     untouched baseline snapshot, so its broken NAMESPACE is bypassed), or
#   - "pkg" to devtools::load_all() the package in the repo root (nothing else
#     attached, so missing namespace imports show up as errors).

suppressPackageStartupMessages(library(Seurat))

# Scripts are run from dev_test/ (pixi's workspace root).
repo_root <- normalizePath("..")

load_spaniche <- function(code_dir) {
  if (identical(code_dir, "pkg")) {
    devtools::load_all(repo_root, quiet = TRUE)
  } else {
    # The original code only works with these attached (its NAMESPACE imports none of them).
    suppressPackageStartupMessages({
      library(tidyverse); library(IntNMF); library(Matrix); library(patchwork); library(scales)
      if (requireNamespace("scSeqComm", quietly = TRUE)) library(scSeqComm)
    })
    pointDistance <<- raster::pointDistance
    for (f in list.files(file.path(code_dir, "R"), full.names = TRUE)) source(f)
    load(file.path(code_dir, "data", "CellChatDB.human.rda"), envir = globalenv())
  }
  invisible(TRUE)
}

example_dir <- function() file.path(repo_root, "tutorial", "example", "ST-colon3")

# Run one step, report OK/FAIL with timing, never abort the script.
step <- function(name, expr) {
  t0 <- Sys.time()
  res <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      message("  [warning] ", conditionMessage(w)); invokeRestart("muffleWarning")
    }),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "step_error")
  )
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (inherits(res, "step_error")) {
    cat(sprintf("[FAIL] %-55s %6.1fs  %s\n", name, dt, gsub("\n", " ", res$msg)))
    return(invisible(NULL))
  }
  cat(sprintf("[ OK ] %-55s %6.1fs\n", name, dt))
  invisible(res)
}
