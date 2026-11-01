# Packages SpaNiche needs that are not packaged on conda-forge/bioconda.
# scSeqComm's other dependencies come from conda (see pixi.toml).
# timeout: the quickSparseM source archive is ~400 MB
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8, timeout = 3600)
if (!requireNamespace("IntNMF", quietly = TRUE)) install.packages("IntNMF")
if (!requireNamespace("add2ggplot", quietly = TRUE)) install.packages("add2ggplot")
if (!requireNamespace("quickSparseM", quietly = TRUE)) {
  remotes::install_gitlab("sysbiobig/quicksparsem", upgrade = "never", dependencies = FALSE)
}
if (!requireNamespace("scSeqComm", quietly = TRUE)) {
  remotes::install_gitlab("sysbiobig/scseqcomm", upgrade = "never", dependencies = FALSE)
}
for (p in c("IntNMF", "add2ggplot", "quickSparseM", "scSeqComm")) cat(p, as.character(packageVersion(p)), "\n")
