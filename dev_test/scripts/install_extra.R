# Packages SpaNiche needs that are not packaged on conda-forge/bioconda.
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)
if (!requireNamespace("IntNMF", quietly = TRUE)) install.packages("IntNMF")
if (!requireNamespace("scSeqComm", quietly = TRUE)) {
  remotes::install_github("sysbiobig/scSeqComm", upgrade = "never", dependencies = TRUE)
}
for (p in c("IntNMF", "scSeqComm")) cat(p, as.character(packageVersion(p)), "\n")
