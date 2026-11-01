# Smoke test of the plotting/downstream functions on a finished tutorial 1a run,
# with only SpaNiche loaded (no tidyverse attached), to check the namespace imports.
# Usage: Rscript scripts/smoke_downstream.R <run_dir>
suppressMessages(devtools::load_all("..", quiet = TRUE))
run <- normalizePath(commandArgs(trailingOnly = TRUE)[1])
example <- normalizePath(file.path("..", "tutorial", "example", "ST-colon3"))
stopifnot(!any(grepl("tidyverse|dplyr|stringr|ggplot2|magrittr", search())))
setwd(tempdir())

check <- function(name, expr) {
  r <- tryCatch({ x <- expr; if (inherits(x, "ggplot") || inherits(x, "patchwork")) invisible(ggplot2::ggplot_build(if (inherits(x, "patchwork")) x[[1]] else x)); "OK" },
                error = function(e) paste("FAIL:", conditionMessage(e)))
  cat(sprintf("%-38s %s\n", name, r))
}

res <- readRDS(file.path(run, "res.list.rds"))
fit <- res$nmf_res
k <- ncol(fit$W)
topic_ct <- as.data.frame(fit$H$H1); rownames(topic_ct) <- paste0("topic", 1:k)
spot_topic <- as.data.frame(fit$W); colnames(spot_topic) <- paste0("topic", 1:k)
moransi <- readRDS(file.path(run, "LRintegratedmatrix_with_moransi.rds"))
pos <- read_tissue_positions(example)
coord <- data.frame(col = pos[rownames(spot_topic), "col"] / sqrt(3), row = pos[rownames(spot_topic), "row"], SB = rownames(spot_topic))

check("plot_topic_lr", plot_topic_lr(t(fit$H$H2)))
check("plot_NMF_components", plot_NMF_components(topic_ct, "topic2"))
check("stat_for_interactionwithmoransi", stat_for_interactionwithmoransi(moransi))
check("distribution_2d (point)", distribution_2d(t(spot_topic), "topic1", coord, thre_to_zero = 0))
check("DEinteractions_twoclusters", {
  d <- DEinteractions_twoclusters(spot_topic, t(res$LRintegratedmatrix))
  stopifnot(nrow(d$all.topic.df) > 0); d
})
