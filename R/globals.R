# Column names used with non-standard evaluation (dplyr, ggplot2), the
# magrittr placeholder, and nmf_modified()'s default arguments, declared so
# R CMD check does not report them as undefined globals.
utils::globalVariables(c(
  ".", "LRinteraction2", "MoransI_observed", "MoransI_p.value", "S_intra",
  "Symbol", "Term", "celltype", "cluster", "fit", "log2FC", "p.adj", "p.value",
  "pathway_name", "pct.1", "pct.2", "pval", "pval_log10_neg", "rho_v", "topic",
  "value", "var_value", "view"
))
