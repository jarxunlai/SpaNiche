#' Calculate Moran's I for spatial ligand-receptor integrated matrix
#'
#' @param LRintegratedmatrix A numeric matrix of integrated LR interaction
#'   scores. Rows correspond to LR interactions and columns correspond to
#'   spatial spots. Row names represent interaction identifiers.
#'   Underscores ('_') in row names will be replaced by hyphens ('-') to ensure
#'   consistency with Seurat's feature naming conventions.
#' @param spatial.seu A Seurat object containing spatial transcriptomics data.
#'   This object must include spatial image information (e.g. Visium) and
#'   spot coordinates. Column names of LRintegratedmatrix must match
#'   the spot names in this object. Works with Seurat v4 and v5 objects.
#' @param cc_interaction A data frame containing ligand–receptor interaction
#'   annotations.
#' @param image Name of the spatial image whose tissue coordinates are used.
#'   Defaults to the first image of \code{spatial.seu}.
#'
#' @details Moran's I is computed with the same definition as
#'   \code{Seurat::FindSpatiallyVariableFeatures(selection.method = "moransi")}
#'   (inverse squared distance weights, normal-approximation p-values as in
#'   \code{ape::Moran.I}), but for all interactions at once with matrix
#'   products instead of one feature at a time.
#'
#' @return A data frame containing Moran's I statistics for spatially variable
#'   ligand–receptor interactions, merged with ligand–receptor annotation
#'   information. The output is sorted by decreasing observed Moran's I.
#' @export
moransi_for_LRintegratedmatrix=function(
    LRintegratedmatrix,
    spatial.seu,
    cc_interaction,
    image = NULL
){
  print("calculating moransi for LR integrated matrix:")

  ### Calculate Moran's I
  LRintegratedmatrix=LRintegratedmatrix[rowSums(LRintegratedmatrix) > 0,]
  # Replace underscores with hyphens to match Seurat's feature naming convention
  rownames(LRintegratedmatrix)=rownames(LRintegratedmatrix) %>% str_replace_all("_","-")
  # if (!identical(colnames(LRintegratedmatrix), colnames(spatial.seu))) {
  #   stop("Column names or order in LRintegratedmatrix are not identical to spatial.seu!")
  # }
  # Ensure that all spatial spots are present in the LR matrix
  stopifnot(all(colnames(spatial.seu) %in% colnames(LRintegratedmatrix)))
  # Reorder LR matrix columns to match the spatial Seurat object
  LRintegratedmatrix=LRintegratedmatrix[,colnames(spatial.seu)]
  # Tissue coordinates as used by Seurat::FindSpatiallyVariableFeatures()
  image = image %||% SeuratObject::Images(spatial.seu)[1]
  tc = SeuratObject::GetTissueCoordinates(spatial.seu[[image]])
  tc = as.matrix(tc[, vapply(tc, is.numeric, logical(1)), drop = FALSE][, 1:2])
  stopifnot(all(colnames(LRintegratedmatrix) %in% rownames(tc)))
  tc = tc[colnames(LRintegratedmatrix), , drop = FALSE]

  # Identify spatially variable LR interactions using Moran's I
  moransi_output_df = .moransi_svf(as.matrix(LRintegratedmatrix), tc)
  moransi_output_df = moransi_output_df %>% arrange(desc(MoransI_observed))
  moransi_output_df$interaction_name_detailed=rownames(moransi_output_df)

  ###
  cc_interaction$interaction_name_new = paste(
    cc_interaction$ligand,cc_interaction$receptor,sep = "->"
    ) %>% str_replace_all("_","-")

  moransi_output_df$interaction_name_new=str_replace(moransi_output_df$interaction_name_detailed," \\(.*","")
  # Ensure consistent ligand->receptor direction
  moransi_output_df$interaction_name_new=moransi_output_df$interaction_name_new %>%
    sapply(
      function(x){
        if(str_detect(x,"<-")){
          tmp1=str_split(x,"<-")[[1]][1]
          tmp2=str_split(x,"<-")[[1]][2]
          tmpnew=paste(tmp2,tmp1,sep = "->")
        }else{
          tmpnew=x
        }
        as.character(tmpnew)
      }
    )

  ### Merge Moran's I results with LR annotation table
  moransi_output_df = cc_interaction %>% inner_join(moransi_output_df,by = "interaction_name_new")
  print("Finished.")
  print("------------------------------------")
  return(moransi_output_df)
}

# Moran's I for every row (feature) of `data` (features x spots), reproducing
# Seurat::FindSpatiallyVariableFeatures(selection.method = "moransi") with ape:
# weights 1/d^2 (zero diagonal), row-normalised; features with zero variance
# are skipped (NA); ranking by p-value then by -|I|; top `nfeatures` flagged.
.moransi_svf = function(data, pos, nfeatures = 2000, chunk = 1000) {
  n = ncol(data)
  w = 1 / as.matrix(dist(pos))^2
  diag(w) = 0
  rs = rowSums(w)
  rs[rs == 0] = 1
  w = w / rs
  s = sum(w)
  S1 = 0.5 * sum((w + t(w))^2)
  S2 = sum((rowSums(w) + colSums(w))^2)
  ei = -1 / (n - 1)

  observed = p.value = rep(NA_real_, nrow(data))
  keep = which(apply(data, 1, var) != 0)
  for (idx in split(keep, ceiling(seq_along(keep) / chunk))) {
    y = data[idx, , drop = FALSE]
    y = y - rowMeans(y)
    v = rowSums(y^2)
    obs = (n / s) * (rowSums(tcrossprod(y, w) * y) / v)
    k = (rowSums(y^4) / n) / (v / n)^2
    sdi = sqrt((n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * s^2) -
                  k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s^2)) /
                 ((n - 1) * (n - 2) * (n - 3) * s^2) - 1 / ((n - 1)^2))
    pv = stats::pnorm(obs, mean = ei, sd = sdi)
    pv = ifelse(obs <= ei, 2 * pv, 2 * (1 - pv))
    observed[idx] = obs
    p.value[idx] = pv
  }

  out = data.frame(MoransI_observed = observed, MoransI_p.value = p.value,
                   row.names = rownames(data))
  ranked = keep[order(p.value[keep], -abs(observed[keep]))]
  out$moransi.spatially.variable = NA
  out$moransi.spatially.variable.rank = NA_integer_
  out$moransi.spatially.variable[ranked] = seq_along(ranked) <= nfeatures
  out$moransi.spatially.variable.rank[ranked] = seq_along(ranked)
  out
}
