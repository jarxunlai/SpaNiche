#' make_distance_matrix
#'
#' @description
#' The function calculates Euclidean distances between spatial spots and converts them into adjacency matrices based on specified distance thresholds.
#'
#' @param spatialdf a data frame or matrix that contains two columns, "col" and "row". Row names should  correspond to spot barcodes.
#' @param dm_type This parameter can be one or more values from c('view1', 'view2'), representing the type of adjacency matrix to be returned. For example, you can choose c('view1') or c('view1', 'view2').
#' @param distance_thre This parameter represents distance, defining the central distance from center to it's first/second-order neighbors. Consistent with dm_type, you can define c(c1) or c(c1, c2), and c1/c2 is distance value.
#' @param digits integer indicating the number of decimal places.
#'
#' @return A named list of adjacency matrices corresponding to each dm_type.
#' Each matrix is normalized so that rows sum to 1 (row-normalization matrix).
#' Matrices are sparse (\code{dgCMatrix}); use \code{as.matrix()} for a dense copy.
#'
#' @import tidyverse
#' @export
make_distance_matrix = function(
    spatialdf,
    dm_type = c("view1", "view2"),
    distance_thre,
    digits = 2
){
  distance_thre = round(distance_thre,digits)

  if (!all(c("col", "row") %in% colnames(spatialdf))) {
    stop("spatialdf must contain columns named 'col' and 'row'")
  }
  spatialdf=spatialdf[,c("col","row")]
  spatialdf=as.matrix(spatialdf)
  if (length(dm_type) != length(distance_thre)) {
    stop("Length of 'distance_thre' must match 'dm_type'")
  }
  dm_type=sort(dm_type)
  distance_thre=sort(distance_thre)
  names(distance_thre) = dm_type

  # Only pairs within the largest threshold can be neighbours, so distances are
  # computed for those pairs only and kept sparse.
  pairs = .neighbor_pairs(spatialdf, max(distance_thre), digits)

  res_list=as.list(1:length(dm_type))
  names(res_list)=dm_type

  for (view in intersect(c("view1", "view2"), dm_type)) {
    thre = distance_thre[view]
    # Same two-step recoding as the dense implementation: <= thre -> 1, then > thre -> 0.
    v = pairs$d
    v[v <= thre] = 1
    v[v > thre] = 0
    keep = pairs$d <= thre & v != 0
    res_list[[view]] = .row_normalized_adjacency(
      pairs$i[keep], pairs$j[keep], nrow(spatialdf), rownames(spatialdf)
    )
  }

  res_list
}
