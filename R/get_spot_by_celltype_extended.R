#' get_spot_by_celltype_extended
#'
#' @description
#' Perform multi-scale spatial smoothing of spot-level cell type abundance
#' for 10x Visium spatial transcriptomics data. Depending on the selected
#' smoothing range (view0 / view1 / view2), the function aggregates cell type
#' information from spatial neighbors and concatenates results into an
#' extended feature matrix.
#' @param spatialdf A data frame or matrix with two columns: \code{col} and
#'   \code{row}, representing spatial coordinates of spots. Row names must be
#'   spot barcodes and consistent with \code{spot_by_celltype}.
#' @param spot_by_celltype A numeric matrix or data frame with spots as rows
#'   and cell types as columns, representing cell type abundance or proportion.
#'   Row names must match those of \code{spatialdf}.
#' @param smoothing_type Character vector specifying the spatial smoothing
#'   levels to include. Supported values are:
#'   \itemize{
#'     \item \code{"view0"}: no smoothing (original matrix)
#'     \item \code{"view1"}: first-order spatial neighbors
#'     \item \code{"view2"}: second-order spatial neighbors
#'   }
#'   Valid combinations are \code{c("view0")},
#'   \code{c("view0","view1")}, and
#'   \code{c("view0","view1","view2")}.
#' @param distance_thre This parameter represents distance, defining the central distance from center to loop1 (view1), and from center to loop2 (view2). For view0, the distance is ignored and should be set to NA. The length of the vector is consistent with the smoothing_type parameter. For example, A typical choice for Visium data is: distance_thre = c(NA,round(2/(3^0.5),2),round(4/(3^0.5),2)).
#' @param digits integer indicating the number of decimal places.
#'
#' @return
#' A numeric matrix with the same number of rows (spots) as the input.
#' Columns correspond to concatenated cell type features from different
#' spatial views. Column names are suffixed with \code{_view1} or
#' \code{_view2} when applicable.
#' @export
get_spot_by_celltype_extended=function(
    spatialdf,
    spot_by_celltype,
    smoothing_type = c('view0', 'view1', 'view2'),
    distance_thre = c(NA,2/(3^0.5),4/(3^0.5)),
    digits = 2
){
  distance_thre = round(distance_thre,digits)

  if (!all(smoothing_type %in% c("view0","view1","view2"))) {
    stop("smoothing_type must be a subset of c('view0','view1','view2').")
  }


  spatialdf = as.matrix(spatialdf)
  spot_by_celltype = as.matrix(spot_by_celltype)
  if (!identical(rownames(spatialdf), rownames(spot_by_celltype))) {
    stop("Row names of spatialdf and spot_by_celltype must be identical.")
  }
  smoothing_type = sort(smoothing_type)

  if (identical(smoothing_type,c('view0', 'view1', 'view2'))) {
    #Make a distance matrix
    distance_matrix = make_distance_matrix(
      spatialdf = spatialdf,
      dm_type = setdiff(smoothing_type,"view0"),
      distance_thre = distance_thre[-1],
      digits = digits
      )

    dst1 = distance_matrix[[smoothing_type[2]]]
    dst2 = distance_matrix[[smoothing_type[3]]]

    spot_by_celltype.dst1=as.matrix(dst1 %*% spot_by_celltype)
    colnames(spot_by_celltype.dst1)=paste0(colnames(spot_by_celltype.dst1),paste0("_",smoothing_type[2]))

    spot_by_celltype.dst2=as.matrix(dst2 %*% spot_by_celltype)
    colnames(spot_by_celltype.dst2)=paste0(colnames(spot_by_celltype.dst2),paste0("_",smoothing_type[3]))

    if (
      identical(rownames(spot_by_celltype),rownames(spot_by_celltype.dst1)) &
      identical(rownames(spot_by_celltype),rownames(spot_by_celltype.dst2))
    ) {
      spot.ct.df=spot_by_celltype %>% cbind(spot_by_celltype.dst1) %>% cbind(spot_by_celltype.dst2)
    }
  }

  if (identical(smoothing_type,c('view0', 'view1'))) {
    #Make a distance matrix
    distance_matrix = make_distance_matrix(
      spatialdf = spatialdf,
      dm_type = setdiff(smoothing_type,"view0"),
      distance_thre = distance_thre[-1],
      digits = digits
      )

    dst1 = distance_matrix[[smoothing_type[2]]]
    spot_by_celltype.dst1=as.matrix(dst1 %*% spot_by_celltype)
    colnames(spot_by_celltype.dst1)=paste0(colnames(spot_by_celltype.dst1),paste0("_",smoothing_type[2]))

    if (
      identical(rownames(spot_by_celltype),rownames(spot_by_celltype.dst1))
    ) {
      spot.ct.df=spot_by_celltype %>% cbind(spot_by_celltype.dst1)
    }
  }

  if (identical(smoothing_type,c("view0"))) {
    spot.ct.df = spot_by_celltype
  }
  if (!exists("spot.ct.df")) {
    stop("Invalid smoothing_type combination.
       Supported: c('view0'), c('view0','view1'), c('view0','view1','view2').")
  }
  spot.ct.df
}
