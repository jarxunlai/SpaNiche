#' Visualize the 2D spatial distribution of a single feature
#'
#' This function visualizes the spatial distribution of a single feature
#' (e.g., gene expression, topic score, or interaction strength) on a 2D
#' coordinate system, typically derived from spatial transcriptomics or
#' imaging-based data. Both point-based and density-based (Nebulosa-style)
#' visualizations are supported.
#'
#' @param feature_matrix A numeric matrix with features as rows and
#'   samples/spots/cells as columns. Row names must contain \code{one_feature},
#'   and column names must match the \code{SB} column in \code{coordinate}.
#' @param one_feature A character string specifying the name of a single feature
#'   (row name of \code{feature_matrix}) to visualize.
#' @param coordinate A data frame or matrix containing spatial coordinates.
#'   It must include three columns:
#'   \itemize{
#'     \item \code{SB}: sample/spot/cell identifier
#'     \item \code{row}: y-axis spatial coordinate
#'     \item \code{col}: x-axis spatial coordinate
#'   }
#' @param thre_to_zero A numeric value between 0 and 1 indicating the quantile
#'   threshold below which expression values will be set to zero. This can be
#'   used to suppress low-level background noise. Values outside [0, 1) will
#'   result in an error.
#' @param color_palette A character vector specifying the color scale.
#'   If a single value is provided, it must be one of \code{"viridis"},
#'   \code{"magma"}, \code{"inferno"}, \code{"plasma"}, or \code{"cividis"}.
#'   If two values are provided, they will be used as the low and high colors
#'   for a continuous gradient.
#' @param visualization_type Character string specifying the visualization type.
#'   Either \code{"point"} (default) for raw point visualization or
#'   \code{"density"} for weighted kernel density estimation using the
#'   \pkg{Nebulosa} framework.
#' @param pt.size Numeric value controlling the size of points.
#' @param axis.title.size Numeric value controlling the font size of axis titles.
#' @param axis.text.size Numeric value controlling the font size of axis tick labels.
#' @param plot.title.size Numeric value controlling the font size of the plot title.
#' @param raster Logical value indicating whether to rasterize the plot using
#'   \pkg{ggrastr}. This is recommended for datasets with a large number of points. By default, it's set to FALSE.
#'
#' @return A \code{ggplot2} object representing the spatial distribution of the
#'   selected feature.
#' @export
distribution_2d=function(
    feature_matrix,
    one_feature,
    coordinate,
    thre_to_zero,
    color_palette = c("lightgrey","blue"),
    visualization_type = "point",
    pt.size = 1,
    axis.title.size = 18,
    axis.text.size = 12,
    plot.title.size = 14,
    raster = F
){

  if (!one_feature %in% rownames(feature_matrix)) {
    stop("Feature '", one_feature, "' not found in feature_matrix.")
  }
  one.feature.exp=feature_matrix[one_feature,] %>% as.data.frame()
  colnames(one.feature.exp) = "exp"
  one.feature.exp$SB=rownames(one.feature.exp)


  n_before <- nrow(one.feature.exp)
  one.feature.exp <- inner_join(one.feature.exp, coordinate, by = "SB")
  n_after <- nrow(one.feature.exp)
  if (n_after < n_before) {
    warning("Some samples were dropped during coordinate join.")
  }
  ### 归0
  if (thre_to_zero >= 0 & thre_to_zero < 1) {
    thre1 = quantile(one.feature.exp$exp,thre_to_zero)
    one.feature.exp$exp[one.feature.exp$exp < thre1] = 0
  } else {
    stop("range for 'thre_to_zero' can only be [0,1).")
  }


  if(visualization_type == "density"){
    if (!requireNamespace("Nebulosa", quietly = TRUE)) {
      stop("visualization_type = 'density' requires the Nebulosa package (Bioconductor).")
    }
    z=utils::getFromNamespace("calculate_density", "Nebulosa")(
      w = one.feature.exp[,1],
      x = one.feature.exp[, c("col", "row")],
      method = "wkde",
      adjust = 1, map = TRUE
    )

    colnames(one.feature.exp)[1] = "original_exp"
    one.feature.exp$exp = z
  }


  p=ggplot(data = one.feature.exp,aes(x=col,y=row,color=exp))+
    geom_point(size=pt.size,shape = 16)+
    labs(title = one_feature)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),

      axis.title = element_text(size = axis.title.size),
      axis.text = element_text(color="black",size = axis.text.size),
      axis.ticks.length = unit(0.15,"cm"),

      plot.title = element_text(hjust = 0.5,size = plot.title.size)
    )


  if (length(color_palette) == 1) {
    p=p+scale_color_viridis_c(option = color_palette)
  } else if (length(color_palette) == 2) {
    p=p+scale_color_gradient(low = color_palette[1],high = color_palette[2])
  }

  # if(raster){
  #   ggrastr::rasterise(p, dpi = 300)
  # }else{
  #   p
  # }
  p <- if (raster) ggrastr::rasterise(p, dpi = 300) else p
  return(p)
}
