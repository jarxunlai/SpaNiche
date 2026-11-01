#' Visualize NMF component contributions by cell type and view
#'
#' This function visualizes the contribution of a single NMF topic across
#' different cell types and views. Two complementary bar plots are generated: one showing
#' topic contributions across cell types within each view, and another
#' showing contributions across views within each cell type. The two plots
#' are combined vertically using \pkg{patchwork}.
#'
#' @param components A numeric data frame or matrix containing NMF component
#'   loadings. Rows correspond to topics and must have row names (e.g.,
#'   \code{topic1}, \code{topic2}). Columns correspond to cell type and view
#'   combinations, formatted as \code{celltype_viewX}.
#' @param one_topic A string that specifies which topic to visualize. Default is "topic1".
#' @param p1.axis.text.x.bottom.size p1.axis.text.x.bottom.size
#' @param p2.axis.text.x.bottom.size p2.axis.text.x.bottom.size
#' @param p1.strip.text.x.size p1.strip.text.x.size
#' @param p2.strip.text.x.size p2.strip.text.x.size
#' @param color_celltype Either NULL or a named vector of colors. If provided, it specifies the colors to use for the different cell types in the plot. If NULL, default ggplot2 colors are used.
#' @param color_view Either NULL or a named vector of colors. If provided, it specifies the colors to use for the different views in the plot. If NULL, default ggplot2 colors are used.
#'
#' @return A combined ggplot2 plot that visualizes the values of the specified topic grouped by celltype and view.
#' @import patchwork
#' @import reshape2
#' @import scales
#' @import RColorBrewer
#' @export
plot_NMF_components = function(
    components,
    one_topic = "topic1",
    p1.axis.text.x.bottom.size = 8,
    p2.axis.text.x.bottom.size = 9,
    p1.strip.text.x.size = 15,
    p2.strip.text.x.size = 9,
    color_celltype = NULL,
    color_view = NULL
){
  if (is.null(rownames(components))) {
    stop("components must have row names corresponding to topic names.")
  }
  components$topic=rownames(components)
  if (!one_topic %in% components$topic) {
    stop(one_topic, " not found in components.")
  }
  components=components%>%reshape2::melt(id.vars=c("topic"))
  colnames(components)[2]="celltype_view"
  components$celltype_view = as.character(components$celltype_view)

  ### modify
  components$view = str_extract(components$celltype_view,"view[0-2].*$")
  # if (any(is.na(components$view))) {
  #   stop("Failed to parse view from column names. Expect format: celltype_viewX")
  # }
  components$view[is.na(components$view)] = "view0"
  components$celltype=str_replace(components$celltype_view,"_view[0-2].*$","")

  each=components%>%dplyr::filter(topic == one_topic)
  ###
  p1=each%>%ggplot(aes(x=celltype,y=value,fill=celltype))+
    geom_bar(stat = "identity")+
    scale_y_continuous(expand = c(0.02,0))+
    facet_grid(.~view)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,color = "black",size = p1.axis.text.x.bottom.size),
      axis.text.y.left = element_text(color = "black",size = 10),
      axis.title = element_text(size = 15),
      strip.text.x = element_text(size = p1.strip.text.x.size),

      legend.position = "top"
    )

  ###
  p2=each%>%ggplot(aes(x=view,y=value,fill=view))+
    geom_bar(stat = "identity")+
    scale_y_continuous(expand = c(0.02,0))+
    facet_grid(.~celltype)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,color = "black",size = p2.axis.text.x.bottom.size),
      axis.text.y.left = element_text(color = "black",size = 10),
      axis.title = element_text(size = 15),
      strip.text.x = element_text(size = p2.strip.text.x.size),

      legend.position = "top"
    )

  ###
  if (is.null(color_celltype)) {
  } else {
    p1 = p1+scale_fill_manual(values = color_celltype)
  }
  if (is.null(color_view)) {
  } else {
    p2 = p2+scale_fill_manual(values = color_view)
  }

  ###
  #p1+p2+plot_layout(widths = c(1, relative_width))
  p1 / p2
}
