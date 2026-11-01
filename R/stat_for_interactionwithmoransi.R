#' stat_for_interactionwithmoransi: Summarize and visualize spatially significant ligand–receptor interactions
#'
#' @description
#' This function filters ligand–receptor (LR) interactions based on
#' Moran's I statistics and visualizes the distribution of spatially
#' significant interactions across functional pathways and annotation
#' categories.
#' @param moransi_df A data frame expected to contain columns such as MoransI_observed, MoransI_p.value, pathway_name, and annotation.
#' @param moransi_thre Threshold for the Moran's I value to filter the interactions.
#' @param p.value_thre Threshold for the p-value to filter the interactions.
#' @param height_ratio Ratio to define the relative heights of the two plots when they are combined.
#'
#' @return
#' A \code{patchwork} object consisting of two vertically stacked bar plots:
#' \enumerate{
#'   \item A bar plot showing the number of spatially significant LR interactions
#'   per pathway.
#'   \item A bar plot showing the number of spatially significant LR interactions
#'   per annotation category.
#' }
#' The function also produces a CSV file as a side effect.
#'
#' @import patchwork
#' @export

stat_for_interactionwithmoransi = function(
    moransi_df,
    moransi_thre = 0.2,
    p.value_thre = 0.05,
    height_ratio = c(5,2)
){

  goodpair=moransi_df %>% dplyr::filter(MoransI_observed > moransi_thre & MoransI_p.value < p.value_thre)
  write.csv(goodpair,file = paste0("goodpair_MoransI_observed_",moransi_thre,"_MoransI_p.value_",p.value_thre,".csv"),quote = T,row.names = F)

  ###
  goodpair$pathway_name=factor(
    goodpair$pathway_name,
    levels = table(goodpair$pathway_name) %>% sort() %>% names()
  )
  p1 = goodpair %>% ggplot(aes(y=pathway_name))+
    geom_bar(fill = "lightblue")+
    theme_minimal()+
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(color="black",size = 12),
      axis.title.x.bottom = element_blank()
    )

  ###
  goodpair$annotation=factor(
    goodpair$annotation,
    levels = table(goodpair$annotation) %>% sort() %>% names()
  )
  p2 = goodpair %>% ggplot(aes(y=annotation))+
    geom_bar(aes(fill = annotation))+
    theme_minimal()+
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(color="black",size = 12),
      legend.position = "none"
    )

  ###
  p1 / p2 + plot_layout(heights = height_ratio)
}
