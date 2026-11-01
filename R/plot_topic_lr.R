#' Plot topic-specific ligand–receptor (LR) interaction heatmap
#'
#' This function visualizes topic-specific ligand–receptor (LR) interactions
#' as a heatmap. For each topic, LR values are first normalized to
#' relative contributions, then filtered by a global quantile threshold.
#' The top LR interactions for each topic are selected and displayed
#' in separate facets.
#'
#' @param topic_lr A numeric data frame or matrix with LR interactions as rows
#'  and topics as columns. Row names must correspond to LR interaction names.
#' @param topic_lr_quantile_thre A numeric value between 0 and 1 specifying the
#'   global quantile threshold used to filter LR values across all topics.
#'   Only LR–topic pairs with values above this quantile will be considered.
#'   Default is \code{0.95}.
#' @param top_num_plot The maximum number of LR interactions to be displayed for each topic in the plot. Default is 10.
#' @param fill_viridis Specifies the "viridis" color palette option to be used for the heatmap. Default is "F".
#' @param axis.title.size The font size of the axis titles. Default is 14.
#' @param axis.text.x.bottom.size The font size for the text on the bottom of the x-axis. Default is 12.
#' @param axis.text.y.left.size The font size for the text on the left of the y-axis. Default is 8.
#' @param strip.text.y.size The font size for the text on the y-axis strip. Default is 12.
#'
#' @return A \code{ggplot2} object representing a faceted heatmap of topic-specific LR interactions.
#' @export
plot_topic_lr = function(
    topic_lr,
    topic_lr_quantile_thre = 0.95,
    top_num_plot = 10,
    fill_viridis = "F",
    axis.title.size = 14,
    axis.text.x.bottom.size = 12,
    axis.text.y.left.size = 8,
    strip.text.y.size = 12
){
  topic_num = ncol(topic_lr)
  topic_lr = topic_lr %>% apply(2,function(x){x*1000/sum(x)}) %>% as.data.frame()
  colnames(topic_lr) = paste0("topic",1:topic_num)
  min_thre = quantile(as.numeric(as.matrix(topic_lr)),topic_lr_quantile_thre)

  topic_lr_copy = topic_lr
  topic_lr$LRinteraction = rownames(topic_lr)
  topic_lr = topic_lr%>%reshape2::melt(id.vars=c("LRinteraction"))
  colnames(topic_lr)[2:3] = c("topic","exp")
  topic_lr$topic = factor(topic_lr$topic,levels = paste0("topic",1:topic_num))

  topic_lr_small = topic_lr[topic_lr$exp > min_thre,]
  topic_lr_small = topic_lr_small %>% arrange(topic,desc(exp))
  write.csv(topic_lr_small,file = paste0("topic_lr_small_",topic_lr_quantile_thre,".csv"),quote = T,row.names = F)
  ### Relatively complete results are saved here;
  ### for visualization, only a subset (top LR per topic) will be plotted.

  topic_lr_plot = topic_lr_small %>% group_by(topic) %>% top_n(top_num_plot,wt = exp)
  topic_lr_plot$LRinteraction2 = make.unique(topic_lr_plot$LRinteraction,sep = "_")
  repeat_lr = topic_lr_plot$LRinteraction2[str_detect(topic_lr_plot$LRinteraction2,"_[0-9]{1,}$")]
  repeat_lr2 = str_remove(repeat_lr,"_[0-9]{1,}$")

  tmp1 = topic_lr_copy[topic_lr_plot$LRinteraction %>% unique(),]
  tmp2 = topic_lr_copy[repeat_lr2,]
  rownames(tmp2) = repeat_lr
  topic_lr_small = rbind(tmp1,tmp2)

  topic_lr_small = topic_lr_small[topic_lr_plot$LRinteraction2,levels(topic_lr_plot$topic)]
  topic_lr_small$LRinteraction2 = rownames(topic_lr_small)
  topic_lr_small = topic_lr_small%>%reshape2::melt(id.vars=c("LRinteraction2"))
  colnames(topic_lr_small)[2:3] = c("topic","exp")
  topic_lr_small$LRinteraction2 = factor(topic_lr_small$LRinteraction2,levels = topic_lr_plot$LRinteraction2 %>% rev())
  topic_lr_small$topic = factor(topic_lr_small$topic,levels = levels(topic_lr_plot$topic))

  lranno = topic_lr_plot[,c("topic","LRinteraction2")]
  colnames(lranno)[1] = "topic_top_lr"
  lranno$topic_top_lr = paste0(lranno$topic_top_lr,"_top_LR")
  topic_lr_small = topic_lr_small %>% inner_join(lranno,by = "LRinteraction2")
  topic_lr_small$topic_top_lr = factor(topic_lr_small$topic_top_lr,levels = paste0(levels(topic_lr_plot$topic),"_top_LR"))

  ###画图
  topic_lr_small %>% ggplot(aes(x=topic,y=LRinteraction2))+
    geom_tile(aes(fill = exp))+
    scale_x_discrete("",expand = c(0,0))+
    scale_y_discrete("LR interaction",expand = c(0,0))+
    scale_fill_viridis_c("relative\nvalue",option = fill_viridis)+
    facet_grid(topic_top_lr~., scales = "free_y", space = "free_y")+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = axis.title.size),
      axis.text.x.bottom = element_text(size = axis.text.x.bottom.size,color = "black",angle = 45,hjust = 1),
      axis.text.y.left = element_text(size = axis.text.y.left.size,color = "black"),
      axis.ticks.length = unit(0.15,"cm"),

      strip.text.y = element_text(size = strip.text.y.size),
    )
}
