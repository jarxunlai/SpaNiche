#' spaniche_nmf_comm: Spatial integrative NMF of cell-type abundance
#' and ligand–receptor interaction matrices, with downstream pathway inference
#'
#' @param spatial.seu ST data is processed as a Seurat object with 'Spatial' as the assay name. Seurat object is normalized by using the 'NormalizeData' function. The order of spots must be consistent with spatialdf.
#' @param spatialdf a data frame or matrix that contains two columns, "col" and "row". Row names must be
#'   spot barcodes and consistent with \code{spot_by_celltype}.
#' @param spot_by_celltype a data frame or matrix representing the abundance of cell types, with spots as rows and cell types as columns. Its row names should be consistent with the row names of 'spatialdf'.
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
#'   As it involves the merging of matrices, the returned matrix's dimension varies depending on the chosen option. The dimensions corresponding to the three options are 1x, 2x, and 3x of the original dimensions, respectively.
#' @param distance_thre This parameter represents distance, defining the central distance from center to loop1 (view1), and from center to loop2 (view2). For view0, the distance is ignored and should be set to NA. The length of the vector is consistent with the smoothing_type parameter. For example, A typical choice for Visium data is: distance_thre = c(NA,round(2/(3^0.5),2),round(4/(3^0.5),2)).
#' @param digits integer indicating the number of decimal places.
#' @param LRDB This parameter represents the reference for cell-cell interactions, which comes from CellChat.
#' @param LR_dm_type This parameter can be one or more values from c('view1', 'view2'), representing the type of adjacency matrix to be returned.
#' @param LR_distance_thre A numeric vector specifying distance thresholds. If only \code{"view1"} is requested, a single value is sufficient, defining the distance from center to first-order neighbors (view1). If \code{"view1"} + \code{"view2"} is requested, two values are required, defining the distance from center to first-order neighbors (view1) and from center to second-order neighbors (view2).
#' @param var_thre Variance threshold
#' @param topic_num An integer specifying the number of components or topics to be extracted. It determines the number of columns in the W matrix and the number of rows in each H matrix.
#' @param defined_weight Either \code{"default"} to use data-driven modality
#'   weights, or a numeric vector of length two specifying user-defined weights
#'   for the cell-type and LR matrices.
#' @param lambda_v Vector of spatial regularization parameters controlling the
#'   strength of Laplacian smoothness.
#' @param sigma_v Gaussian kernel parameter. Used in the computation of the Laplacian matrix, influencing how spatial information is incorporated into the factorization.
#' @param maxiter The maximum number of iterations allowed for the factorization process. This acts as a stopping criterion to prevent the algorithm from running indefinitely.
#' @param st.count Convergence counter. If the change in the reconstruction error remains below epsilon for st.count consecutive iterations, the algorithm stops, assuming it has converged.
#' @param epsilon Threshold on the relative change in the objective function used to assess convergence. If the relative change in error is less than epsilon for st.count consecutive iterations, convergence is assumed.
#' @param topic_lr_quantile_thre A threshold for quantile selection. Only data above this quantile will be considered for visualization. Default is 0.95.
#' @param spot_topic.quantile.thre Quantile threshold defining topic-high and
#'   topic-low spot groups (default: top and bottom 20 %).
#' @param S_intra_thre S_intra_thre
#' @param term_topn Number of top enriched functional terms retained for visualization.
#' @param num_core Number of CPU cores used for parallel computation in downstream scSeqComm analysis.
#'   Each worker is limited to one BLAS thread so the total stays at \code{num_core}.
#' @param assay Assay holding the expression data. Defaults to \code{"Spatial"}
#'   if present, otherwise the default assay. Seurat v4 and v5 objects are both supported.
#' @param n_threads Number of threads the BLAS/OpenMP libraries may use while
#'   this function runs (restored afterwards). Defaults to
#'   \code{getOption("SpaNiche.n_threads", 1L)}; \code{NULL} leaves the current
#'   setting untouched (OpenBLAS/MKL then use every core).
#'
#' @details Requires the \pkg{scSeqComm} package (\url{https://gitlab.com/sysbiobig/scseqcomm}).
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{LRintegratedmatrix}: Filtered spatial LR interaction matrix.
#'     \item \code{nmf_res}: Results from spatially regularized integrative NMF.
#'     \item \code{go_res}: GO enrichment results for topic-specific signaling. Only terms with a p-value < 0.01 are retained.
#'     \item \code{go_plot}: Visualization of enriched functional pathways. No p-value filtering is applied; the top "term_topn" terms are plotted.
#'   }
#' @import Seurat
#' @import IntNMF
#' @import scales
#' @export
spaniche_nmf_comm = function(
    ## Basic inputs
  spatial.seu,
  spatialdf,
  spot_by_celltype,
  smoothing_type = c('view0', 'view1', 'view2'),
  distance_thre = c(NA,2/(3^0.5),4/(3^0.5)),
  digits = 2,
  ## Parameters related to ligand–receptor interactions
  LRDB = SpaNiche::CellChatDB.human,
  LR_dm_type = c("view1", "view2"),
  LR_distance_thre = c(2/(3^0.5),4/(3^0.5)),
  var_thre = c(0.2,0.99),
  ## Parameters related to joint matrix factorization
  topic_num = 15,
  defined_weight = "default",
  lambda_v =c(0.5, 1, 2),
  sigma_v = c(0.5, 1, 1.5),
  maxiter = 200,
  st.count = 20,
  epsilon = 1e-04,
  # Downstream cellular responses
  topic_lr_quantile_thre = 0.95,
  spot_topic.quantile.thre = 0.2,
  S_intra_thre = 0.8,
  term_topn = 20,
  num_core = 1,
  assay = NULL,
  n_threads = getOption("SpaNiche.n_threads", 1L)
){
  if (!requireNamespace("scSeqComm", quietly = TRUE)) {
    stop("spaniche_nmf_comm() requires the scSeqComm package: ",
         "remotes::install_gitlab('sysbiobig/scseqcomm')", call. = FALSE)
  }
  old_threads = .set_threads(n_threads)
  on.exit(.restore_threads(old_threads), add = TRUE)
  assay = .resolve_assay(spatial.seu, assay)
  distance_thre = round(distance_thre,digits)
  LR_distance_thre = round(LR_distance_thre,digits)


  # Ensure consistent spot ordering
  if(identical(rownames(spatialdf),colnames(spatial.seu))){
    print("The spot order is consistent between spatialdf and spatial.seu.")
  } else {
    spatialdf = spatialdf[colnames(spatial.seu),]
    print("The spot order between spatialdf and spatial.seu is inconsistent. Adjust the spot order of spatialdf.")
  }

  if(identical(rownames(spatialdf),rownames(spot_by_celltype))){
    print("The spot order is consistent between spatialdf and spot_by_celltype.")
  } else {
    spot_by_celltype=spot_by_celltype[rownames(spatialdf),]
    print("The spot order between spatialdf and spot_by_celltype is inconsistent. Adjust the spot order of spot_by_celltype.")
  }

  # Construct spatially extended cell-type enrichment matrix
  spot.ct.df = get_spot_by_celltype_extended(
    spatialdf = spatialdf,
    spot_by_celltype = spot_by_celltype,
    smoothing_type = smoothing_type,
    distance_thre = distance_thre,
    digits = digits
  )
  write.csv(spot.ct.df,file = "spot_by_celltype_extended.csv",quote = F,row.names = T)

  ### Part II: Construction of ligand–receptor integrated matrix
  step1_res=prepare_for_LRintegratedmatrix(spatial.seu,LRDB,assay = assay)
  distance_matrix = make_distance_matrix_LR(spatialdf = spatialdf,dm_type = LR_dm_type,distance_thre = LR_distance_thre,digits = digits)
  LRelementmatrix = step1_res[[1]]
  cc_interaction = step1_res[[2]]

  # Receptor expression within each spot
  rmat = LRelementmatrix %>% t() %>% .[,cc_interaction$receptor]
  # Ligand expression within each spot
  lmat = LRelementmatrix %>% t() %>% .[,cc_interaction$ligand]
  # NOTE: This section may be generalized in the future; currently assumes four LR configurations (m1–m4).
  m1=generate_LRintegratedmatrix(whoisneighbor = distance_matrix$view1,neighbor_geneexp = rmat,center_geneexp = lmat,center_is_ligand = T)
  m2=generate_LRintegratedmatrix(whoisneighbor = distance_matrix$view1,neighbor_geneexp = lmat,center_geneexp = rmat,center_is_ligand = F)
  colnames(m1)=paste0(colnames(m1)," (view1)")
  colnames(m2)=paste0(colnames(m2)," (view1)")

  m3=generate_LRintegratedmatrix(whoisneighbor = distance_matrix$view2,neighbor_geneexp = rmat,center_geneexp = lmat,center_is_ligand = T)
  m4=generate_LRintegratedmatrix(whoisneighbor = distance_matrix$view2,neighbor_geneexp = lmat,center_geneexp = rmat,center_is_ligand = F)
  colnames(m3)=paste0(colnames(m3)," (view2)")
  colnames(m4)=paste0(colnames(m4)," (view2)")

  LRintegratedmatrix=m1 %>% cbind(m2) %>% cbind(m3) %>% cbind(m4)
  LRintegratedmatrix[is.na(LRintegratedmatrix)] = 0
  saveRDS(LRintegratedmatrix,file = "LRintegratedmatrix.rds")
  rm(list = paste0("m",1:4))

  ### Compute Moran's I and annotate ligand–receptor interactions ###
  moransi_output_df = moransi_for_LRintegratedmatrix(t(LRintegratedmatrix),spatial.seu,cc_interaction)
  saveRDS(moransi_output_df,file = "LRintegratedmatrix_with_moransi.rds")
  print("Moran's index calculation completed")

  ### Part III: Joint factorization of cell-type and LR interaction matrices
  LRintegratedmatrix = LRintegratedmatrix[,colSums(LRintegratedmatrix) > 0]
  colnames(LRintegratedmatrix)=colnames(LRintegratedmatrix) %>% str_replace_all("_","-")
  LRintegratedmatrix=LRintegratedmatrix[spot.ct.df %>% rownames(),]

  # Variance-based filtering of LR interactions
  var.df1 = LRintegratedmatrix %>% apply(2, function(x){var(x)}) %>% as.data.frame()
  colnames(var.df1) = "var_value"
  var.df1$interaction = rownames(var.df1)
  var.df1 = var.df1 %>% arrange(var_value)
  var.df1$index = 1:length(var.df1$var_value)
  var_thre_1 = quantile(var.df1$var_value,var_thre[1])
  var_thre_2 = quantile(var.df1$var_value,var_thre[2])
  var.df1 = var.df1 %>% dplyr::filter(var_value > var_thre_1 & var_value < var_thre_2)
  used_interaction = var.df1$interaction
  used_interaction = intersect(colnames(LRintegratedmatrix),used_interaction)
  LRintegratedmatrix=LRintegratedmatrix[,used_interaction]
  print("Filtering the LR integrated matrix based on variance and selecting some LR pairs.")

  # to matrix
  spot.ct.df = as.matrix(spot.ct.df)
  LRintegratedmatrix = as.matrix(LRintegratedmatrix)
  # Estimate initial modality weights prior to running IntNMF
  dat <- list(spot.ct.df,LRintegratedmatrix)
  theta_v = define_initial_weight(dat)
  # The function nmf.mnnals requires the samples to be on rows and variables on columns.
  fit <- nmf.mnnals(dat=dat,k=topic_num,maxiter=50,st.count=10,n.ini=3,ini.nndsvd=TRUE,seed=TRUE,wt = theta_v)
  # Determine refined modality weights
  if (length(defined_weight) == 2 & all(defined_weight > 0)) {
    rho_v = define_good_weight(X = dat,Wzero = fit$W,Hzero = fit$H,theta_v = theta_v,weight_input = defined_weight)
  } else if (length(defined_weight) == 1 & defined_weight == "default") {
    rho_v = theta_v
  } else {
    stop("defined_weight error")
  }
  print("Weights for different matrices have been determined, now starting joint factorization of the two matrices:")

  # Joint matrix factorization
  myfit = nmf_modified(
    dat = dat,
    Wzero = fit$W,
    Hzero = fit$H,
    wt = rho_v,
    k=topic_num,
    lambda = lambda_v,
    sigma = sigma_v,
    maxiter = maxiter,
    st.count = st.count,
    spatialdf = as.matrix(spatialdf[,c("row","col")]),
    epsilon = epsilon
  )
  saveRDS(myfit,file = "nmf_modified_fit.rds")
  print("joint factorization has been completed.")

  ### Part IV: Linking latent topics to downstream cellular responses
  spot_topic = as.data.frame(myfit$W)
  colnames(spot_topic) = paste0("topic",1:topic_num)
  allSB = rownames(spot_topic)

  gene_expr_matrix = .get_layer_data(spatial.seu, "data", assay)
  gene_expr_matrix = as.matrix(gene_expr_matrix)
  gene_expr_matrix = gene_expr_matrix[rowSums(gene_expr_matrix) > 0,]

  tmpplot = plot_topic_lr(
    topic_lr = t(myfit$H$H2),
    topic_lr_quantile_thre = topic_lr_quantile_thre
  )
  topic_lr_small_file = dir(getwd(),"topic_lr_small_.*csv")
  topic_lr_small = read.csv(topic_lr_small_file)

  # Color palette for topic visualization
  color_cluster = hue_pal()(topic_num)
  names(color_cluster) = paste0("topic",1:topic_num)

  # Objects to be exported
  plot.list=list()
  all.go = data.frame()

  print("Exploring pathway changes due to cell interactions:")

  for (ti in colnames(spot_topic)) {
    tmpst = spot_topic[,ti]
    tmpquan = quantile(tmpst,c(spot_topic.quantile.thre,1 - spot_topic.quantile.thre))
    low.index = which(tmpst <= tmpquan[1])
    high.index = which(tmpst > tmpquan[2])
    lowSB = allSB[low.index]
    highSB = allSB[high.index]


    ########## Load scRNA-seq data
    tmp_gene_expr_matrix = gene_expr_matrix[,c(highSB,lowSB)]
    tmp_anno = data.frame(
      "Cell_ID" = c(highSB,lowSB),
      "Cluster_ID" = c(
        rep(paste(ti,"high",sep = "_"),length(highSB)),
        rep(paste(ti,"low",sep = "_"),length(lowSB))
      )
    )
    ########## Ligand-receptor pairs
    tmp_topic_lr = topic_lr_small %>% dplyr::filter(topic == ti)
    tmp_topic_lr$lr = str_remove(tmp_topic_lr$LRinteraction," .*$")
    tmp_topic_lr$r = ifelse(str_detect(tmp_topic_lr$lr,">"),
                            str_remove(tmp_topic_lr$lr,"^.*>"),
                            str_remove(tmp_topic_lr$lr,"<.*$"))
    tmp_topic_lr$l = ifelse(str_detect(tmp_topic_lr$lr,">"),
                            str_remove(tmp_topic_lr$lr,"->.*$"),
                            str_remove(tmp_topic_lr$lr,"^.*<-"))
    tmp_topic_lr=tmp_topic_lr[,c("l","r")] %>% unique()
    colnames(tmp_topic_lr) = c("ligand","receptor")

    tmp_topic_lr$ligand = str_replace(tmp_topic_lr$ligand,"^HLA-","xxx") %>% str_replace("-",",") %>% str_replace("xxx","HLA-")
    tmp_topic_lr$receptor = str_replace(tmp_topic_lr$receptor,"-",",")
    LR_db = tmp_topic_lr
    ########## Transcriptional regulatory networks
    TF_TG_db <- scSeqComm::TF_TG_TRRUSTv2_HTRIdb_RegNetwork_High
    ########## Receptor-Transcription factor a-priori association
    TF_PPR <- scSeqComm::TF_PPR_KEGG_human


    ########## Identify and quantify intercellular and intracellular signaling
    scSeqComm_res <- .with_single_thread_workers(scSeqComm::scSeqComm_analyze(gene_expr = tmp_gene_expr_matrix,
                                       cell_metadata = tmp_anno,
                                       inter_signaling = F,
                                       LR_pairs_DB = LR_db,
                                       TF_reg_DB = TF_TG_db,
                                       R_TF_association = TF_PPR,
                                       N_cores = num_core,
                                       DEmethod = "wilcoxon"),
                                       n_workers = num_core)


    ### Downstream cellular responses induced by ligand–receptor interactions
    topic_high_comm = dplyr::filter(
      scSeqComm_res$comm_results,
      cluster == paste(ti,"high",sep = "_") & S_intra >= S_intra_thre
      )
    topic_high_comm = as.data.frame(topic_high_comm) #20260426
    if (dim(topic_high_comm)[1] < 2) {next} #20231018

    # GO analysis of topic_high communication
    geneUniverse <- unique(unlist(scSeqComm_res$TF_reg_DB_scrnaseq))
    cell_functional_response <- scSeqComm::scSeqComm_GO_analysis(
      results_signaling = topic_high_comm,
      geneUniverse = geneUniverse,
      method = "general")

    # plot
    cell_functional_response$pval = as.numeric(cell_functional_response$pval)
    cell_functional_response = cell_functional_response %>% arrange(pval)
    cell_functional_response$pval_log10_neg = -log10(cell_functional_response$pval)
    cell_functional_response$cluster = ti
    all.go = rbind(all.go,cell_functional_response %>% dplyr::filter(pval < 0.01))
    tmpres = cell_functional_response %>% slice_head(n = term_topn)

    tmpres=tmpres%>%arrange(pval_log10_neg)
    tmpres$Term=factor(tmpres$Term,levels = tmpres$Term)

    tmpbar = tmpres %>% ggplot(aes(x=Term,y=pval_log10_neg))+
      geom_hline(yintercept = -log10(0.01),color = "black",alpha=0.7)+
      geom_bar(stat="identity",alpha=0.8,aes(fill=cluster))+
      geom_text(mapping = aes(x=Term,y=0,label=Term),hjust=0)+
      scale_x_discrete("")+
      scale_y_continuous("-log10(p value)",expand = c(0.02,0))+
      scale_fill_manual(values = color_cluster)+
      coord_flip()+
      labs(title = ti)+
      theme_bw()+
      theme(
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x.bottom = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 20)
      )

    index=which(ti == colnames(spot_topic))
    plot.list[[get("index")]]=tmpbar
  }
  print("Completed.")

  output = list(
    LRintegratedmatrix=LRintegratedmatrix,
    nmf_res=myfit,
    go_res=all.go, # p-values filtered with threshold pval < 0.01
    go_plot=plot.list # no p-value filtering, only the top "term_topn" terms are plotted
  )

  print("Output all results.")
  return(output)
}
