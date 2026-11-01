#' spaniche_nmf_celltype_and_lr: Joint factorization of cell-type abundance  and ligand-receptor interaction matrices
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
#' @param assay Assay holding the expression data. Defaults to \code{"Spatial"}
#'   if present, otherwise the default assay. Seurat v4 and v5 objects are both supported.
#' @param n_threads Number of threads the BLAS/OpenMP libraries may use while
#'   this function runs (restored afterwards). Defaults to
#'   \code{getOption("SpaNiche.n_threads", 1L)}; \code{NULL} leaves the current
#'   setting untouched (OpenBLAS/MKL then use every core).
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{LRintegratedmatrix}: Filtered ligand–receptor interaction matrix.
#'     \item \code{nmf_res}: Results from spatially regularized integrative NMF.
#'   }
#' @import Seurat
#' @import tidyverse
#' @import IntNMF
#' @import scales
#' @export
spaniche_nmf_celltype_and_lr = function(
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
  assay = NULL,
  n_threads = getOption("SpaNiche.n_threads", 1L)
){
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


  output = list(
    LRintegratedmatrix=LRintegratedmatrix,
    nmf_res=myfit
  )

  print("Output all results.")
  return(output)
}
