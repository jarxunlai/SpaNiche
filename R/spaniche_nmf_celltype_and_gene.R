#' spaniche_nmf_celltype_and_gene: Joint factorization of cell-type abundance  and gene expression matrices
#'
#' This function performs spatially regularized integrative non-negative
#' matrix factorization by jointly decomposing a spot-by-cell-type
#' abundance matrix and a spot-by-gene expression matrix.
#'
#' Spatial neighborhood information is incorporated by constructing
#' multi-view representations (view0, view1, view2) based on spot coordinates,
#' followed by graph-regularized joint matrix factorization.
#'
#' The following parameters define the basic inputs of the analysis.
#' @param spatial.seu A Seurat object containing spatial transcriptomics data.
#' The order of spots must be consistent with spatialdf.
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
#' Parameters related to gene expression.
#' In practice, genes may need to be pre-selected (e.g., highly variable genes).
#' @param ngenes The number of highly variable genes.
#' Parameters related to joint matrix factorization.
#' @param topic_num An integer specifying the number of components or topics to be extracted. It determines the number of columns in the W matrix and the number of rows in each H matrix.
#' @param defined_weight Either \code{"default"} to use data-driven modality
#'   weights, or a numeric vector of length two specifying user-defined weights
#'   for the cell-type and gene expression matrices.
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
#'     \item \code{nmf_res}: Results from spatially regularized integrative NMF,
#'     including factor matrices, convergence diagnostics, and cluster assignments.
#'   }
#' @import Seurat
#' @import IntNMF
#' @import scales
#' @export
spaniche_nmf_celltype_and_gene = function(
  spatial.seu,
  spatialdf,
  spot_by_celltype,
  smoothing_type = c('view0', 'view1', 'view2'),
  distance_thre = c(NA,2/(3^0.5),4/(3^0.5)),
  digits = 2,
  ngenes = 2000,
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

  # Construct spatially extended matrices
  spot.ct.df = get_spot_by_celltype_extended(
    spatialdf = spatialdf,
    spot_by_celltype = spot_by_celltype,
    smoothing_type = smoothing_type,
    distance_thre = distance_thre,
    digits = digits
  )
  write.csv(spot.ct.df,file = "spot_by_celltype_extended.csv",quote = F,row.names = T)

  ### Part II: Extract gene expression matrix
  spatial.seu <- NormalizeData(spatial.seu, assay = assay, normalization.method = "LogNormalize", scale.factor = 10000)
  spatial.seu <- FindVariableFeatures(spatial.seu, assay = assay, selection.method = "vst", nfeatures = ngenes)
  genedf <- .get_layer_data(spatial.seu, "data", assay)[VariableFeatures(spatial.seu, assay = assay),] %>% as.matrix() %>% t()

  # Ensure consistent spot ordering
  if(identical(rownames(spatialdf),rownames(genedf))){
    print("The spot order is consistent between spatialdf and genedf.")
  } else {
    genedf=genedf[rownames(spatialdf),]
    print("The spot order between spatialdf and genedf is inconsistent. Adjust the spot order of genedf.")
  }

  # extended matrices
  spot.gene.df = get_spot_by_celltype_extended(
    spatialdf = spatialdf,
    spot_by_celltype = genedf,
    smoothing_type = smoothing_type,
    distance_thre = distance_thre,
    digits = digits
  )

  ### Part III: Joint factorization of cell-type abundance and gene expression matrices

  # df to matrix
  spot.ct.df = as.matrix(spot.ct.df)
  spot.gene.df = as.matrix(spot.gene.df)
  # Estimate initial modality weights prior to running IntNMF
  dat <- list(spot.ct.df,spot.gene.df)
  theta_v = define_initial_weight(dat)
  # The function nmf.mnnals requires the samples to be on rows and variables on columns.
  fit <- nmf.mnnals(dat=dat,k=topic_num,maxiter=50,st.count=10,n.ini=3,ini.nndsvd=TRUE,seed=TRUE,wt = theta_v)
  # Determine more refined modality weights
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
    nmf_res=myfit
  )

  print("Output all results.")
  return(output)
}
