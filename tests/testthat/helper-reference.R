# Reference implementations: the original dense SpaNiche code (v1.1.0), kept
# verbatim except raster::pointDistance(x, lonlat = FALSE) -> as.matrix(dist(x)),
# which computes the same Euclidean distances. Used as test oracles.

ref_make_distance_matrix = function(spatialdf, dm_type = c("view1", "view2"), distance_thre, digits = 2){
  distance_thre = round(distance_thre,digits)
  spatialdf=as.matrix(spatialdf[,c("col","row")])
  dm_type=sort(dm_type)
  distance_thre=sort(distance_thre)
  names(distance_thre) = dm_type
  dst <- as.matrix(dist(spatialdf)); dimnames(dst) <- NULL
  dst=round(dst,digits)
  res_list=as.list(1:length(dm_type))
  names(res_list)=dm_type
  for (view in intersect(c("view1", "view2"), dm_type)) {
    thre1 = distance_thre[view]
    d=dst
    d[d <= thre1] = 1
    d[d > thre1] = 0
    colnames(d)=rownames(spatialdf)
    rownames(d)=rownames(spatialdf)
    coef1=rowSums(d)^(-1)
    res_list[[view]] = apply(d, 2, function(x){ x*coef1})
  }
  res_list
}

ref_make_distance_matrix_LR = function(spatialdf, dm_type = c("view1", "view2"), distance_thre, digits = 2){
  distance_thre = round(distance_thre,digits)
  spatialdf=as.matrix(spatialdf[,c("col","row")])
  dm_type=sort(dm_type)
  distance_thre=sort(distance_thre)
  dst <- as.matrix(dist(spatialdf)); dimnames(dst) <- NULL
  dst=round(dst,digits)
  res_list=as.list(1:length(dm_type))
  names(res_list)=dm_type
  if ("view1" %in% dm_type) {
    thre1 = distance_thre[1]
    d=dst
    d[d <= thre1] = 1
    d[d > thre1] = 0
    colnames(d)=rownames(spatialdf); rownames(d)=rownames(spatialdf)
    coef1=rowSums(d)^(-1)
    res_list$view1 = apply(d, 2, function(x){ x*coef1})
  }
  if ("view2" %in% dm_type) {
    thre1 = distance_thre[1]; thre2 = distance_thre[2]
    d=dst
    d[d == 0] = 1
    d[d > 0 & d <= thre1] = 0
    d[d > thre1 & d <= thre2] = 1
    d[d > thre2] = 0
    colnames(d)=rownames(spatialdf); rownames(d)=rownames(spatialdf)
    coef2=rowSums(d)^(-1)
    res_list$view2 = apply(d, 2, function(x){ x*coef2})
  }
  res_list
}

ref_generate_LRintegratedmatrix = function(whoisneighbor, neighbor_geneexp, center_geneexp){
  (whoisneighbor %*% neighbor_geneexp) * center_geneexp
}

ref_nmf_modified = function(
  dat = dat,
  Wzero = fit$W,
  Hzero = fit$H,
  wt = rho_v,
  k=15,
  lambda = 1,
  sigma = 0.5,
  maxiter = 50,
  st.count = 20,
  spatialdf = as.matrix(spatialdf[,c("row","col")]),
  epsilon = 1e-04
  ){
  # I would like to express my gratitude to TJJ and CZH
  # for their invaluable assistance in completing the core part of this code.
  # I consulted with them on numerous mathematical issues
  # and am deeply appreciative of their insights and explanations.
  # Additionally, I would like to thank the R package IntNMF,
  # as our preliminary concepts were inspired by it.


  #2024.07
  #Should a random seed be set during this process?

  # Basic variables
  n = nrow(Wzero)
  n.dat = length(Hzero)
  if (n.dat != length(wt)) stop("Number of weights must match number of data")

  # Since multiple values of lambda and sigma are explored, results are stored as lists.
  min.f.WH <- NULL
  W.list = NULL
  for (i in 1:n.dat) {assign(paste("H",i,".list",sep = ""),NULL)}
  convergence.list = NULL

  for (lambdai in lambda) {
    for (sigmai in sigma) {
      count = 0
      iter = 1
      convergence = NULL
      W=Wzero
      H=Hzero

      # Compute adjacency, degree, and Laplacian matrices using a Gaussian kernel based on spatial coordinates
      laplacian_res = ref_laplacian(
        spa_mat = spatialdf,sigma = sigmai
      )
      laplacian_matrix = laplacian_res$laplacian
      degree_matrix = laplacian_res$degree
      adj_matrix = laplacian_res$adj

      while ((iter < maxiter) & (count < st.count)) {

        #step1: H->W
        Aminus=matrix(0,nrow = n,ncol = k)
        Aplus=matrix(0,nrow = n,ncol = k)
        for (i in 1:n.dat) {
          Aminus = Aminus + 2*wt[i] * dat[[i]] %*% t(H[[i]])
          Aplus = Aplus + 2 * wt[i] * W %*% H[[i]] %*% t(H[[i]])
        }

        # right_term = 2*lambdai*laplacian_matrix %*% W
        # Bminus = matrix(0,nrow = n,ncol = k)
        # Bplus = matrix(0,nrow = n,ncol = k)
        # Bminus = ifelse(right_term < 0,-right_term,0)
        # Bplus = ifelse(right_term >= 0,right_term,0)
        Bminus = 2*lambdai* adj_matrix %*% W
        Bplus = 2*lambdai* degree_matrix %*% W
        W = W * ((Aminus+Bminus) / pmax(Aplus+Bplus,.Machine$double.eps))

        #step2: W->H
        for (i in 1:n.dat) {
          Aminus=matrix(0,nrow = k,ncol = ncol(Hzero[[i]]))
          Aplus=matrix(0,nrow = k,ncol = ncol(Hzero[[i]]))
          Aminus = 2 * wt[i] * t(W) %*% dat[[i]]
          Aplus = 2 * wt[i] * t(W) %*% W %*% H[[i]]
          Aplus = pmax(Aplus,.Machine$double.eps)
          H[[i]] = H[[i]] * (Aminus / Aplus)
        }

        # Compute the following for each iteration
        # (i) abs.diff : Change in reconstruction error
        # (ii) f.WH    : Objective function
        # (iii) difference_level
        if (iter==1) {
          for (i in 1:n.dat) {
            assign(paste("d",i,".old",sep=""),W %*% H[[i]])
          }

          f.WH <- 0
          for (i in 1:n.dat) {
            f.WH <- f.WH + wt[i]*sum(
              (dat[[i]] - eval(parse(text=paste("d",i,".old",sep=""))))^2
            )
          }
          f.WH <- f.WH + lambdai * sum(diag(t(W) %*% laplacian_matrix %*% W))

          loss.previous = f.WH
          loss.current = NA
          difference_level = NA
          abs.diff = NA
        } else {
          for (i in 1:n.dat) {
            assign(paste("d",i,".new",sep=""),W %*% H[[i]])
          }

          abs.diff <- 0
          for (i in 1:n.dat) {
            abs.diff <- abs.diff + sum(
              abs(
                eval(parse(text=paste("d",i,".new",sep=""))) - eval(parse(text=paste("d",i,".old",sep="")))
              )
            ) / sum(
              eval(parse(text=paste("d",i,".old",sep="")))
            )
          }

          for (i in 1:n.dat) {
            assign(paste("d",i,".old",sep=""),eval(parse(text=paste("d",i,".new",sep=""))))
          }

          f.WH <- 0
          for (i in 1:n.dat) {
            f.WH <- f.WH + wt[i]*sum(
              (dat[[i]] - eval(parse(text=paste("d",i,".new",sep=""))))^2
            )
          }
          f.WH <- f.WH + lambdai * sum(diag(t(W) %*% laplacian_matrix %*% W))

          loss.current = f.WH
          difference_level = abs(loss.current - loss.previous) / loss.previous
          loss.previous = loss.current
        }

        if (!is.na(difference_level) & difference_level < epsilon) {count <- count + 1} else {count <- 0}

        #Store convergence-related metrics
        convergence <- rbind(convergence,c(iter,count,ifelse(!is.na(difference_level) & difference_level < epsilon,1,0),difference_level,abs.diff,f.WH))
        iter <- iter + 1
      } # End of while loop

      colnames(convergence) <- c("iter","count","stability","difference_level","abs.diff","f.WH")
      convergence = as.data.frame(convergence)
      convergence$lambda = lambdai
      convergence$sigma = sigmai

      # Store all candidate W, H, and convergence results
      min.f.WH <- c(min.f.WH,f.WH)
      W.list <- c(W.list,list(W))
      for (i in 1:n.dat) {
        assign(paste("H",i,".list",sep=""),
               c(
                 eval(parse(text=paste("H",i,".list",sep=""))),list(H[[i]])
               )
        )
      }
      convergence.list <- c(convergence.list,list(convergence))

    }
  } # End of lambda and sigma loops
  H.list <- NULL
  for (i in 1:n.dat) H.list <- c(H.list,list(eval(parse(text=paste("H",i,".list",sep="")))))
  names(H.list) <- paste("H",1:n.dat,".list",sep="")


  # Select the factorization (W and H) that leads to the lowest objective function 'f.WH'
  convergence <- convergence.list[[which.min(min.f.WH)]]
  W <- W.list[[which.min(min.f.WH)]]
  H <- NULL
  for (i in 1:n.dat){
    assign(paste("H",i,sep=""),H.list[[i]][[which.min(min.f.WH)]])
    H <- c(H,list(eval(parse(text=paste("H",i,sep="")))))
  }
  names(H) <- paste("H",1:n.dat,sep="")

  # Compute cluster membership
  clusters <- apply(W,1,which.max)
  # Output the following results
  output <- list(W=W,H=H,convergence=convergence,min.f.WH=min.f.WH,clusters=clusters)
  return(output)
}

ref_laplacian = function(spa_mat, sigma = 1){
  dst2=as.matrix(dist(spa_mat)^2)
  adj_matrix = exp(-1 * dst2 / (2*sigma^2))
  degree_matrix <- diag(rowSums(adj_matrix))
  list(laplacian = degree_matrix - adj_matrix, degree = degree_matrix, adj = adj_matrix)
}

# A small Visium-like hexagonal grid in SpaNiche's (row, col / sqrt(3)) units.
hex_grid = function(nrow = 12, ncol = 24) {
  g = expand.grid(row = 0:(nrow - 1), col = 0:(ncol - 1))
  g = g[(g$row + g$col) %% 2 == 0, ]
  g$col = g$col / sqrt(3)
  rownames(g) = sprintf("spot%03d", seq_len(nrow(g)))
  g
}
