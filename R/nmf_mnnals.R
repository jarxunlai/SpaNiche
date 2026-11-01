# Initial integrative NMF, numerically identical to IntNMF::nmf.mnnals()
# (same seed, same sequence of random numbers, same NNLS solvers from IntNMF,
# same stopping rule), returning the W, H, convergence and min.f.WH that
# SpaNiche uses.
#
# IntNMF::nmf.mnnals() also builds n x n connectivity and consensus matrices
# in every iteration; for a few thousand spots that bookkeeping dominates its
# run time. Here the stopping rule "the connectivity matrix did not change"
# is evaluated on the cluster labels instead, and no consensus matrix is
# returned.
.nmf_mnnals = function(dat, k, maxiter = 200, st.count = 20, n.ini = 30,
                       ini.nndsvd = TRUE, seed = TRUE, wt = rep(1, length(dat))) {
  if (seed) set.seed(12345)
  n.dat = length(dat)
  n = nrow(dat[[1]])
  if (n.dat != length(wt)) stop("Number of weights must match number of data")
  for (i in 1:n.dat) {
    if (!all(dat[[i]] >= 0)) stop("All values must be positive. There are -ve entries in dat", i)
  }
  nndsvd = utils::getFromNamespace(".nndsvd.internal", "IntNMF")
  H_fcnnls = utils::getFromNamespace("H.fcnnls", "IntNMF")
  W_fcnnls = utils::getFromNamespace("W.fcnnls", "IntNMF")

  min.f.WH = NULL
  best = NULL
  for (j in 1:n.ini) {
    abs.diff = NA
    if (j == 1 && ini.nndsvd) {
      tmp.H = lapply(dat, function(d) nndsvd(d, k, flag = 0)$H)
      W = t(W_fcnnls(x = tmp.H, y = dat, weight = wt)$coef)
    } else {
      W = matrix(stats::runif(n * k, min = min(unlist(dat)), max = max(unlist(dat))), nrow = n)
    }

    clust.old = NULL
    count = 1
    iter = 1
    convergence = NULL
    while ((iter < maxiter) && (count < st.count)) {
      W = sweep(W, 2, pmax(sqrt(colSums(W^2)), .Machine$double.eps), "/")
      H = lapply(dat, function(d) H_fcnnls(W, d)$coef)
      W = t(W_fcnnls(x = H, y = dat, weight = wt)$coef)

      WH.new = lapply(1:n.dat, function(i) W %*% H[[i]])
      if (iter > 1) {
        abs.diff = 0
        for (i in 1:n.dat) abs.diff = abs.diff + sum(abs(WH.new[[i]] - WH.old[[i]])) / sum(WH.old[[i]])
      }
      WH.old = WH.new
      f.WH = 0
      for (i in 1:n.dat) f.WH = f.WH + wt[i] * sum((dat[[i]] - WH.new[[i]])^2)
      f.WH = sqrt(f.WH)

      clust.mem = apply(W, 1, which.max)
      stable = .same_partition(clust.mem, clust.old)
      if (stable) count = count + 1 else count = 0
      convergence = rbind(convergence, c(iter, count, ifelse(stable, 1, 0), abs.diff, f.WH))
      clust.old = clust.mem
      iter = iter + 1
    }
    colnames(convergence) = c("iter", "count", "stability", "abs.diff", "f.WH")
    min.f.WH = c(min.f.WH, f.WH)
    # which.min() keeps the first of tied minima, so only a strictly smaller value replaces
    if (is.null(best) || f.WH < best$f.WH) {
      best = list(W = W, H = H, convergence = convergence, f.WH = f.WH)
    }
  }

  names(best$H) = paste("H", 1:n.dat, sep = "")
  list(W = best$W, H = best$H, convergence = best$convergence,
       min.f.WH = min.f.WH, clusters = apply(best$W, 1, which.max))
}

# TRUE when two labelings define the same partition, i.e. the same
# co-membership (connectivity) matrix: every class of `a` maps to exactly one
# class of `b` and vice versa.
.same_partition = function(a, b) {
  if (is.null(b) || length(a) != length(b)) return(FALSE)
  n_pairs = nrow(unique(cbind(a, b)))
  n_pairs == length(unique(a)) && n_pairs == length(unique(b))
}
