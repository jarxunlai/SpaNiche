# The optimised implementations must reproduce the original dense ones.

test_that("make_distance_matrix (sparse) equals the original dense version", {
  g = hex_grid()
  for (thre in list(c(2 / sqrt(3), 4 / sqrt(3)), c(1.5, 2.2), c(1, 3))) {
    new = make_distance_matrix(g, c("view1", "view2"), thre)
    ref = ref_make_distance_matrix(g, c("view1", "view2"), thre)
    expect_s4_class(new$view1, "dgCMatrix")
    for (v in c("view1", "view2")) {
      expect_equal(as.matrix(new[[v]]), ref[[v]], tolerance = 1e-14)
    }
  }
  # single view
  new = make_distance_matrix(g, "view1", 2 / sqrt(3))
  ref = ref_make_distance_matrix(g, "view1", 2 / sqrt(3))
  expect_equal(as.matrix(new$view1), ref$view1, tolerance = 1e-14)
})

test_that("make_distance_matrix_LR (sparse) equals the original dense version", {
  g = hex_grid()
  for (thre in list(c(2 / sqrt(3), 4 / sqrt(3)), c(1.5, 2.2), c(1.2, 2))) {
    new = make_distance_matrix_LR(g, c("view1", "view2"), thre)
    ref = ref_make_distance_matrix_LR(g, c("view1", "view2"), thre)
    for (v in c("view1", "view2")) {
      expect_equal(as.matrix(new[[v]]), ref[[v]], tolerance = 1e-14)
    }
  }
  # thre1 < 1: view1 is empty, and the recoding puts each spot itself into view2
  expect_warning(new <- make_distance_matrix_LR(g, c("view1", "view2"), c(0.9, 2)), "no neighbours")
  ref = ref_make_distance_matrix_LR(g, c("view1", "view2"), c(0.9, 2))
  expect_equal(as.matrix(new$view2), ref$view2, tolerance = 1e-14)
  expect_true(all(diag(as.matrix(new$view2)) > 0))
})

test_that("spots without neighbours give zero rows with a warning (dense code gave NaN)", {
  g = hex_grid()
  # thresholds below 1: the dense recoding (<= t -> 1, then > t -> 0) empties every row
  expect_warning(new <- make_distance_matrix(g, "view1", 0.5), "no neighbours")
  ref = ref_make_distance_matrix(g, "view1", 0.5)
  expect_true(all(is.nan(ref$view1)))
  expect_true(all(as.matrix(new$view1) == 0))
})

test_that("neighbour search handles rounding at the threshold and block edges", {
  set.seed(1)
  pts = matrix(runif(2 * 3000, 0, 40), ncol = 2)
  rownames(pts) = paste0("p", 1:3000)
  pairs = SpaNiche:::.neighbor_pairs(pts, 1.15, digits = 2, block = 257L)
  d = round(as.matrix(dist(pts)), 2)
  expected = which(d <= 1.15, arr.ind = TRUE)
  expect_equal(nrow(pairs), nrow(expected))
  expect_setequal(paste(pairs$i, pairs$j), paste(expected[, 1], expected[, 2]))
})

test_that("generate_LRintegratedmatrix gives the same values with a sparse neighbour matrix", {
  g = hex_grid()
  n = nrow(g)
  set.seed(2)
  genes = c("A", "B", "C", "D")
  expr = matrix(rexp(n * 4), n, dimnames = list(rownames(g), genes))
  receptors = c("A", "B", "A", "C", "A")  # repeated genes, as in real LR tables
  ligands = c("D", "D", "B", "C", "C")
  rmat = expr[, receptors]
  lmat = expr[, ligands]
  nb = make_distance_matrix_LR(g, c("view1", "view2"), c(2 / sqrt(3), 4 / sqrt(3)))
  ref_nb = ref_make_distance_matrix_LR(g, c("view1", "view2"), c(2 / sqrt(3), 4 / sqrt(3)))
  out = suppressMessages(capture.output(
    new <- generate_LRintegratedmatrix(nb$view2, rmat, lmat, center_is_ligand = TRUE)
  ))
  ref = ref_generate_LRintegratedmatrix(ref_nb$view2, rmat, lmat)
  expect_equal(unname(new), unname(ref), tolerance = 1e-12)
  expect_equal(colnames(new), paste0(ligands, "->", receptors, " (CisL, NisR)"))
})

test_that("get_spot_by_celltype_extended returns the original dense result", {
  g = hex_grid()
  set.seed(3)
  ct = matrix(runif(nrow(g) * 3), ncol = 3, dimnames = list(rownames(g), c("T", "B", "M")))
  new = get_spot_by_celltype_extended(g, ct)
  ref_dm = ref_make_distance_matrix(g, c("view1", "view2"), c(2 / sqrt(3), 4 / sqrt(3)))
  ref = cbind(ct, ref_dm$view1 %*% ct, ref_dm$view2 %*% ct)
  expect_true(is.matrix(new))
  expect_equal(unname(new), unname(ref), tolerance = 1e-12)
  expect_equal(colnames(new), c(colnames(ct), paste0(colnames(ct), "_view1"), paste0(colnames(ct), "_view2")))
})

test_that("nmf_modified gives the same factorisation as the original implementation", {
  g = hex_grid(10, 20)
  n = nrow(g)
  k = 4
  set.seed(4)
  dat = list(matrix(rexp(n * 6), n), matrix(rexp(n * 30), n))
  Wzero = matrix(runif(n * k), n)
  Hzero = list(H1 = matrix(runif(k * 6), k), H2 = matrix(runif(k * 30), k))
  args = list(dat = dat, Wzero = Wzero, Hzero = Hzero, wt = c(1, 0.5), k = k,
              lambda = c(0.5, 2), sigma = c(0.5, 1), maxiter = 60, st.count = 5,
              spatialdf = as.matrix(g[, c("row", "col")]), epsilon = 1e-4)
  new = do.call(nmf_modified, args)
  ref = do.call(ref_nmf_modified, args)
  expect_equal(new$W, ref$W, tolerance = 1e-8)
  expect_equal(new$H, ref$H, tolerance = 1e-8)
  expect_equal(new$min.f.WH, ref$min.f.WH, tolerance = 1e-10)
  expect_equal(new$convergence$iter, ref$convergence$iter)
  expect_equal(new$convergence$f.WH, ref$convergence$f.WH, tolerance = 1e-10)
  expect_equal(new$clusters, ref$clusters)
})

test_that("vectorised Moran's I matches ape::Moran.I with Seurat's weights", {
  skip_if_not_installed("ape")
  g = hex_grid(8, 16)
  n = nrow(g)
  set.seed(5)
  smooth = outer(sin(g$row / 2), rep(1, 6)) + matrix(rnorm(n * 6, sd = 0.3), n)
  x = t(cbind(smooth, matrix(rexp(n * 6), n), 1))  # last feature is constant
  rownames(x) = paste0("f", 1:13)
  pos = as.matrix(g[, c("row", "col")]) * 37.5  # coordinates' scale does not matter
  res = SpaNiche:::.moransi_svf(x, pos, nfeatures = 5)

  w = 1 / as.matrix(dist(pos))^2
  diag(w) = 0
  for (f in 1:12) {
    a = ape::Moran.I(x[f, ], w)
    expect_equal(res$MoransI_observed[f], a$observed, tolerance = 1e-12)
    expect_equal(res$MoransI_p.value[f], a$p.value, tolerance = 1e-10)
  }
  expect_true(is.na(res$MoransI_observed[13]))
  ord = order(res$MoransI_p.value[1:12], -abs(res$MoransI_observed[1:12]))
  expect_equal(res$moransi.spatially.variable.rank[ord], 1:12)
  expect_equal(sum(res$moransi.spatially.variable, na.rm = TRUE), 5)
})

test_that(".same_partition is equivalent to comparing co-membership matrices", {
  conn = function(x) outer(x, x, "==")
  set.seed(6)
  for (r in 1:200) {
    a = sample(1:4, 12, replace = TRUE)
    b = if (r %% 3 == 0) c(5, 9, 7, 2)[a] else sample(1:4, 12, replace = TRUE)  # relabelled copy or random
    expect_identical(SpaNiche:::.same_partition(a, b), all(conn(a) == conn(b)))
  }
  expect_false(SpaNiche:::.same_partition(1:3, NULL))
})

test_that(".nmf_mnnals reproduces IntNMF::nmf.mnnals exactly", {
  set.seed(7)
  n = 120
  dat = list(matrix(rexp(n * 8), n), matrix(rexp(n * 40), n))
  for (args in list(
    list(k = 4, maxiter = 50, st.count = 10, n.ini = 3, ini.nndsvd = TRUE, wt = c(1, 0.3)),
    list(k = 3, maxiter = 40, st.count = 3, n.ini = 4, ini.nndsvd = FALSE, wt = c(2, 1))
  )) {
    args = c(list(dat = dat, seed = TRUE), args)
    new = do.call(SpaNiche:::.nmf_mnnals, args)
    ref = do.call(IntNMF::nmf.mnnals, args)
    expect_identical(new$W, ref$W)
    expect_identical(new$H, ref$H)
    expect_identical(new$convergence, ref$convergence)
    expect_identical(new$min.f.WH, ref$min.f.WH)
    expect_identical(new$clusters, ref$clusters)
    expect_true(any(ref$convergence[, "stability"] == 1))  # the stopping rule was exercised
  }
})
