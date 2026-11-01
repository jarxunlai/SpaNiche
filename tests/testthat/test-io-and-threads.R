positions_df = function() {
  data.frame(
    barcode = c("AAA-1", "CCC-1", "GGG-1"),
    in_tissue = c(1L, 1L, 0L),
    array_row = c(0L, 1L, 2L),
    array_col = c(0L, 1L, 2L),
    pxl_row_in_fullres = c(100L, 200L, 300L),
    pxl_col_in_fullres = c(110L, 210L, 310L)
  )
}

expect_positions = function(pos) {
  expect_equal(rownames(pos), c("AAA-1", "CCC-1", "GGG-1"))
  expect_equal(colnames(pos), c("tissue", "row", "col", "imagerow", "imagecol"))
  expect_equal(pos$row, c(0, 1, 2))
  expect_equal(pos$imagecol, c(110, 210, 310))
}

test_that("read_tissue_positions reads Space Ranger 1.x (no header) csv", {
  d = file.path(withr::local_tempdir(), "spatial")
  dir.create(d)
  write.table(positions_df(), file.path(d, "tissue_positions_list.csv"),
              sep = ",", row.names = FALSE, col.names = FALSE)
  expect_positions(read_tissue_positions(dirname(d)))
})

test_that("read_tissue_positions reads Space Ranger >= 2.0 csv with header", {
  d = file.path(withr::local_tempdir(), "spatial")
  dir.create(d)
  write.csv(positions_df(), file.path(d, "tissue_positions.csv"), row.names = FALSE)
  expect_positions(read_tissue_positions(dirname(d)))
  expect_positions(read_tissue_positions(d))  # the spatial/ folder itself
})

test_that("read_tissue_positions reads Visium HD parquet", {
  skip_if_not_installed("arrow")
  d = file.path(withr::local_tempdir(), "spatial")
  dir.create(d)
  arrow::write_parquet(positions_df(), file.path(d, "tissue_positions.parquet"))
  expect_positions(read_tissue_positions(dirname(d)))
})

test_that(".set_threads limits BLAS threads and .restore_threads restores them", {
  before = RhpcBLASctl::blas_get_num_procs()
  old = SpaNiche:::.set_threads(1L)
  expect_equal(RhpcBLASctl::blas_get_num_procs(), 1L)
  SpaNiche:::.restore_threads(old)
  expect_equal(RhpcBLASctl::blas_get_num_procs(), before)
  expect_null(SpaNiche:::.set_threads(NULL))
  expect_error(SpaNiche:::.set_threads(0), "positive integer")
})

test_that(".get_layer_data works for the installed Seurat object version", {
  m = Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 2, 3), x = c(5, 3, 1), dims = c(3, 3),
                           dimnames = list(c("g1", "g2", "g3"), c("c1", "c2", "c3")))
  obj = SeuratObject::CreateSeuratObject(m, assay = "Spatial")
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    expect_error(SpaNiche:::.get_layer_data(obj, "data"), "NormalizeData")
  }
  obj = suppressWarnings(Seurat::NormalizeData(obj, verbose = FALSE))
  d = SpaNiche:::.get_layer_data(obj, "data")
  expect_equal(dim(d), c(3L, 3L))
  expect_equal(unname(d["g1", "c1"]), log1p(5 / 5 * 1e4))
  expect_equal(SpaNiche:::.resolve_assay(obj), "Spatial")
})
