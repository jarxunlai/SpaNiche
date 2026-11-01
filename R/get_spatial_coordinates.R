#' Spot coordinates of a Visium Seurat object in SpaNiche's spatialdf format
#'
#' @description
#' Returns the array coordinates (\code{row}, \code{col}) of every spot, which
#' is the \code{spatialdf} input of the SpaNiche analysis functions. Works with
#' Seurat v4 objects (\code{VisiumV1} images, which store array coordinates)
#' and Seurat v5 objects (\code{VisiumV2} images, which only store pixel
#' centroids; the array coordinates are then read from the Space Ranger
#' output given in \code{data.dir}).
#'
#' @param object A Seurat object created with \code{Load10X_Spatial()}.
#' @param data.dir The Space Ranger output directory that was passed to
#'   \code{Load10X_Spatial()} (the folder containing \code{spatial/}). Only
#'   needed when the image does not store array coordinates (Seurat v5).
#' @param image Name of the image to use. Defaults to the first image.
#' @param hex_to_euclidean If \code{TRUE} (default), divide \code{col} by
#'   \eqn{\sqrt{3}} so that Euclidean distances on (row, col) match the
#'   hexagonal Visium layout, as expected by the default distance thresholds
#'   \code{2/sqrt(3)} and \code{4/sqrt(3)}.
#'
#' @return A data frame with one row per spot (row names are barcodes, in the
#'   order of \code{colnames(object)}) and columns \code{tissue}, \code{row},
#'   \code{col}, \code{imagerow}, \code{imagecol}.
#' @export
get_spatial_coordinates = function(
    object,
    data.dir = NULL,
    image = NULL,
    hex_to_euclidean = TRUE
){
  image = image %||% SeuratObject::Images(object)[1]
  if (is.na(image)) stop("The Seurat object contains no spatial image.")
  img = object[[image]]

  if (methods::.hasSlot(img, "coordinates") &&
      all(c("row", "col") %in% colnames(img@coordinates))) {
    coords = img@coordinates
  } else if (!is.null(data.dir)) {
    coords = read_tissue_positions(data.dir)
  } else {
    stop("Image '", image, "' (class ", class(img)[1], ") does not store array ",
         "row/col coordinates. Pass data.dir = the Space Ranger output folder ",
         "used with Load10X_Spatial().")
  }

  spots = intersect(colnames(object), SeuratObject::Cells(img))
  missing = setdiff(spots, rownames(coords))
  if (length(missing) > 0) {
    stop(length(missing), " spots of the image are missing from the tissue positions, e.g. ",
         paste(utils::head(missing, 3), collapse = ", "))
  }
  coords = coords[spots, c("tissue", "row", "col", "imagerow", "imagecol")]
  if (hex_to_euclidean) coords$col = coords$col / sqrt(3)
  coords
}

#' Read Space Ranger tissue positions
#'
#' @description
#' Reads the spot positions written by Space Ranger, whichever file format the
#' version produced: \code{tissue_positions_list.csv} (Space Ranger 1.x, no
#' header), \code{tissue_positions.csv} (Space Ranger >= 2.0, with header) or
#' \code{tissue_positions.parquet} (Visium HD, requires the \pkg{arrow} package).
#'
#' @param data.dir The Space Ranger output directory (containing
#'   \code{spatial/}), or the \code{spatial/} directory itself.
#'
#' @return A data frame with barcodes as row names and columns \code{tissue},
#'   \code{row}, \code{col}, \code{imagerow}, \code{imagecol}.
#' @export
read_tissue_positions = function(data.dir){
  candidates = file.path(
    rep(c(file.path(data.dir, "spatial"), data.dir), each = 3),
    c("tissue_positions.parquet", "tissue_positions.csv", "tissue_positions_list.csv")
  )
  file = candidates[file.exists(candidates)][1]
  if (is.na(file)) stop("No tissue_positions file found in ", data.dir)

  std_names = c("tissue", "row", "col", "imagerow", "imagecol")
  if (grepl("\\.parquet$", file)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Reading ", basename(file), " requires the arrow package.")
    }
    df = as.data.frame(arrow::read_parquet(file))
  } else {
    df = utils::read.csv(file, header = basename(file) == "tissue_positions.csv",
                         check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (ncol(df) < 6) stop("Unexpected format of ", file)
  # Space Ranger always writes: barcode, in_tissue, array_row, array_col,
  # pxl_row_in_fullres, pxl_col_in_fullres
  out = df[, 2:6]
  colnames(out) = std_names
  rownames(out) = df[[1]]
  out
}
