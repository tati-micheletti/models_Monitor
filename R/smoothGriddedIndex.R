#' Spatially smooth a gridded index (blend cell edges into a gradient)
#'
#' Applies a Gaussian-weighted moving-window average across neighbouring
#' cells, one year-layer at a time (`terra::focal()` applies to every
#' layer of a multi-layer SpatRaster automatically). This is the spatial
#' analogue of `smoothIndex()` (which smooths a single national series
#' over time); here the smoothing is over space, within one year --
#' useful for showing a regional pattern as a continuous gradient rather
#' than blocky grid cells, which reads less like "this specific area is
#' called out" and more like an environmental gradient (see the project
#' discussion on regional maps that don't point at specific
#' administrative regions).
#'
#' @param griddedIndex SpatRaster from `computeGriddedCombinedIndex()`
#'   (or `aggregateSpeciesToGrid()`), one or more layers.
#' @param smoothRadiusFactor Numeric. The Gaussian kernel's standard
#'   deviation, as a multiple of the raster's own cell size. Larger =
#'   smoother/blurrier. Default 1.5.
#' @return SpatRaster, same dimensions/layers as `griddedIndex`, smoothed.
#'   Renormalized near NA boundaries (coastline/edge cells with a partial
#'   neighborhood) -- confirmed necessary: without it, `focal()`'s Gaussian
#'   weights no longer sum to 1 once NA neighbors are dropped, biasing
#'   edge cells toward zero and pushing the smoothed range outside the
#'   raw data's actual range, which a real weighted average never should.
smoothGriddedIndex <- function(griddedIndex, smoothRadiusFactor = 1.5) {
  cellSize <- terra::res(griddedIndex)[1]
  w <- terra::focalMat(griddedIndex, d = smoothRadiusFactor * cellSize, type = "Gauss")

  weightedSum <- terra::focal(griddedIndex, w = w, fun = "sum", na.rm = TRUE, na.policy = "omit")
  availableWeight <- terra::focal(!is.na(griddedIndex), w = w, fun = "sum", na.rm = TRUE, na.policy = "omit")

  weightedSum / availableWeight
}
