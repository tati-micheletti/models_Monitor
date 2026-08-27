#' Write two single-layer rasters to one two-layer file, working around a
#' terra stacking bug
#'
#' Stacking a derived raster (e.g. `prob >= threshold`) together with its
#' source raster and writing directly causes the source layer's values
#' to zero out in some terra versions. Writing each layer to its own
#' temp file first, then re-reading and combining, avoids it.
#'
#' @param layer1 SpatRaster, first layer.
#' @param layer2 SpatRaster, second layer.
#' @param layerNames Character vector of length 2, output layer names.
#' @param outPath Character. Final output file path.
#' @return Invisibly, `outPath`.
writeTwoLayerRaster <- function(layer1, layer2, layerNames, outPath) {
  tmp1 <- tempfile(fileext = ".tif")
  tmp2 <- tempfile(fileext = ".tif")
  terra::writeRaster(layer1, tmp1, overwrite = TRUE)
  terra::writeRaster(layer2, tmp2, overwrite = TRUE)

  rOut <- c(terra::rast(tmp1), terra::rast(tmp2))
  names(rOut) <- layerNames
  terra::writeRaster(rOut, outPath, overwrite = TRUE)
  file.remove(tmp1, tmp2)

  invisible(outPath)
}
