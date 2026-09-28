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
  terra::writeRaster(combineTwoLayerRaster(layer1, layer2, layerNames), outPath, overwrite = TRUE)
  invisible(outPath)
}

#' Combine two single-layer rasters into one two-layer SpatRaster in memory,
#' working around the same terra stacking bug `writeTwoLayerRaster()` does
#'
#' Used where the caller needs the combined raster OBJECT (e.g. to pass
#' through `reproducible::Cache()`, whose cached return value has no fixed
#' output path of its own) rather than writing straight to a final path --
#' see `predictBRTToRaster()`. Forces the combined raster fully into memory
#' (`values<-` with a full replacement) before deleting the temp files it
#' was built from, so the returned object stays valid afterward.
#'
#' @param layer1 SpatRaster, first layer.
#' @param layer2 SpatRaster, second layer.
#' @param layerNames Character vector of length 2, output layer names.
#' @return SpatRaster, two layers, in-memory.
combineTwoLayerRaster <- function(layer1, layer2, layerNames) {
  tmp1 <- tempfile(fileext = ".tif")
  tmp2 <- tempfile(fileext = ".tif")
  terra::writeRaster(layer1, tmp1, overwrite = TRUE)
  terra::writeRaster(layer2, tmp2, overwrite = TRUE)

  rOut <- c(terra::rast(tmp1), terra::rast(tmp2))
  names(rOut) <- layerNames
  terra::values(rOut) <- terra::values(rOut)
  file.remove(tmp1, tmp2)

  rOut
}
