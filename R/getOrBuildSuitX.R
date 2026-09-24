#' Get (from cache) or build the complete-case suitability matrix for one
#' species/year, for use as `bootstrapMetaModelTrend()`'s `newX`
#'
#' Building this matrix means `loadSuitability()`-ing and `resample()`-ing
#' all three scales onto the shared reference grid -- and that resample is
#' expensive at this raster's real extent (still the uncropped, full-Europe
#' bounding box; see improvements.md), on the order of tens of minutes
#' across a full prediction-year range for one species, confirmed by a real
#' timed dry run. `metaModel()`'s own prediction cache
#' (`isValidPredictionRaster()`) only skips re-*predicting*, not this
#' *input*-side reconstruction, so without a cache of its own, a cache-hit
#' year would silently pay the full resample cost again on every re-run of
#' the trend bootstrap (e.g. re-tuning `nBootTrend`, or a later rerun for a
#' newly-added species) -- exactly the kind of repeated, avoidable cost
#' `isValidRasterFile()`-style caching exists to prevent elsewhere in this
#' pipeline.
#'
#' @param spClean Character. Species name with spaces replaced by underscores.
#' @param year Integer.
#' @param modelDirs Named list with `europe`, `landscape`, `habitat` prediction directories.
#' @param refRaster SpatRaster. Shared 200m reference grid -- pass a version
#'   already cropped to (approximately) Germany's extent if you have one:
#'   `loadSuitability()`'s `resample()` cost scales with the target grid's
#'   size, and this reference grid is otherwise still the uncropped,
#'   full-Europe bounding box (~738M cells, ~1.9% real German data), which
#'   is what made the resample step take ~85 minutes for one species in a
#'   real timed dry run. Cropping to Germany first (a cheap extent
#'   subsetting operation, not an interpolation) before this function's own
#'   resample calls is the fix; it doesn't need to be an exact national
#'   mask, just a much smaller bounding box.
#' @param suitCols Character vector of the three suitability column names,
#'   in the order the ridge model was trained on.
#' @param outputDir Character. Directory to cache the extracted matrix in
#'   (the meta-model's own output directory).
#' @return Numeric matrix (cells x 3), or `NULL` if any of the three scales'
#'   suitability rasters is missing for this species/year.
getOrBuildSuitX <- function(spClean, year, modelDirs, refRaster, suitCols, outputDir) {
  cachePath <- file.path(outputDir, paste0(spClean, "_meta_suitX_", year, ".rds"))

  if (isValidCachedRDS(cachePath)) {
    return(readRDS(cachePath))
  }

  rClim <- loadSuitability("climate", spClean, year, modelDirs, refRaster)
  rLand <- loadSuitability("landscape", spClean, year, modelDirs, refRaster)
  rHab <- loadSuitability("habitat", spClean, year, modelDirs, refRaster)
  if (is.null(rClim) || is.null(rLand) || is.null(rHab)) return(NULL)

  suitStack <- c(rClim, rLand, rHab)
  names(suitStack) <- suitCols

  # na.rm = TRUE is essential, not cosmetic: this stack is still the
  # uncropped full-Europe extent (~738M cells, ~1.9% real German data), so
  # na.rm = FALSE would materialize a data.frame with ~98% NA rows just to
  # discard them a line later.
  predDf <- as.data.frame(suitStack, xy = FALSE, na.rm = TRUE)
  suitX <- as.matrix(predDf[, suitCols])

  saveRDS(suitX, cachePath)
  suitX
}
