#' Build regional (gridded + smoothed) combined-index maps at several
#' resolutions, for direct comparison
#'
#' Ties together `computeGriddedCombinedIndex()` and `smoothGriddedIndex()`
#' for each requested cell size, writing both the raw gridded raster and
#' its smoothed version for each -- so e.g. a 20km and a 50km version can
#' be compared side by side, raw and smoothed, rather than committing to
#' one resolution.
#'
#' @param species Character vector of Latin species names.
#' @param years Integer vector of years.
#' @param baselineYear Integer. Index baseline.
#' @param metaDir Character. `metaModel()`'s output directory.
#' @param outputDir Character. Directory to save outputs in.
#' @param cellSizesM Numeric vector. Grid cell sizes to compute, in meters
#'   (default 10km, 20km and 50km).
#' @param smoothRadiusFactor Numeric. Passed to `smoothGriddedIndex()`.
#' @param countryBoundary SpatVector or NULL, passed through to
#'   `computeGriddedCombinedIndex()`/`aggregateSpeciesToGrid()`. Strongly
#'   recommended (e.g. GADM level-0 Germany, as `tools/maskToGermany.R`
#'   already fetches) -- without it, cells straddling the border can
#'   quietly average in non-German source data; see
#'   `aggregateSpeciesToGrid()`'s docstring.
#' @return Invisibly, a named list (by cell size, as character e.g.
#'   `"20000"`) of `list(raw = SpatRaster, smoothed = SpatRaster)`.
#'   Everything is also written to `outputDir` as
#'   `regional_index_<cellSizeKm>km_raw.tif` /
#'   `regional_index_<cellSizeKm>km_smoothed.tif`.
computeRegionalIndex <- function(species, years, baselineYear, metaDir, outputDir,
                                  cellSizesM = c(10000, 20000, 50000), smoothRadiusFactor = 1.5,
                                  countryBoundary = NULL) {
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  result <- list()

  for (cellSizeM in cellSizesM) {
    cellKm <- cellSizeM / 1000
    message("Computing regional index at ", cellKm, "km resolution...")

    raw <- computeGriddedCombinedIndex(species, years, baselineYear, metaDir, cellSizeM,
                                        countryBoundary = countryBoundary)
    message("  Smoothing...")
    smoothed <- smoothGriddedIndex(raw, smoothRadiusFactor = smoothRadiusFactor)

    terra::writeRaster(raw, file.path(outputDir, paste0("regional_index_", cellKm, "km_raw.tif")),
                        overwrite = TRUE)
    terra::writeRaster(smoothed, file.path(outputDir, paste0("regional_index_", cellKm, "km_smoothed.tif")),
                        overwrite = TRUE)
    message("  Saved -> regional_index_", cellKm, "km_{raw,smoothed}.tif")

    result[[as.character(cellSizeM)]] <- list(raw = raw, smoothed = smoothed)
  }

  invisible(result)
}
