#' Aggregate a species' area-mean occurrence probability to a coarser grid
#'
#' Same underlying data as `extractMetaProbSeries()`, but instead of
#' collapsing each year's map to one country-wide number, aggregates it to
#' a grid of `cellSizeM`-sided cells (each cell's value = mean `meta_prob`
#' of the 200m pixels inside it) -- the building block for a regional
#' (rather than national) index. See `computeGriddedCombinedIndex()`.
#'
#' @param species Character. Latin species name.
#' @param years Integer vector of years to extract.
#' @param metaDir Character. Directory holding metaModel()'s output.
#' @param cellSizeM Numeric. Target grid cell size in meters (e.g. 20000
#'   for a 20km grid). Must be a multiple of the source raster's own
#'   resolution or close to it; `terra::aggregate()`'s `fact` is rounded
#'   to the nearest whole number of source cells.
#' @return SpatRaster, one layer per year with real data (years with a
#'   missing/unreadable file are silently dropped, with a warning), named
#'   by year.
aggregateSpeciesToGrid <- function(species, years, metaDir, cellSizeM) {
  spClean <- gsub(" ", "_", species)

  layers <- lapply(years, function(yr) {
    f <- file.path(metaDir, paste0(spClean, "_meta_suitability_", yr, ".tif"))
    if (!file.exists(f)) return(NULL)
    r <- tryCatch(terra::rast(f), error = function(e) NULL)
    if (is.null(r) || !("meta_prob" %in% names(r))) return(NULL)
    fact <- max(1, round(cellSizeM / terra::res(r)[1]))
    terra::aggregate(r[["meta_prob"]], fact = fact, fun = "mean", na.rm = TRUE)
  })

  valid <- !sapply(layers, is.null)
  if (!all(valid)) {
    warning(species, ": missing/unreadable for years ", paste(years[!valid], collapse = ", "))
  }
  if (!any(valid)) return(NULL)

  stack <- terra::rast(layers[valid])
  names(stack) <- as.character(years[valid])
  stack
}
