#' Extract a species' area-mean occurrence probability time series
#'
#' Reads `<species>_meta_suitability_<year>.tif` (the `metaModel()` output,
#' `meta_prob` + `binary` bands) for one species across a year range and
#' collapses each year's raster to a single scalar: the spatial mean of
#' `meta_prob` across Germany. This scalar plays the role of a
#' Zbinden-et-al.-style population index value for that species/year --
#' built from occurrence probability rather than a count-based abundance
#' index, since that's what this pipeline's meta-model produces.
#'
#' @param species Character. Latin species name (with spaces, as used
#'   elsewhere in the pipeline).
#' @param years Integer vector of years to extract.
#' @param metaDir Character. Directory holding metaModel()'s output
#'   (`outputs/<runName>/models_Monitor/meta`).
#' @return Named numeric vector, one value per year in `years` (name = as.character(year));
#'   `NA` for any year whose file is missing or unreadable.
extractMetaProbSeries <- function(species, years, metaDir) {
  spClean <- gsub(" ", "_", species)

  vals <- vapply(years, function(yr) {
    f <- file.path(metaDir, paste0(spClean, "_meta_suitability_", yr, ".tif"))
    if (!file.exists(f)) return(NA_real_)
    r <- tryCatch(terra::rast(f), error = function(e) NULL)
    if (is.null(r) || !("meta_prob" %in% names(r))) return(NA_real_)
    as.numeric(terra::global(r[["meta_prob"]], "mean", na.rm = TRUE)[1, 1])
  }, numeric(1))

  names(vals) <- as.character(years)
  vals
}
