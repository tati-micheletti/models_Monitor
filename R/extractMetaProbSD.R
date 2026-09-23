#' Extract a species' spatial standard error of the mean, per year
#'
#' Companion to `extractMetaProbSeries()` -- same files, same year loop, but
#' pulls the spatial standard deviation and valid-cell count of `meta_prob`
#' instead of just the mean, so a per-species-year variance is available for
#' the DDA/PECBMS-style analytical combined-index CI (see
#' `computeCombinedIndexAnalytical()`). Treats grid cells as independent
#' samples for `SE = sd / sqrt(n)` -- a simplification, since neighbouring
#' cells are spatially autocorrelated in reality, so this SE is a lower
#' bound on the true uncertainty, not an exact figure.
#'
#' @param species Character. Latin species name.
#' @param years Integer vector of years to extract.
#' @param metaDir Character. Directory holding metaModel()'s output.
#' @return Named numeric vector of SE(mean) values, one per year (name =
#'   as.character(year)); `NA` for any year whose file is missing/unreadable.
extractMetaProbSD <- function(species, years, metaDir) {
  spClean <- gsub(" ", "_", species)

  vals <- vapply(years, function(yr) {
    f <- file.path(metaDir, paste0(spClean, "_meta_suitability_", yr, ".tif"))
    if (!file.exists(f)) return(NA_real_)
    r <- tryCatch(terra::rast(f), error = function(e) NULL)
    if (is.null(r) || !("meta_prob" %in% names(r))) return(NA_real_)

    stats <- terra::global(r[["meta_prob"]], c("sd", "notNA"), na.rm = TRUE)
    sdVal <- stats[1, "sd"]
    n <- stats[1, "notNA"]
    if (is.na(sdVal) || is.na(n) || n < 2) return(NA_real_)
    sdVal / sqrt(n)
  }, numeric(1))

  names(vals) <- as.character(years)
  vals
}
