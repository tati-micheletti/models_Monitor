#' Extract the three scales' suitability scores at occurrence locations
#'
#' @param spPa data.frame with `x`, `y` columns (one year's occurrence rows).
#' @param spClean Character. Species name with spaces replaced by underscores.
#' @param year Integer.
#' @param modelDirs Named list with `europe`, `landscape`, `habitat` prediction directories.
#' @param refRaster SpatRaster. Shared 200m reference grid.
#' @param idCols Character vector of ID/metadata columns from `spPa` to keep.
#' @return data.frame with `idCols`, the three `<scale>_mean_prob` columns,
#'   and `year`; or NULL if any scale's suitability is unavailable.
extractSuitability <- function(spPa, spClean, year, modelDirs, refRaster, idCols) {
  rClim <- loadSuitability("climate", spClean, year, modelDirs, refRaster)
  rLand <- loadSuitability("landscape", spClean, year, modelDirs, refRaster)
  rHab <- loadSuitability("habitat", spClean, year, modelDirs, refRaster)

  if (is.null(rClim) || is.null(rLand) || is.null(rHab)) {
    missing <- c(if (is.null(rClim)) "climate", if (is.null(rLand)) "landscape",
                 if (is.null(rHab)) "habitat")
    warning("Missing suitability for year ", year, ": ", paste(missing, collapse = ", "))
    return(NULL)
  }

  suitStack <- c(rClim, rLand, rHab)
  coords <- as.matrix(spPa[, c("x", "y")])
  envVals <- terra::extract(suitStack, coords)

  cbind(spPa[, idCols], envVals, year = year)
}
