#' Load one scale's suitability raster for one species/year, resampled to
#' the shared 200m reference grid
#'
#' Every scale is unconditionally resampled onto `refRaster` -- habitat-scale
#' extent is not stable across years (the underlying land use source data
#' changes processing extent between releases), so habitat can no longer be
#' assumed to already match the reference grid just because it's nominally
#' "200m." The reference grid is a static DEM-derived raster (not a
#' species-specific prediction), so no cells within Germany are lost by this.
#'
#' @param scale Character. One of "climate", "landscape", "habitat".
#' @param spClean Character. Species name with spaces replaced by underscores.
#' @param year Integer.
#' @param modelDirs Named list with `europe`, `landscape`, `habitat` prediction directories.
#' @param refRaster SpatRaster. Shared 200m reference grid.
#' @return SpatRaster (one layer, named `<scale>_mean_prob`), or NULL if
#'   the source prediction file doesn't exist.
loadSuitability <- function(scale, spClean, year, modelDirs, refRaster) {
  f <- switch(scale,
              climate = file.path(modelDirs$europe, paste0(spClean, "_pred_EU_", year, ".tif")),
              landscape = file.path(modelDirs$landscape, paste0(spClean, "_pred_landscape_", year, ".tif")),
              habitat = file.path(modelDirs$habitat, paste0(spClean, "_pred_habitat_", year, ".tif")))

  if (!file.exists(f)) return(NULL)

  r <- tryCatch(terra::rast(f)[["mean_prob"]], error = function(e) NULL)
  if (is.null(r)) return(NULL)

  r <- terra::resample(r, refRaster, method = "bilinear")
  names(r) <- paste0(scale, "_mean_prob")
  r
}
