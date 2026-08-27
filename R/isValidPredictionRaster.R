#' Check a saved prediction raster exists and isn't a silent all-zero failure
#'
#' Beyond basic file-openability, checks that `layerName` (when present)
#' has a non-zero maximum -- a BRT/ridge prediction that silently wrote
#' all-zero values is a bug, not a valid cached result, and should be
#' recomputed rather than treated as a cache hit.
#'
#' @param path Character. Path to a raster file.
#' @param layerName Character. Layer to sanity-check, default "mean_prob".
#' @return Logical.
isValidPredictionRaster <- function(path, layerName = "mean_prob") {
  if (!file.exists(path)) return(FALSE)
  tryCatch({
    r <- terra::rast(path)
    if (terra::nlyr(r) == 0) return(FALSE)
    if (layerName %in% names(r)) {
      return(terra::minmax(r[[layerName]])[2] > 0)
    }
    TRUE
  }, error = function(e) FALSE)
}
