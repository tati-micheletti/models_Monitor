#' Determine the spatial CV block size for one species
#'
#' Uses `blockCV::cv_spatial_autocor()` to estimate the natural spatial
#' autocorrelation range, then caps/floors it to sane bounds. Follows
#' Wiedenroth et al.
#'
#' @param sfOcc `sf` object of presence/absence points.
#' @param maxBlockSizeM Numeric. Maximum block size in metres.
#' @param minBlockSizeM Numeric. Minimum block size in metres (default 0,
#'   i.e. no floor -- pass a floor explicitly where the original scripts
#'   used one, e.g. European scale).
#' @return Numeric, the final block size in metres.
determineBlockSize <- function(sfOcc, maxBlockSizeM, minBlockSizeM = 0) {
  cvBlocksize <- tryCatch({
    blockCV::cv_spatial_autocor(x = sfOcc, column = "occurrence", plot = FALSE)
  }, error = function(e) {
    warning("cv_spatial_autocor failed: ", e$message, " -- using maxBlockSizeM")
    list(range = maxBlockSizeM)
  })

  cvBlSize <- cvBlocksize$range

  if (cvBlSize > maxBlockSizeM) {
    message("Block size (", round(cvBlSize / 1000), "km) exceeds cap -- capping at ",
            maxBlockSizeM / 1000, "km")
    cvBlSize <- maxBlockSizeM
  }
  if (cvBlSize < minBlockSizeM) {
    message("Block size (", round(cvBlSize / 1000), "km) below minimum -- setting to ",
            minBlockSizeM / 1000, "km")
    cvBlSize <- minBlockSizeM
  }

  message("Final block size: ", round(cvBlSize / 1000), "km")
  cvBlSize
}
