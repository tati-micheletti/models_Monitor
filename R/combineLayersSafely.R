#' Safely combine several single-layer SpatRasters into one multi-layer raster
#'
#' Generalizes `models_Monitor`'s `combineTwoLayerRaster()` (2-layer case) to
#' N layers: directly `c()`-ing several DERIVED (computed, not read-from-
#' disk) rasters can silently zero out earlier layers in some terra
#' versions/session states -- see `writeTwoLayerRaster.R`'s own docstring.
#' Writing each layer to its own temp file first, then re-reading and
#' combining, avoids it. Found via a real failure chaining
#' `computeChangeMaps()` into the same SpaDES session as `models_Monitor`
#' (it worked fine called standalone in a fresh R session, which is how
#' the pre-existing root-level `runIndex.R` script always ran it -- see
#' DECISIONS.md's 2026-09-28 "runIndex_Monitor" entry).
#'
#' NOTE: this exact function is deliberately duplicated verbatim across
#' modules (same filename) so each has no load-time dependency on another
#' module being loaded. Keep all copies byte-identical --
#' tools/check_duplicated_functions.R checks this.
#'
#' @param layers List of single-layer SpatRasters, each already named via
#'   `names(r) <- "..."`.
#' @return SpatRaster, all layers, in-memory, names taken from each input
#'   layer's own `names()`.
combineLayersSafely <- function(layers) {
  layerNames <- vapply(layers, names, character(1))

  tmpFiles <- vapply(layers, function(r) {
    f <- tempfile(fileext = ".tif")
    terra::writeRaster(r, f, overwrite = TRUE)
    f
  }, character(1))

  rOut <- terra::rast(tmpFiles)
  names(rOut) <- layerNames
  terra::values(rOut) <- terra::values(rOut)
  file.remove(tmpFiles)

  rOut
}
