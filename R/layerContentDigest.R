#' Content digest of selected layers of a covariate stack (species-level cache keys)
#'
#' `reproducible::Cache()` keys on the CONTENT of a raster, but a file-backed subset of layers is digested as the whole file (tested 2026-10-21: changing an unused layer of the file
#' re-triggered the function), so subsetting a stack is not enough to make a species' cache independent of layers it does not use. Instead the function is called with
#' `omitArgs = "covStack"` and `.cacheExtra = list(covLayers = layerContentDigest(covStack, layers))`: the key then depends on the content of exactly the layers the species uses (tested:
#' an unused layer changes nothing, a used one recomputes only the species that use it).
#'
#' If `file` (the file the stack was read from) is given, the digests of ALL its layers are computed once and kept in a sidecar `<file>.layerdigests.rds`, valid as long as size and
#' modification time of the file are the same, so repeated calls (every species, every task) do not read the layers again.
#'
#' @param covStack SpatRaster (names = layer names without year suffix).
#' @param layers Character. Layers to include (names missing in `covStack` are an error).
#' @param file Character or NULL. File of the stack, for the sidecar.
#' @return One character string (the digests of the layers, in the order of `layers`).
layerContentDigest <- function(covStack, layers, file = NULL) {
  layers <- setdiff(layers, c("x", "y"))
  idx <- layerIndex(covStack, layers)
  dg <- function(i) digest::digest(terra::values(covStack[[i]]), algo = "xxhash64")
  if (!is.null(file) && file.exists(file)) {
    side <- paste0(file, ".layerdigests.rds"); info <- file.info(file); stamp <- paste(info$size, as.numeric(info$mtime))
    tab <- if (file.exists(side)) tryCatch(readRDS(side), error = function(e) NULL) else NULL
    if (is.null(tab) || !identical(tab$stamp, stamp) || !identical(tab$names, names(covStack))) {
      tab <- list(stamp = stamp, names = names(covStack), digests = stats::setNames(vapply(seq_along(names(covStack)), dg, character(1)), names(covStack)))
      try(saveRDS(tab, side), silent = TRUE)
    }
    return(paste(unname(tab$digests[idx]), collapse = ""))
  }
  paste(vapply(idx, dg, character(1)), collapse = "")
}

#' Positions of `layers` in a covariate stack. The layer names of some stacks carry the year of their source ("grassland_2010", "built_up_2006"; the landscape stack of
#' `loadCovariates()`), the predictor names of `speciesConfig_predictors.csv` do not: a name is matched exactly first, then without a trailing `_<year>`.
layerIndex <- function(covStack, layers) {
  nm <- names(covStack); base <- gsub("_[0-9]{4}$", "", nm)
  idx <- match(layers, nm); noHit <- is.na(idx); idx[noHit] <- match(layers[noHit], base)
  if (anyNA(idx)) stop("layer(s) not in the covariate stack: ", paste(layers[is.na(idx)], collapse = ", "), " (stack: ", paste(nm, collapse = ", "), ")", call. = FALSE)
  idx
}

#' The stack of one species: only the layers it uses
#'
#' @param covStack SpatRaster of all layers. @param layers Character vector of the species' predictors (from `speciesConfig_predictors.csv`, scale column), or NULL (= all layers).
#' @param sp,scale Character, for the message. @return SpatRaster of the species' layers (x and y are coordinates, not layers).
speciesCovStack <- function(covStack, layers, sp = "", scale = "") {
  if (is.null(layers)) return(covStack)
  layers <- setdiff(layers, c("x", "y"))
  idx <- tryCatch(layerIndex(covStack, layers), error = function(e)
    stop("[PREDICTOR ERROR] ", sp, ": ", conditionMessage(e), " -- requested at the ", scale, " scale (speciesConfig_predictors.csv)", call. = FALSE))
  covStack[[idx]]
}
