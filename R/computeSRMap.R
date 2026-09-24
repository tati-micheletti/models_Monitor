#' Species richness map for one year, from meta-model outputs
#'
#' Two variants, both computed (pick whichever fits the audience):
#' - `SR_prob`: sum of per-species occurrence probability across species --
#'   no threshold needed, propagates each species' own uncertainty
#'   directly (the "expected richness" approach standard in the stacked-SDM
#'   richness literature).
#' - `SR_binary`: sum of each species' already-thresholded `binary` band
#'   (from `metaModel()`/`predictRidgeToRaster()`, using each species' own
#'   CV-optimized threshold) -- the more traditional "number of species
#'   predicted present" map.
#'
#' @param species Character vector of Latin species names to include.
#' @param year Integer. Year to compute the richness map for.
#' @param metaDir Character. Directory holding metaModel()'s output.
#' @return A 2-layer SpatRaster (`SR_prob`, `SR_binary`), or `NULL` if no
#'   species had a valid file for `year`.
computeSRMap <- function(species, year, metaDir) {
  probLayers <- list()
  binLayers <- list()

  for (sp in species) {
    spClean <- gsub(" ", "_", sp)
    f <- file.path(metaDir, paste0(spClean, "_meta_suitability_", year, ".tif"))
    if (!file.exists(f)) next
    r <- tryCatch(terra::rast(f), error = function(e) NULL)
    if (is.null(r)) next
    if ("meta_prob" %in% names(r)) probLayers[[sp]] <- r[["meta_prob"]]
    if ("binary" %in% names(r)) binLayers[[sp]] <- r[["binary"]]
  }

  if (length(probLayers) == 0) {
    warning("No valid species rasters found for year ", year, " -- returning NULL")
    return(NULL)
  }

  srProb <- sum(terra::rast(probLayers), na.rm = TRUE)
  names(srProb) <- "SR_prob"

  srBinary <- sum(terra::rast(binLayers), na.rm = TRUE)
  names(srBinary) <- "SR_binary"

  c(srProb, srBinary)
}
