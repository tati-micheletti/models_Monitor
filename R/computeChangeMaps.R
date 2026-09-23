#' Prevalence-change and gain/loss maps between two years
#'
#' Per Zbinden et al. (2005)'s own caution -- a combined/community index
#' can mask species-specific heterogeneity (one species' strong increase
#' can offset another's strong decline), so per-species indices must
#' always be examined alongside it, never the combined index in
#' isolation -- this returns both per-species maps and their community
#' aggregate, for every comparison.
#'
#' Three map types, all cell-by-cell raster algebra on existing
#' `metaModel()` output -- no new model fitting:
#' - **Prevalence change** (continuous, no threshold): `deltaP = P_current - P_ref`.
#' - **Gain/loss** (categorical, from the already-thresholded `binary`
#'   band): 0 = stable absent, 1 = loss, 2 = gain, 3 = stable present.
#' - **Stable/increase/decrease** (categorical, from `deltaP` and
#'   `changeThresh`): -1 = decrease, 0 = stable, 1 = increase.
#'
#' @param species Character vector of Latin species names.
#' @param yearRef Integer. The earlier ("reference") year.
#' @param yearCurrent Integer. The later ("current") year.
#' @param metaDir Character. Directory holding metaModel()'s output.
#' @param changeThresh Numeric. Minimum absolute `deltaP` to call a cell
#'   "increase"/"decrease" rather than "stable". Default 0.05 -- a
#'   judgment call, state it explicitly in any methods write-up.
#' @return List with:
#'   - `perSpecies`: named list (by species) of 3-layer SpatRasters
#'     (`deltaP`, `gainLoss`, `stableIncDec`), or `NULL` for a species
#'     missing either year.
#'   - `community`: a 4-layer SpatRaster (`meanDeltaP`, `netGainLoss`,
#'     `gainCount`, `lossCount`), aggregated across all species with valid
#'     data for both years.
computeChangeMaps <- function(species, yearRef, yearCurrent, metaDir, changeThresh = 0.05) {
  perSpecies <- list()
  deltaPLayers <- list()
  signedLayers <- list()
  gainLayers <- list()
  lossLayers <- list()

  readBands <- function(sp, yr) {
    spClean <- gsub(" ", "_", sp)
    f <- file.path(metaDir, paste0(spClean, "_meta_suitability_", yr, ".tif"))
    if (!file.exists(f)) return(NULL)
    tryCatch(terra::rast(f), error = function(e) NULL)
  }

  for (sp in species) {
    rRef <- readBands(sp, yearRef)
    rCur <- readBands(sp, yearCurrent)

    if (is.null(rRef) || is.null(rCur) ||
        !all(c("meta_prob", "binary") %in% names(rRef)) ||
        !all(c("meta_prob", "binary") %in% names(rCur))) {
      perSpecies[[sp]] <- NULL
      next
    }

    deltaP <- rCur[["meta_prob"]] - rRef[["meta_prob"]]
    names(deltaP) <- "deltaP"

    gainLoss <- rRef[["binary"]] + 2 * rCur[["binary"]]
    names(gainLoss) <- "gainLoss"

    stableIncDec <- terra::classify(
      deltaP,
      rcl = matrix(c(-Inf, -changeThresh, -1,
                      -changeThresh, changeThresh, 0,
                      changeThresh, Inf, 1),
                    ncol = 3, byrow = TRUE)
    )
    names(stableIncDec) <- "stableIncDec"

    perSpecies[[sp]] <- c(deltaP, gainLoss, stableIncDec)

    deltaPLayers[[sp]] <- deltaP
    signedLayers[[sp]] <- (gainLoss == 2) - (gainLoss == 1)  # +1 gain, -1 loss, 0 else
    gainLayers[[sp]] <- (gainLoss == 2)
    lossLayers[[sp]] <- (gainLoss == 1)
  }

  if (length(deltaPLayers) == 0) {
    warning("No species had valid data for both ", yearRef, " and ", yearCurrent,
            " -- community layer is NULL")
    community <- NULL
  } else {
    meanDeltaP <- mean(terra::rast(deltaPLayers), na.rm = TRUE)
    names(meanDeltaP) <- "meanDeltaP"
    netGainLoss <- sum(terra::rast(signedLayers), na.rm = TRUE)
    names(netGainLoss) <- "netGainLoss"
    gainCount <- sum(terra::rast(gainLayers), na.rm = TRUE)
    names(gainCount) <- "gainCount"
    lossCount <- sum(terra::rast(lossLayers), na.rm = TRUE)
    names(lossCount) <- "lossCount"
    community <- c(meanDeltaP, netGainLoss, gainCount, lossCount)
  }

  list(perSpecies = perSpecies, community = community)
}
