#' Load the habitat (200m) covariate stack for one year
#'
#' Habitat-scale equivalent of `loadCovariates()` (which is landscape/1km
#' scale) -- loads year-specific land use, the CORINE land cover
#' snapshot, and the static DEM derivatives at habitat resolution,
#' resamples everything to the land use reference grid, and backfills
#' the hedges layer for years where it is unavailable.
#'
#' NOTE: this expects `landuse_<year>_habitat_<scaleLabel>.tif` and
#' `landcover_<corineYear>_habitat_<scaleLabel>.tif` to already exist under
#' `habitatDir` -- produced by dataPrep_Monitor's `prepareLanduse()`/
#' `prepareLandcover()`. The `<scaleLabel>` suffix (e.g. `scale_02`) is
#' `basename(habitatDir)` itself -- a second, redundant safety layer
#' beyond the containing folder name (2026-09-28).
#'
#' NOTE: this exact function is deliberately duplicated verbatim in
#' models_Monitor (same filename) so that module has no load-time
#' dependency on dataPrep_Monitor being loaded. Keep both copies
#' byte-identical -- tools/check_duplicated_functions.R checks this.
#'
#' @param year Integer. Target year for covariate loading.
#' @param habitatDir Path to the habitat-scale outputs folder.
#' @return SpatRaster stack with all habitat covariates, or NULL if
#'   required files are missing.
loadHabitatCovariates <- function(year, habitatDir) {

  corineYr <- corineYear(year)
  message("Loading habitat covariates for ", year, " (CORINE: ", corineYr, ")")

  # Resolution appended to every filename below (second safety layer beyond
  # the containing scaleLabel()-named folder, `habitatDir` itself) -- see
  # aggregateAndSave.R/computeLanduse.R/computeLandcover.R (the writers).
  resSuffix <- basename(habitatDir)

  luFile <- file.path(habitatDir, paste0("landuse_", year, "_habitat_", resSuffix, ".tif"))
  lcFile <- file.path(habitatDir, paste0("landcover_", corineYr, "_habitat_", resSuffix, ".tif"))
  # solar_radiation intentionally excluded (see DECISIONS.md, 2026-09-26 --
  # dropped as a predictor for every species, not well-scaled/possibly
  # capturing noise from other unmodeled factors)
  demFiles <- c(file.path(habitatDir, paste0("elevation_habitat_", resSuffix, ".tif")),
                file.path(habitatDir, paste0("slope_habitat_", resSuffix, ".tif")))

  missing <- c(luFile, lcFile, demFiles)[!file.exists(c(luFile, lcFile, demFiles))]
  if (length(missing) > 0) {
    warning("Missing habitat covariate files: ", paste(missing, collapse = ", "))
    return(NULL)
  }

  lu <- terra::rast(luFile)
  lc <- terra::rast(lcFile)
  elev <- terra::rast(demFiles[1])
  slope <- terra::rast(demFiles[2])

  # Resample all to land use reference grid (resolves origin offsets
  # between layers, same fix as landscape scale in loadCovariates())
  luRef <- lu[[1]]
  lc <- terra::resample(lc, luRef, method = "bilinear")
  elev <- terra::resample(elev, luRef, method = "bilinear")
  slope <- terra::resample(slope, luRef, method = "bilinear")

  covStack <- c(lu, lc, elev, slope)

  # Strip year/CORINE suffixes so column names match the training data
  names(covStack) <- gsub("_\\d{4}$", "", names(covStack))

  # Hedges backfilling (200m scale) -- same rule as loadCovariates():
  # available 2017-2021, 2024-2025; missing 2005-2016 backfilled from
  # 2017; missing 2022-2023 backfilled from 2021.
  hedgeCol <- names(covStack)[grepl("^hedges", names(covStack))]

  if (length(hedgeCol) > 0) {
    hedgeVals <- terra::values(covStack[[hedgeCol]])
    if (all(is.na(hedgeVals))) {
      refYear <- if (year <= 2016) 2017L else
        if (year %in% c(2022, 2023)) 2021L else
          NULL

      if (!is.null(refYear)) {
        message("Hedges NA for ", year, " -- backfilling from ", refYear)
        luRefYr <- terra::rast(file.path(habitatDir,
                                          paste0("landuse_", refYear, "_habitat_", resSuffix, ".tif")))
        hedgeRef <- luRefYr[["hedges"]]
        hedgeRef <- terra::resample(hedgeRef, covStack[[1]], method = "bilinear")
        names(hedgeRef) <- hedgeCol
        covStack[[hedgeCol]] <- hedgeRef
        message("Hedges backfilled from ", refYear)
      } else {
        message("Hedges NA for ", year, " -- no backfill rule defined")
      }
    }
  }

  message("Covariate stack: ", terra::nlyr(covStack), " layers")
  covStack
}
