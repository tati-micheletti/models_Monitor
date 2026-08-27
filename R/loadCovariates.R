#' Load the landscape covariate stack for one year
#'
#' Loads year-specific land use, the CORINE land cover snapshot, and the
#' static DEM derivatives, resamples everything to the land use reference
#' grid, and backfills the hedges layer for years where it is unavailable.
#'
#' NOTE: this expects `landuse_<year>_landscape.tif` and
#' `landcover_<corineYear>_landscape.tif` to already exist under
#' `landscapeDir` -- produced by dataPrep_Monitor's `prepareLanduse()`/
#' `prepareLandcover()`.
#'
#' NOTE: this exact function is deliberately duplicated verbatim in
#' models_Monitor (same filename) so that module has no load-time
#' dependency on dataPrep_Monitor being loaded. Keep both copies
#' byte-identical -- tools/check_duplicated_functions.R checks this.
#'
#' @param year Integer. Target year for covariate loading.
#' @param landscapeDir Path to the landscape-scale outputs folder.
#' @param habitatDir Path to the habitat-scale outputs folder (reserved
#'   for future habitat-scale use; unused for now).
#' @return SpatRaster stack with all landscape covariates, or NULL if
#'   required files are missing.
loadCovariates <- function(year, landscapeDir, habitatDir) {

  corineYr <- corineYear(year)
  message("Loading covariates for ", year, " (CORINE: ", corineYr, ")")

  # Land use (year-specific, 14 categories)
  luFile <- file.path(landscapeDir, paste0("landuse_", year, "_landscape.tif"))
  if (!file.exists(luFile)) {
    warning("Land use file missing for ", year, ": ", luFile)
    return(NULL)
  }
  lu <- terra::rast(luFile)

  # Land cover (CORINE snapshot, 3 categories)
  lcFile <- file.path(landscapeDir, paste0("landcover_", corineYr, "_landscape.tif"))
  if (!file.exists(lcFile)) {
    warning("Land cover file missing for CORINE ", corineYr, ": ", lcFile)
    return(NULL)
  }
  lc <- terra::rast(lcFile)

  # DEM derivatives (static)
  elev <- terra::rast(file.path(landscapeDir, "elevation_landscape.tif"))
  slope <- terra::rast(file.path(landscapeDir, "slope_landscape.tif"))
  solar <- terra::rast(file.path(landscapeDir, "solar_radiation_landscape.tif"))

  names(lu) <- paste0(names(lu), "_", year)
  names(lc) <- paste0(names(lc), "_", corineYr)

  # Resample to land use reference grid
  # Resolves ~73m origin offset between land cover/DEM and land use grids
  # (both EPSG:3035 at 1km but different grid origins due to different
  # source data)
  luRef <- lu[[1]]
  lc <- terra::resample(lc, luRef, method = "bilinear")
  elev <- terra::resample(elev, luRef, method = "bilinear")
  slope <- terra::resample(slope, luRef, method = "bilinear")
  solar <- terra::resample(solar, luRef, method = "bilinear")

  # Hedges back/forward-filling
  # Hedge data availability:
  #   2017-2022, 2024-2025: real data
  #   2005-2016: backfill from 2017 (earliest available)
  #   2022-2023: forwardfill from 2021 (nearest preceding year)
  # Always use nearest PRECEDING year where possible to avoid using future
  # landscape data for past occurrences. Exception: 2005-2016 use 2017
  # since no earlier data exists.
  hedgeCol <- names(lu)[grepl("^hedges", names(lu))]

  if (length(hedgeCol) > 0) {
    hedgeVals <- terra::values(lu[[hedgeCol]])

    if (all(is.na(hedgeVals))) {
      refYear <- if (year <= 2016) 2017L else
        if (year %in% c(2022, 2023)) 2021L else
          NULL

      if (!is.null(refYear)) {
        message("Hedges NA for ", year, " -- backfilling from ", refYear)

        luRefYr <- terra::rast(file.path(landscapeDir,
                                          paste0("landuse_", refYear, "_landscape.tif")))
        hedgeRef <- luRefYr[["hedges"]]
        hedgeRef <- terra::resample(hedgeRef, lu[[1]], method = "bilinear")
        names(hedgeRef) <- hedgeCol
        lu[[hedgeCol]] <- hedgeRef

        message("Hedges backfilled from ", refYear, " (mean: ",
                round(mean(terra::values(hedgeRef), na.rm = TRUE), 4), ")")
      } else {
        message("Hedges NA for ", year, " -- no backfill rule defined")
      }
    }
  }

  covStack <- c(lu, lc, elev, slope, solar)
  message("Covariate stack: ", terra::nlyr(covStack), " layers")
  covStack
}
