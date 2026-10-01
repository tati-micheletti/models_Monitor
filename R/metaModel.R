#' Train the ridge regression meta-model and predict to all years
#'
#' Combines the three scale-specific BRT suitability predictions
#' (climate/landscape/habitat, all resampled to a shared 200m grid) via
#' ridge regression (`glmnet`, alpha = 0). Trains on pooled habitat
#' occurrence (`habitatYears`, 2022-2025), block-cross-validates using the
#' `foldID` column already attached to `inputsDataGerHabitat` (sourced
#' from inputs_Monitor -- no spatial-join reconstruction needed, unlike
#' the original script, and this makes it robust to spatialBlocking
#' having been real or mimicked upstream), then predicts onto
#' `predictionYears` (the full German modeling period -- 2005-2021 are
#' hindcasts, assuming stable scale weighting over time).
#'
#' @param inputsDataGerHabitat Named list (by species) with `data`, from
#'   `sim$inputsData$gerHabitat` (inputs_Monitor). Must include a `year` column.
#' @param habitatYears Integer vector, or named list (species -> integer
#'   vector), of years with real habitat occurrence data (the meta-model's
#'   training years). A named list lets each species train on its own
#'   real-data window (e.g. Buteo buteo/Sturnus vulgaris's real MhB
#'   point-count data is negligible before ~2020, while other species
#'   genuinely span a wider range) -- see `resolveYearsPerSpecies()` in
#'   `sharedSpeciesConfig.R`. A flat vector applies the same years to every
#'   species (backward compatible).
#' @param predictionYears Integer vector of years to predict onto.
#' @param modelDirs Named list with `europe`, `landscape`, `habitat` prediction directories.
#' @param refRaster SpatRaster. Shared 200m reference grid (a static
#'   DEM-derived raster from dataPrep_Monitor -- deliberately NOT any
#'   species' habitat prediction, so the reference grid never depends on
#'   model output ordering).
#' @param outputDir Character. Directory to save model/performance/prediction outputs in.
#' @param nBootTrend Integer. If > 0, also bootstrap the ridge fit
#'   `nBootTrend` times (see `bootstrapMetaModelTrend()`) to quantify
#'   model-fitting uncertainty in the per-year area-mean trend -- distinct
#'   from, and in addition to, spatial-averaging precision. Default 0 (off):
#'   existing runs are unaffected unless this is explicitly requested, since
#'   it adds `nBootTrend` extra ridge refits per species.
#' @param gadmCacheDir Character. Directory to cache the GADM Germany
#'   boundary in (only used, and only fetched, when `nBootTrend > 0`) --
#'   see `bootRefRaster` below.
#' @param cachePath Character, or NULL (default). Directory for
#'   `reproducible::Cache()`'s per-species(-year) cache -- e.g.
#'   `cachePath(sim)`, a stable location shared across runs (NOT the
#'   per-run timestamped `outputDir`), so an unchanged species/config
#'   reuses its fit across separate `runMe.R` invocations instead of
#'   refitting from scratch every run -- matching `modelGerHabitat()`/
#'   `modelGerLandscape()`/`modelEurope()`'s own convention. NULL falls
#'   back to a temp directory, for standalone/test calls.
#' @return Named list (by species) with `modelPath`, `perfPath`, `varimpPath`,
#'   `predictions` (named by year), `perf`, `varimp`, and (if `nBootTrend > 0`)
#'   `trendBootPath`/`trendBoot`.
metaModel <- function(inputsDataGerHabitat, habitatYears, predictionYears, modelDirs,
                       refRaster, outputDir, nBootTrend = 0, gadmCacheDir = NULL,
                       cachePath = NULL) {

  if (is.null(cachePath)) cachePath <- file.path(tempdir(), "birdMonitor_cache")
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  suitCols <- c("climate_mean_prob", "landscape_mean_prob", "habitat_mean_prob")
  idCols <- c("AREA_NATIONAL_CODE", "latin_name", "occurrence", "x", "y", "foldID")

  # A Germany-cropped copy of refRaster, used ONLY for getOrBuildSuitX()'s
  # resample() calls below -- refRaster itself (and every prediction this
  # function writes) stays at its original, uncropped, full-Europe extent,
  # unchanged from before nBootTrend existed. Cropping first is what makes
  # that resample affordable: it scales with the target grid's size, and a
  # real timed dry run showed ~85 minutes for one species at the
  # uncropped extent (~738M cells, ~1.9% real German data) -- Germany's
  # real extent is a small fraction of that.
  bootRefRaster <- refRaster
  if (nBootTrend > 0) {
    germany <- geodata::gadm(country = "DEU", level = 0,
                              path = if (is.null(gadmCacheDir)) tempdir() else gadmCacheDir)
    germanyProj <- terra::project(germany, terra::crs(refRaster))
    bootRefRaster <- terra::crop(refRaster, germanyProj)
  }

  result <- list()

  for (sp in names(inputsDataGerHabitat)) {
    spClean <- gsub(" ", "_", sp)
    message("\n  == ", sp, " ==============================")

    spPaAll <- inputsDataGerHabitat[[sp]]$data

    outModel <- file.path(outputDir, paste0(spClean, "_ridge_meta.rds"))
    outPerf <- file.path(outputDir, paste0(spClean, "_perf_meta.rds"))
    outVarimp <- file.path(outputDir, paste0(spClean, "_varimp_meta.rds"))

    habitatYearsForSp <- if (is.list(habitatYears)) habitatYears[[sp]] else habitatYears

    message("Extracting suitability scores (training years ",
            paste(range(habitatYearsForSp), collapse = "-"), ")...")

    trainList <- lapply(habitatYearsForSp, function(yr) {
      spYr <- spPaAll[spPaAll$year == yr, ]
      if (nrow(spYr) == 0) return(NULL)
      message("Year ", yr, ": ", nrow(spYr), " records")
      extractSuitability(spYr, spClean, yr, modelDirs, refRaster, idCols)
    })
    trainList <- trainList[!sapply(trainList, is.null)]

    if (length(trainList) == 0) {
      warning("No training data for ", sp, " -- skipping")
      next
    }

    trainDf <- do.call(rbind, trainList)
    trainDf <- trainDf[stats::complete.cases(trainDf[, suitCols]), ]

    nPres <- sum(trainDf$occurrence == 1)
    nAbs <- sum(trainDf$occurrence == 0)
    message("Training data: ", nrow(trainDf), " (", nPres, " pres / ", nAbs, " abs)")

    if (nPres < 10 || nAbs < 10) {
      warning("Too few records for ", sp, " -- skipping")
      next
    }

    X <- as.matrix(trainDf[, suitCols])
    y <- trainDf$occurrence

    ridgeM <- reproducible::Cache(
      fitRidgeOneSpecies, X = X, y = y, suitCols = suitCols,
      cachePath = cachePath, userTags = c("metaModel", "fit", spClean))
    coefs <- stats::coef(ridgeM$model, s = ridgeM$lambda)
    message("Coefficients: intercept = ", round(coefs[1], 4),
            " | climate = ", round(coefs[2], 4),
            " | landscape = ", round(coefs[3], 4),
            " | habitat = ", round(coefs[4], 4))
    saveRDS(ridgeM, outModel)
    message("Model saved -> ", outModel)

    varimp <- reproducible::Cache(
      computeVariableImportance, trainDf = trainDf, suitCols = suitCols, lambda = ridgeM$lambda,
      cachePath = cachePath, userTags = c("metaModel", "varimp", spClean))
    saveRDS(varimp, outVarimp)
    message(sprintf("Importance: climate = %.3f | landscape = %.3f | habitat = %.3f",
                     varimp$imp_climate, varimp$imp_landscape, varimp$imp_habitat))

    metaPerf <- reproducible::Cache(
      evalRidgeOneSpecies, X = X, y = y, foldID = trainDf$foldID,
      cachePath = cachePath, userTags = c("metaModel", "eval", spClean))
    saveRDS(metaPerf, outPerf)
    message("Performance: AUC = ", round(metaPerf$AUC, 3), " | TSS = ", round(metaPerf$TSS, 3),
            " | D2 = ", round(metaPerf$D2, 3))
    if (metaPerf$AUC < 0.7) warning("AUC < 0.7 for ", sp, " -- interpret with caution.")

    message("Predicting onto German grid for ", length(predictionYears), " years...")
    message("(years before ", min(habitatYearsForSp), " are hindcasts -- trained on ",
            paste(range(habitatYearsForSp), collapse = "-"), " relationships)")

    spPredFiles <- list()
    newXByYear <- list()

    for (yr in predictionYears) {
      outTif <- file.path(outputDir, paste0(spClean, "_meta_suitability_", yr, ".tif"))

      if (isValidPredictionRaster(outTif, layerName = "meta_prob")) {
        message("Year ", yr, ": cache hit")
        spPredFiles[[as.character(yr)]] <- outTif
        if (nBootTrend > 0) {
          suitX <- getOrBuildSuitX(spClean, yr, modelDirs, bootRefRaster, suitCols, outputDir)
          if (!is.null(suitX)) newXByYear[[as.character(yr)]] <- suitX
        }
        next
      }

      rClim <- loadSuitability("climate", spClean, yr, modelDirs, refRaster)
      rLand <- loadSuitability("landscape", spClean, yr, modelDirs, refRaster)
      rHab <- loadSuitability("habitat", spClean, yr, modelDirs, refRaster)

      if (is.null(rClim) || is.null(rLand) || is.null(rHab)) {
        message("Year ", yr, ": missing suitability -- skipping")
        next
      }

      suitStack <- c(rClim, rLand, rHab)
      names(suitStack) <- suitCols

      predRaster <- reproducible::Cache(
        predictRidgeToRaster, suitStack = suitStack, suitCols = suitCols,
        ridgeModel = ridgeM$model, lambda = ridgeM$lambda, thresh = metaPerf$thresh,
        cachePath = cachePath, userTags = c("metaModel", "predict", spClean, as.character(yr)))
      terra::writeRaster(predRaster, outTif, overwrite = TRUE)
      message("Year ", yr, ": saved -> ", basename(outTif))
      spPredFiles[[as.character(yr)]] <- outTif

      if (nBootTrend > 0) {
        suitXPath <- file.path(outputDir, paste0(spClean, "_meta_suitX_", yr, ".rds"))
        predDf <- as.data.frame(suitStack, xy = FALSE, na.rm = TRUE)
        suitX <- as.matrix(predDf[, suitCols])
        saveRDS(suitX, suitXPath)
        newXByYear[[as.character(yr)]] <- suitX
      }
    }

    trendBootPath <- NULL
    trendBoot <- NULL
    if (nBootTrend > 0 && length(newXByYear) > 0) {
      trendBootPath <- file.path(outputDir, paste0(spClean, "_meta_trend_boot.rds"))
      if (isValidCachedRDS(trendBootPath)) {
        message("Loading cached trend bootstrap...")
        trendBoot <- readRDS(trendBootPath)
      } else {
        message("Bootstrapping ridge fit (", nBootTrend, " reps) for trend uncertainty...")
        trendBoot <- bootstrapMetaModelTrend(X, y, newXByYear, nBoot = nBootTrend)
        saveRDS(trendBoot, trendBootPath)
        message("Trend bootstrap saved -> ", trendBootPath)
      }
    }

    result[[sp]] <- list(modelPath = outModel, perfPath = outPerf, varimpPath = outVarimp,
                          predictions = spPredFiles, perf = metaPerf, varimp = varimp,
                          trendBootPath = trendBootPath, trendBoot = trendBoot)
    message("  == Done: ", sp, " ==============================")
  }

  result
}

#' Fit one species' ridge meta-model (the actual computation `metaModel()`
#' wraps in `reproducible::Cache()`)
#'
#' Pulled into its own function so Cache()'s digest covers exactly `X`/`y`
#' -- i.e. it changes (and only that species refits) when the underlying
#' scale-level suitability extraction actually changes, not on every new
#' `runMe.R` run/output folder. Mirrors `fitBRTOneSpecies()`
#' (`modelGerLandscape.R`) -- same rationale, same reason it's Cache()'d
#' by the caller rather than internally.
#'
#' @param X Numeric matrix, the 3 scales' suitability columns.
#' @param y Numeric vector, 0/1 occurrence.
#' @param suitCols Character vector of the 3 suitability column names
#'   (stored alongside the model for `predictRidgeToRaster()`'s own use).
#' @return List with `model` (a `cv.glmnet()` fit), `lambda`
#'   (`lambda.1se`), and `suit_cols`.
fitRidgeOneSpecies <- function(X, y, suitCols) {
  message("Training ridge regression (alpha = 0)...")
  set.seed(42)
  ridgeCv <- glmnet::cv.glmnet(x = X, y = y, family = "binomial", alpha = 0,
                                nfolds = 10, standardize = TRUE)
  list(model = ridgeCv, lambda = ridgeCv$lambda.1se, suit_cols = suitCols)
}

#' Evaluate one species' fitted ridge meta-model via block cross-validation
#'
#' Same Cache()-ing rationale as `fitRidgeOneSpecies()` -- digest covers
#' `X`/`y`/`foldID`, so a changed input table naturally triggers
#' re-evaluation too. Mirrors `evalBRTOneSpecies()` (`modelGerLandscape.R`).
#'
#' @param X Numeric matrix, the 3 scales' suitability columns.
#' @param y Numeric vector, 0/1 occurrence.
#' @param foldID Integer vector, spatial block-CV fold assignment.
#' @return The `evalSDM()` row (AUC/TSS/Kappa/Sens/Spec/PCC/D2/thresh).
evalRidgeOneSpecies <- function(X, y, foldID) {
  message("Running block cross-validation...")
  cvPred <- blockCVPredictRidge(X, y, foldID)
  validIdx <- !is.na(cvPred)
  evalSDM(y[validIdx], cvPred[validIdx])
}
