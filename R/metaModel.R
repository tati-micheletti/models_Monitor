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
#' @param habitatYears Integer vector of years with real habitat occurrence
#'   data (the meta-model's training years).
#' @param predictionYears Integer vector of years to predict onto.
#' @param modelDirs Named list with `europe`, `landscape`, `habitat` prediction directories.
#' @param refRaster SpatRaster. Shared 200m reference grid (a static
#'   DEM-derived raster from dataPrep_Monitor -- deliberately NOT any
#'   species' habitat prediction, so the reference grid never depends on
#'   model output ordering).
#' @param outputDir Character. Directory to save model/performance/prediction outputs in.
#' @return Named list (by species) with `modelPath`, `perfPath`, `varimpPath`,
#'   `predictions` (named by year), `perf`, and `varimp`.
metaModel <- function(inputsDataGerHabitat, habitatYears, predictionYears, modelDirs,
                       refRaster, outputDir) {

  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  suitCols <- c("climate_mean_prob", "landscape_mean_prob", "habitat_mean_prob")
  idCols <- c("AREA_NATIONAL_CODE", "latin_name", "occurrence", "x", "y", "foldID")

  result <- list()

  for (sp in names(inputsDataGerHabitat)) {
    spClean <- gsub(" ", "_", sp)
    message("\n  == ", sp, " ==============================")

    spPaAll <- inputsDataGerHabitat[[sp]]$data

    outModel <- file.path(outputDir, paste0(spClean, "_ridge_meta.rds"))
    outPerf <- file.path(outputDir, paste0(spClean, "_perf_meta.rds"))
    outVarimp <- file.path(outputDir, paste0(spClean, "_varimp_meta.rds"))

    message("Extracting suitability scores (training years ",
            paste(range(habitatYears), collapse = "-"), ")...")

    trainList <- lapply(habitatYears, function(yr) {
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

    if (isValidCachedRDS(outModel)) {
      message("Loading cached ridge model...")
      ridgeM <- readRDS(outModel)
    } else {
      message("Training ridge regression (alpha = 0)...")
      set.seed(42)
      ridgeCv <- glmnet::cv.glmnet(x = X, y = y, family = "binomial", alpha = 0,
                                    nfolds = 10, standardize = TRUE)
      ridgeM <- list(model = ridgeCv, lambda = ridgeCv$lambda.1se, suit_cols = suitCols)

      coefs <- stats::coef(ridgeCv, s = ridgeCv$lambda.1se)
      message("Coefficients: intercept = ", round(coefs[1], 4),
              " | climate = ", round(coefs[2], 4),
              " | landscape = ", round(coefs[3], 4),
              " | habitat = ", round(coefs[4], 4))

      saveRDS(ridgeM, outModel)
      message("Model saved -> ", outModel)
    }

    if (isValidCachedRDS(outVarimp)) {
      message("Loading cached variable importance...")
      varimp <- readRDS(outVarimp)
    } else {
      message("Computing variable importance (leave-one-out)...")
      varimp <- computeVariableImportance(trainDf, suitCols, ridgeM$lambda)
      saveRDS(varimp, outVarimp)
      message(sprintf("Importance: climate = %.3f | landscape = %.3f | habitat = %.3f",
                       varimp$imp_climate, varimp$imp_landscape, varimp$imp_habitat))
    }

    if (isValidCachedRDS(outPerf)) {
      message("Loading cached performance...")
      metaPerf <- readRDS(outPerf)
    } else {
      message("Running block cross-validation...")
      cvPred <- blockCVPredictRidge(X, y, trainDf$foldID)
      validIdx <- !is.na(cvPred)
      metaPerf <- evalSDM(y[validIdx], cvPred[validIdx])
      saveRDS(metaPerf, outPerf)

      message("Performance: AUC = ", round(metaPerf$AUC, 3), " | TSS = ", round(metaPerf$TSS, 3),
              " | D2 = ", round(metaPerf$D2, 3))
      if (metaPerf$AUC < 0.7) warning("AUC < 0.7 for ", sp, " -- interpret with caution.")
    }

    message("Predicting onto German grid for ", length(predictionYears), " years...")
    message("(years before ", min(habitatYears), " are hindcasts -- trained on ",
            paste(range(habitatYears), collapse = "-"), " relationships)")

    spPredFiles <- list()

    for (yr in predictionYears) {
      outTif <- file.path(outputDir, paste0(spClean, "_meta_suitability_", yr, ".tif"))

      if (isValidPredictionRaster(outTif, layerName = "meta_prob")) {
        message("Year ", yr, ": cache hit")
        spPredFiles[[as.character(yr)]] <- outTif
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

      predictRidgeToRaster(suitStack, suitCols, ridgeM$model, ridgeM$lambda, metaPerf$thresh, outTif)
      message("Year ", yr, ": saved -> ", basename(outTif))
      spPredFiles[[as.character(yr)]] <- outTif
    }

    result[[sp]] <- list(modelPath = outModel, perfPath = outPerf, varimpPath = outVarimp,
                          predictions = spPredFiles, perf = metaPerf, varimp = varimp)
    message("  == Done: ", sp, " ==============================")
  }

  result
}
