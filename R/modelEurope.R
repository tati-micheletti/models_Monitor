#' Train the European climate BRT and predict onto all rolling-window
#' climatologies
#'
#' For each species: trains a BRT on `inputsData[[sp]]$data` (inputs_Monitor's
#' final europe-scale table -- already has the resolved predictor columns,
#' `occurrence`, and `foldID`), block-cross-validates using that table's own
#' `foldID` column, then predicts onto every `bioclim_<start>-<end>.tif` in
#' `climateTargetYears` (one prediction raster per target year -- these feed
#' the meta-model as the climate suitability covariate).
#'
#' @param inputsData Named list (by species) with `data`/`predictors`, from
#'   `sim$inputsData$europe` (inputs_Monitor).
#' @param climateTargetYears Integer vector of target years to predict onto.
#' @param climateWindowLength Integer. Rolling window length in years.
#' @param climateOutputDir Character. Directory of bioclim_<start>-<end>.tif files.
#' @param outputDir Character. Directory to save model/performance/prediction outputs in.
#' @param initialLR Numeric. Default starting learning rate for `optimizeBRT()`,
#'   used when a species has neither a `perSpeciesLR` override nor a
#'   previously-persisted converged LR (see `resolveStartingLR()`).
#' @param perSpeciesLR Named numeric vector/list, or NULL (default). Per-species
#'   starting-LR overrides, keyed by species Latin name.
#' @return Named list (by species) with `modelPath`, `perfPath`, `predictions`
#'   (named by year), and `perf` (the evalSDM() row).
modelEurope <- function(inputsData, climateTargetYears, climateWindowLength,
                         climateOutputDir, outputDir, initialLR = 0.01, perSpeciesLR = NULL) {

  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  predictionFiles <- lapply(climateTargetYears, function(yr) {
    startYr <- yr - (climateWindowLength - 1)
    fname <- paste0("bioclim_", startYr, "-", yr, ".tif")
    fpath <- file.path(climateOutputDir, fname)
    list(year = yr, path = fpath, exists = file.exists(fpath))
  })

  nMissing <- sum(!sapply(predictionFiles, function(x) x$exists))
  if (nMissing > 0) {
    missingYrs <- sapply(predictionFiles[!sapply(predictionFiles, function(x) x$exists)],
                          function(x) x$year)
    warning(nMissing, " climatology files missing for years: ", paste(missingYrs, collapse = ", "))
  }

  result <- list()

  for (sp in names(inputsData)) {
    spClean <- gsub(" ", "_", sp)
    message("\n  == ", sp, " ==============================")

    spPa <- inputsData[[sp]]$data
    predSel <- inputsData[[sp]]$predictors

    outModel <- file.path(outputDir, paste0(spClean, "_BRT_EU.rds"))
    outPerf <- file.path(outputDir, paste0(spClean, "_perf_EU.rds"))

    message("Records: ", nrow(spPa), " (", sum(spPa$occurrence == 1), " pres / ",
            sum(spPa$occurrence == 0), " abs)")
    message("Predictors (", length(predSel), "): ", paste(predSel, collapse = ", "))

    if (isValidCachedRDS(outModel)) {
      message("Loading cached BRT model...")
      brtM <- readRDS(outModel)
    } else {
      startingLR <- resolveStartingLR(sp, spClean, defaultLR = initialLR, lrStateDir = outputDir,
                                       perSpeciesLR = perSpeciesLR, lrStateSuffix = "_EU")
      message("Training BRT (optimising learning rate, starting from ", startingLR, ")...")
      brtM <- optimizeBRT(spPa, predSel, "occurrence", startingLR)
      if (is.null(brtM)) {
        warning("Skipping ", sp, " at Europe scale -- optimizeBRT() gave up (see its own ",
                "warning above for why). No model saved; this species will simply be ",
                "absent from Europe-scale results until its data issue is fixed.")
        next
      }
      persistConvergedLR(brtM, spClean, lrStateDir = outputDir, lrStateSuffix = "_EU")
      saveRDS(brtM, outModel)
      message("Model saved -> ", outModel)
    }

    if (isValidCachedRDS(outPerf)) {
      message("Loading cached performance metrics...")
      brtPerf <- readRDS(outPerf)
    } else {
      message("Running block cross-validation...")
      cvPred <- blockCVPredictBRT(spPa, predSel, "occurrence", spPa$foldID, brtM)
      brtPerf <- evalSDM(spPa$occurrence, cvPred)
      saveRDS(brtPerf, outPerf)

      message("Performance: AUC = ", round(brtPerf$AUC, 3), " | TSS = ", round(brtPerf$TSS, 3),
              " | D2 = ", round(brtPerf$D2, 3))
      if (brtPerf$AUC < 0.7) warning("AUC < 0.7 for ", sp, " -- interpret with caution.")

      predOcc <- gbm::predict.gbm(brtM, spPa[, predSel], n.trees = brtM$gbm.call$best.trees,
                                   type = "response")
      exclDev <- explDeviance(spPa$occurrence, predOcc)
      message("D2 at occurrence locations: ", round(exclDev, 3))
      saveRDS(data.frame(species = sp, D2 = exclDev),
              file.path(outputDir, paste0(spClean, "_expl_dev_EU.rds")))
    }

    message("Predicting onto ", length(predictionFiles), " climatologies...")
    spPredFiles <- list()

    for (pf in predictionFiles) {
      if (!pf$exists) {
        message("Year ", pf$year, ": climatology missing -- skipping")
        next
      }

      outTif <- file.path(outputDir, paste0(spClean, "_pred_EU_", pf$year, ".tif"))

      if (isValidPredictionRaster(outTif)) {
        message("Year ", pf$year, ": cache hit")
        spPredFiles[[as.character(pf$year)]] <- outTif
        next
      }

      bioclimYr <- terra::rast(pf$path)
      predictBRTToRaster(bioclimYr, predSel, brtM, brtPerf$thresh, outTif)
      message("Year ", pf$year, ": saved -> ", basename(outTif))
      spPredFiles[[as.character(pf$year)]] <- outTif
    }

    result[[sp]] <- list(modelPath = outModel, perfPath = outPerf,
                          predictions = spPredFiles, perf = brtPerf)
    message("  == Done: ", sp, " ==============================")
  }

  result
}
