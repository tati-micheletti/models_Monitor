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
#' @param cachePath Character, or NULL (default). Directory for
#'   `reproducible::Cache()`'s per-species(-year) cache -- e.g.
#'   `cachePath(sim)`, a stable location shared across runs (NOT the
#'   per-run timestamped output folder). NULL falls back to a temp
#'   directory, for standalone/test calls.
#' @return Named list (by species) with `modelPath`, `perfPath`, `predictions`
#'   (named by year), and `perf` (the evalSDM() row).
modelEurope <- function(inputsData, climateTargetYears, climateWindowLength,
                         climateOutputDir, outputDir, initialLR = 0.01, perSpeciesLR = NULL,
                         cachePath = NULL) {

  if (is.null(cachePath)) cachePath <- file.path(tempdir(), "birdMonitor_cache")

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

    startingLR <- resolveStartingLR(sp, spClean, defaultLR = initialLR, lrStateDir = outputDir,
                                     perSpeciesLR = perSpeciesLR, lrStateSuffix = "_EU")
    brtM <- reproducible::Cache(
      fitBRTOneSpecies, sp = sp, spPa = spPa, predSel = predSel, startingLR = startingLR,
      cachePath = cachePath, userTags = c("modelEurope", "fit", spClean))

    if (is.null(brtM)) {
      warning("Skipping ", sp, " at Europe scale -- optimizeBRT() gave up (see its own ",
              "warning above for why). No model saved; this species will simply be ",
              "absent from Europe-scale results until its data issue is fixed.")
      next
    }
    persistConvergedLR(brtM, spClean, lrStateDir = outputDir, lrStateSuffix = "_EU")
    saveRDS(brtM, outModel)
    message("Model saved -> ", outModel)

    evalResult <- reproducible::Cache(
      evalBRTOneSpecies, sp = sp, spPa = spPa, predSel = predSel, brtM = brtM,
      cachePath = cachePath, userTags = c("modelEurope", "eval", spClean))
    brtPerf <- evalResult$perf
    saveRDS(brtPerf, outPerf)

    message("Performance: AUC = ", round(brtPerf$AUC, 3), " | TSS = ", round(brtPerf$TSS, 3),
            " | D2 = ", round(brtPerf$D2, 3))
    if (brtPerf$AUC < 0.7) warning("AUC < 0.7 for ", sp, " -- interpret with caution.")
    message("D2 at occurrence locations: ", round(evalResult$explDev, 3))
    saveRDS(data.frame(species = sp, D2 = evalResult$explDev),
            file.path(outputDir, paste0(spClean, "_expl_dev_EU.rds")))

    message("Predicting onto ", length(predictionFiles), " climatologies...")
    spPredFiles <- list()

    for (pf in predictionFiles) {
      if (!pf$exists) {
        message("Year ", pf$year, ": climatology missing -- skipping")
        next
      }

      outTif <- file.path(outputDir, paste0(spClean, "_pred_EU_", pf$year, ".tif"))

      bioclimYr <- terra::rast(pf$path)
      predRaster <- reproducible::Cache(
        predictBRTToRaster, covStack = bioclimYr, predictors = predSel, brtModel = brtM,
        thresh = brtPerf$thresh, cachePath = cachePath,
        userTags = c("modelEurope", "predict", spClean, as.character(pf$year)))
      terra::writeRaster(predRaster, outTif, overwrite = TRUE)
      message("Year ", pf$year, ": saved -> ", basename(outTif))
      spPredFiles[[as.character(pf$year)]] <- outTif
    }

    result[[sp]] <- list(modelPath = outModel, perfPath = outPerf,
                          predictions = spPredFiles, perf = brtPerf)
    message("  == Done: ", sp, " ==============================")
  }

  result
}
