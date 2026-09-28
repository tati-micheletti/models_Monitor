#' Train the German landscape BRT (1km) and predict onto the full German grid
#'
#' For each species: trains ONE BRT (pooled across years) on
#' `inputsData[[sp]]$data` (inputs_Monitor's final landscape-scale table),
#' block-cross-validates using that table's own `foldID` column, then
#' predicts onto `predictionYears` (the full German modeling period).
#'
#' Each species' own landscape resolution is resolved from `resolutionConfig`
#' (falling back to `sharedResolutionM` when unset), so covariates and
#' outputs are read/written to that species' own `scale_X` folder -- species
#' sharing a resolution also share one cached covariate stack per year (e.g.
#' Milvus milvus at a coarser 15km landscape window than the shared 1km
#' default -- see DECISIONS.md's 2026-09-28 entry).
#'
#' @param inputsData Named list (by species) with `data`/`predictors`, from
#'   `sim$inputsData$gerLandscape` (inputs_Monitor).
#' @param predictionYears Integer vector of years to predict onto.
#' @param processedRoot Character. `predictors/processed` directory (without
#'   the scale_X leaf -- each species' own leaf is appended internally).
#' @param outputRoot Character. Directory to save model/performance/
#'   prediction outputs in (without the scale_X leaf).
#' @param resolutionConfig Named list, species -> scale -> resolution (m), or
#'   NULL. See `models_Monitor`'s parameter of the same name.
#' @param sharedResolutionM Numeric. Shared default landscape resolution (m),
#'   used for any species absent from `resolutionConfig`.
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
modelGerLandscape <- function(inputsData, predictionYears, processedRoot, outputRoot,
                               resolutionConfig = NULL, sharedResolutionM,
                               initialLR = 0.08, perSpeciesLR = NULL, cachePath = NULL) {

  if (is.null(cachePath)) cachePath <- file.path(tempdir(), "birdMonitor_cache")
  result <- list()
  covStacksByResolution <- list()

  for (sp in names(inputsData)) {
    spClean <- gsub(" ", "_", sp)
    message("\n  == ", sp, " ==============================")

    resM <- resolveResolutionM(sp, "landscape", resolutionConfig, sharedResolutionM)
    landscapeOutputDir <- file.path(processedRoot, scaleLabel(resM))
    outputDir <- file.path(outputRoot, scaleLabel(resM))
    dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

    resKey <- as.character(resM)
    if (is.null(covStacksByResolution[[resKey]])) {
      covStacksByResolution[[resKey]] <- stats::setNames(
        lapply(predictionYears, function(yr) {
          covStack <- loadCovariates(yr, landscapeOutputDir, NULL)
          if (!is.null(covStack)) names(covStack) <- gsub("_\\d{4}$", "", names(covStack))
          covStack
        }),
        as.character(predictionYears))
    }
    covStacksByYear <- covStacksByResolution[[resKey]]

    spPa <- inputsData[[sp]]$data
    predSel <- inputsData[[sp]]$predictors

    outModel <- file.path(outputDir, paste0(spClean, "_BRT_landscape.rds"))
    outPerf <- file.path(outputDir, paste0(spClean, "_perf_landscape.rds"))

    message("Records: ", nrow(spPa), " (", sum(spPa$occurrence == 1), " pres / ",
            sum(spPa$occurrence == 0), " abs)")
    message("Predictors (", length(predSel), "): ", paste(predSel, collapse = ", "))

    startingLR <- resolveStartingLR(sp, spClean, defaultLR = initialLR, lrStateDir = outputDir,
                                     perSpeciesLR = perSpeciesLR, lrStateSuffix = "_landscape")
    brtM <- reproducible::Cache(
      fitBRTOneSpecies, sp = sp, spPa = spPa, predSel = predSel, startingLR = startingLR,
      cachePath = cachePath, userTags = c("modelGerLandscape", "fit", spClean))

    if (is.null(brtM)) {
      warning("Skipping ", sp, " at landscape scale -- optimizeBRT() gave up (see its own ",
              "warning above for why). No model saved; this species will simply be ",
              "absent from landscape-scale results until its data issue is fixed.")
      next
    }
    persistConvergedLR(brtM, spClean, lrStateDir = outputDir, lrStateSuffix = "_landscape")
    saveRDS(brtM, outModel)
    message("Model saved -> ", outModel)

    evalResult <- reproducible::Cache(
      evalBRTOneSpecies, sp = sp, spPa = spPa, predSel = predSel, brtM = brtM,
      cachePath = cachePath, userTags = c("modelGerLandscape", "eval", spClean))
    brtPerf <- evalResult$perf
    saveRDS(brtPerf, outPerf)

    message("Performance: AUC = ", round(brtPerf$AUC, 3), " | TSS = ", round(brtPerf$TSS, 3),
            " | D2 = ", round(brtPerf$D2, 3))
    if (brtPerf$AUC < 0.7) warning("AUC < 0.7 for ", sp, " -- interpret with caution.")
    message("D2 at occurrence locations: ", round(evalResult$explDev, 3))
    saveRDS(data.frame(species = sp, D2 = evalResult$explDev),
            file.path(outputDir, paste0(spClean, "_expl_dev_landscape.rds")))

    message("Predicting onto German grid for ", length(predictionYears), " years...")
    spPredFiles <- list()

    for (yr in predictionYears) {
      outTif <- file.path(outputDir, paste0(spClean, "_pred_landscape_", yr, ".tif"))

      covStack <- covStacksByYear[[as.character(yr)]]
      if (is.null(covStack)) {
        warning("Year ", yr, ": covariates unavailable -- skipping")
        next
      }

      # x/y are never actual covStack layers -- predictBRTToRaster() derives
      # them from the raster's own cell coordinates (see DECISIONS.md,
      # 2026-09-26), so they're never "missing" here.
      missingPreds <- setdiff(predSel, c(names(covStack), "x", "y"))
      if (length(missingPreds) > 0) {
        warning("Year ", yr, ": missing predictors: ", paste(missingPreds, collapse = ", "))
        next
      }

      predRaster <- reproducible::Cache(
        predictBRTToRaster, covStack = covStack, predictors = predSel, brtModel = brtM,
        thresh = brtPerf$thresh, cachePath = cachePath,
        userTags = c("modelGerLandscape", "predict", spClean, as.character(yr)))
      terra::writeRaster(predRaster, outTif, overwrite = TRUE)
      message("Year ", yr, ": saved -> ", basename(outTif))
      spPredFiles[[as.character(yr)]] <- outTif
    }

    result[[sp]] <- list(modelPath = outModel, perfPath = outPerf,
                          predictions = spPredFiles, perf = brtPerf)
    message("  == Done: ", sp, " ==============================")
  }

  result
}

#' Fit one species' BRT (the actual computation `modelGerLandscape()`,
#' `modelGerHabitat()` and `modelEurope()` each wrap in `reproducible::Cache()`)
#'
#' Pulled into its own function so Cache()'s digest covers exactly
#' `spPa`/`predSel`/`startingLR` -- i.e. it changes (and only that
#' species refits) when inputs_Monitor's output for this species changes,
#' or its BRT starting learning rate changes.
#'
#' @param sp Character. Species Latin name (for the training message only).
#' @param spPa data.frame. This species' model-ready table (inputs_Monitor's output).
#' @param predSel Character vector. Resolved predictor columns.
#' @param startingLR Numeric. Resolved starting learning rate.
#' @return A fitted `gbm.step()` model, or NULL if `optimizeBRT()` gave up.
fitBRTOneSpecies <- function(sp, spPa, predSel, startingLR) {
  message("Training BRT (optimising learning rate, starting from ", startingLR, ")...")
  optimizeBRT(spPa, predSel, "occurrence", startingLR)
}

#' Evaluate one species' fitted BRT via block cross-validation
#'
#' Same Cache()-ing rationale as `fitBRTOneSpecies()` -- digest covers
#' `spPa`/`predSel`/`brtM`, so a refit (different `brtM`) or a changed
#' input table naturally triggers re-evaluation too.
#'
#' @param sp Character. Species Latin name (unused internally; kept for a
#'   stable, self-documenting Cache() digest).
#' @param spPa data.frame. This species' model-ready table.
#' @param predSel Character vector. Resolved predictor columns.
#' @param brtM A fitted `gbm.step()` model.
#' @return List with `perf` (the `evalSDM()` row) and `explDev` (deviance
#'   explained at occurrence locations).
evalBRTOneSpecies <- function(sp, spPa, predSel, brtM) {
  message("Running block cross-validation...")
  cvPred <- blockCVPredictBRT(spPa, predSel, "occurrence", spPa$foldID, brtM)
  brtPerf <- evalSDM(spPa$occurrence, cvPred)

  predOcc <- gbm::predict.gbm(brtM, spPa[, predSel], n.trees = brtM$gbm.call$best.trees,
                               type = "response")
  exclDev <- explDeviance(spPa$occurrence, predOcc)

  list(perf = brtPerf, explDev = exclDev)
}
