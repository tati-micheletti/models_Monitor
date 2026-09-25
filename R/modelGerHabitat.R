#' Train the German habitat BRT (200m) and predict onto the full German grid
#'
#' For each species: trains ONE BRT (pooled across the 2022-2025 years that
#' habitat occurrence data covers) on `inputsData[[sp]]$data` (inputs_Monitor's
#' final habitat-scale table), block-cross-validates using that table's own
#' `foldID` column, then predicts onto `predictionYears` (by default the full
#' German modeling period, 2005-2025 -- extrapolating beyond the training
#' years' covariate relationships, matching Wiedenroth et al.).
#'
#' @param inputsData Named list (by species) with `data`/`predictors`, from
#'   `sim$inputsData$gerHabitat` (inputs_Monitor).
#' @param predictionYears Integer vector of years to predict onto.
#' @param habitatOutputDir Character. Directory of habitat-scale covariate rasters.
#' @param outputDir Character. Directory to save model/performance/prediction outputs in.
#' @param initialLR Numeric. Default starting learning rate for `optimizeBRT()`,
#'   used when a species has neither a `perSpeciesLR` override nor a
#'   previously-persisted converged LR (see `resolveStartingLR()`).
#' @param perSpeciesLR Named numeric vector/list, or NULL (default). Per-species
#'   starting-LR overrides, keyed by species Latin name.
#' @return Named list (by species) with `modelPath`, `perfPath`, `predictions`
#'   (named by year), and `perf` (the evalSDM() row).
modelGerHabitat <- function(inputsData, predictionYears, habitatOutputDir, outputDir,
                             initialLR = 0.08, perSpeciesLR = NULL) {

  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  result <- list()

  for (sp in names(inputsData)) {
    spClean <- gsub(" ", "_", sp)
    message("\n  == ", sp, " ==============================")

    spPa <- inputsData[[sp]]$data
    predSel <- inputsData[[sp]]$predictors

    outModel <- file.path(outputDir, paste0(spClean, "_BRT_habitat.rds"))
    outPerf <- file.path(outputDir, paste0(spClean, "_perf_habitat.rds"))

    message("Records: ", nrow(spPa), " (", sum(spPa$occurrence == 1), " pres / ",
            sum(spPa$occurrence == 0), " abs)")
    message("Predictors (", length(predSel), "): ", paste(predSel, collapse = ", "))

    if (isValidCachedRDS(outModel)) {
      message("Loading cached BRT model...")
      brtM <- readRDS(outModel)
    } else {
      startingLR <- resolveStartingLR(sp, spClean, defaultLR = initialLR, lrStateDir = outputDir,
                                       perSpeciesLR = perSpeciesLR, lrStateSuffix = "_habitat")
      message("Training BRT (optimising learning rate, starting from ", startingLR, ")...")
      brtM <- optimizeBRT(spPa, predSel, "occurrence", startingLR)
      persistConvergedLR(brtM, spClean, lrStateDir = outputDir, lrStateSuffix = "_habitat")
      saveRDS(brtM, outModel)
      message("Model saved -> ", outModel)
    }

    if (isValidCachedRDS(outPerf)) {
      message("Loading cached performance...")
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
              file.path(outputDir, paste0(spClean, "_expl_dev_habitat.rds")))
    }

    message("Predicting onto German habitat grid for ", length(predictionYears), " years...")
    spPredFiles <- list()

    for (yr in predictionYears) {
      outTif <- file.path(outputDir, paste0(spClean, "_pred_habitat_", yr, ".tif"))

      if (isValidPredictionRaster(outTif)) {
        message("Year ", yr, ": cache hit")
        spPredFiles[[as.character(yr)]] <- outTif
        next
      }

      covStack <- loadHabitatCovariates(yr, habitatOutputDir)
      if (is.null(covStack)) {
        warning("Year ", yr, ": covariates unavailable -- skipping")
        next
      }

      missingPreds <- setdiff(predSel, names(covStack))
      if (length(missingPreds) > 0) {
        warning("Year ", yr, ": missing predictors: ", paste(missingPreds, collapse = ", "))
        next
      }

      predictBRTToRaster(covStack, predSel, brtM, brtPerf$thresh, outTif)
      message("Year ", yr, ": saved -> ", basename(outTif))
      spPredFiles[[as.character(yr)]] <- outTif
    }

    result[[sp]] <- list(modelPath = outModel, perfPath = outPerf,
                          predictions = spPredFiles, perf = brtPerf)
    message("  == Done: ", sp, " ==============================")
  }

  result
}
