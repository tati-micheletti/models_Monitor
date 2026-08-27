#' Train the German landscape BRT (1km) and predict onto the full German grid
#'
#' For each species: trains ONE BRT (pooled across years) on
#' `inputsData[[sp]]$data` (inputs_Monitor's final landscape-scale table),
#' block-cross-validates using that table's own `foldID` column, then
#' predicts onto `predictionYears` (the full German modeling period).
#'
#' @param inputsData Named list (by species) with `data`/`predictors`, from
#'   `sim$inputsData$gerLandscape` (inputs_Monitor).
#' @param predictionYears Integer vector of years to predict onto.
#' @param landscapeOutputDir Character. Directory of landscape-scale covariate rasters.
#' @param outputDir Character. Directory to save model/performance/prediction outputs in.
#' @param initialLR Numeric. Starting learning rate for `optimizeBRT()`.
#' @return Named list (by species) with `modelPath`, `perfPath`, `predictions`
#'   (named by year), and `perf` (the evalSDM() row).
modelGerLandscape <- function(inputsData, predictionYears, landscapeOutputDir, outputDir,
                               initialLR = 0.08) {

  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  result <- list()

  for (sp in names(inputsData)) {
    spClean <- gsub(" ", "_", sp)
    message("\n  == ", sp, " ==============================")

    spPa <- inputsData[[sp]]$data
    predSel <- inputsData[[sp]]$predictors

    outModel <- file.path(outputDir, paste0(spClean, "_BRT_landscape.rds"))
    outPerf <- file.path(outputDir, paste0(spClean, "_perf_landscape.rds"))

    message("Records: ", nrow(spPa), " (", sum(spPa$occurrence == 1), " pres / ",
            sum(spPa$occurrence == 0), " abs)")
    message("Predictors (", length(predSel), "): ", paste(predSel, collapse = ", "))

    if (isValidCachedRDS(outModel)) {
      message("Loading cached BRT model...")
      brtM <- readRDS(outModel)
    } else {
      message("Training BRT (optimising learning rate)...")
      brtM <- optimizeBRT(spPa, predSel, "occurrence", initialLR)
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
              file.path(outputDir, paste0(spClean, "_expl_dev_landscape.rds")))
    }

    message("Predicting onto German grid for ", length(predictionYears), " years...")
    spPredFiles <- list()

    for (yr in predictionYears) {
      outTif <- file.path(outputDir, paste0(spClean, "_pred_landscape_", yr, ".tif"))

      if (isValidPredictionRaster(outTif)) {
        message("Year ", yr, ": cache hit")
        spPredFiles[[as.character(yr)]] <- outTif
        next
      }

      covStack <- loadCovariates(yr, landscapeOutputDir, NULL)
      if (is.null(covStack)) {
        warning("Year ", yr, ": covariates unavailable -- skipping")
        next
      }
      names(covStack) <- gsub("_\\d{4}$", "", names(covStack))

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
