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
#' @param cachePath Character, or NULL (default). Directory for
#'   `reproducible::Cache()`'s per-species(-year) cache -- e.g.
#'   `cachePath(sim)`, a stable location shared across runs (NOT the
#'   per-run timestamped output folder). NULL falls back to a temp
#'   directory, for standalone/test calls.
#' @return Named list (by species) with `modelPath`, `perfPath`, `predictions`
#'   (named by year), and `perf` (the evalSDM() row).
modelGerHabitat <- function(inputsData, predictionYears, habitatOutputDir, outputDir,
                             initialLR = 0.08, perSpeciesLR = NULL, cachePath = NULL) {

  if (is.null(cachePath)) cachePath <- file.path(tempdir(), "birdMonitor_cache")
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  result <- list()

  covStacksByYear <- stats::setNames(
    lapply(predictionYears, function(yr) loadHabitatCovariates(yr, habitatOutputDir)),
    as.character(predictionYears))

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

    startingLR <- resolveStartingLR(sp, spClean, defaultLR = initialLR, lrStateDir = outputDir,
                                     perSpeciesLR = perSpeciesLR, lrStateSuffix = "_habitat")
    brtM <- reproducible::Cache(
      fitBRTOneSpecies, sp = sp, spPa = spPa, predSel = predSel, startingLR = startingLR,
      cachePath = cachePath, userTags = c("modelGerHabitat", "fit", spClean))

    if (is.null(brtM)) {
      warning("Skipping ", sp, " at habitat scale -- optimizeBRT() gave up (see its own ",
              "warning above for why). No model saved; this species will simply be ",
              "absent from habitat-scale results until its data issue is fixed.")
      next
    }
    persistConvergedLR(brtM, spClean, lrStateDir = outputDir, lrStateSuffix = "_habitat")
    saveRDS(brtM, outModel)
    message("Model saved -> ", outModel)

    evalResult <- reproducible::Cache(
      evalBRTOneSpecies, sp = sp, spPa = spPa, predSel = predSel, brtM = brtM,
      cachePath = cachePath, userTags = c("modelGerHabitat", "eval", spClean))
    brtPerf <- evalResult$perf
    saveRDS(brtPerf, outPerf)

    message("Performance: AUC = ", round(brtPerf$AUC, 3), " | TSS = ", round(brtPerf$TSS, 3),
            " | D2 = ", round(brtPerf$D2, 3))
    if (brtPerf$AUC < 0.7) warning("AUC < 0.7 for ", sp, " -- interpret with caution.")
    message("D2 at occurrence locations: ", round(evalResult$explDev, 3))
    saveRDS(data.frame(species = sp, D2 = evalResult$explDev),
            file.path(outputDir, paste0(spClean, "_expl_dev_habitat.rds")))

    message("Predicting onto German habitat grid for ", length(predictionYears), " years...")
    spPredFiles <- list()

    for (yr in predictionYears) {
      outTif <- file.path(outputDir, paste0(spClean, "_pred_habitat_", yr, ".tif"))

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
        userTags = c("modelGerHabitat", "predict", spClean, as.character(yr)))
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
