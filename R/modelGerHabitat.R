#' Train the German habitat BRT (200m) and predict onto the full German grid
#'
#' For each species: trains ONE BRT (pooled across the 2022-2025 years that
#' habitat occurrence data covers) on `inputsData[[sp]]$data` (inputs_Monitor's
#' final habitat-scale table), block-cross-validates using that table's own
#' `foldID` column, then predicts onto `predictionYears` (by default the full
#' German modeling period, 2005-2025 -- extrapolating beyond the training
#' years' covariate relationships, matching Wiedenroth et al.).
#'
#' Each species' own habitat resolution is resolved from `resolutionConfig`
#' (falling back to `sharedResolutionM` when unset), so covariates and
#' outputs are read/written to that species' own `scale_X` folder -- species
#' sharing a resolution also share one cached covariate stack per year.
#'
#' @param inputsData Named list (by species) with `data`/`predictors`, from
#'   `sim$inputsData$gerHabitat` (inputs_Monitor).
#' @param predictionYears Integer vector of years to predict onto.
#' @param processedRoot Character. `predictors/processed` directory (without
#'   the scale_X leaf -- each species' own leaf is appended internally).
#' @param outputRoot Character. Directory to save model/performance/
#'   prediction outputs in (without the scale_X leaf).
#' @param resolutionConfig Named list, species -> scale -> resolution (m), or
#'   NULL. See `models_Monitor`'s parameter of the same name.
#' @param sharedResolutionM Numeric. Shared default habitat resolution (m),
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
#' @param keepStacksInMemory Logical. FALSE (default): each year's covariate stack is loaded when that year is
#'   predicted and dropped right after. TRUE: all years' stacks stay in memory and are re-used across species. One year
#'   at 200 m is ~2.5 GB (measured), so keeping 21 years needs ~50 GB -- far more than a cluster task has (the first EVE
#'   habitat array died out of memory at 13 GB because of exactly this). Only worth it on a big local machine.
#' @param yearChunk NULL (default): predict every year in this call. Otherwise `c(i, n)`: predict only the i-th of n
#'   equal chunks of `predictionYears`. The cluster runs a species as n SEPARATE R processes (one chunk each) because, on
#'   EVE, process memory kept growing from year to year until the task was killed even with 24 GB; a fresh process every
#'   few years resets it. The model fit and its evaluation are cached, so each chunk starts quickly.
#' @return Named list (by species) with `modelPath`, `perfPath`, `predictions`
#'   (named by year), and `perf` (the evalSDM() row).
modelGerHabitat <- function(inputsData, predictionYears, processedRoot, outputRoot,
                             resolutionConfig = NULL, sharedResolutionM,
                             initialLR = 0.08, perSpeciesLR = NULL, cachePath = NULL,
                             keepStacksInMemory = FALSE, yearChunk = NULL) {

  if (is.null(cachePath)) cachePath <- file.path(tempdir(), "birdMonitor_cache")
  result <- list()
  stackCache <- list()

  for (sp in names(inputsData)) {
    spClean <- gsub(" ", "_", sp)
    message("\n  == ", sp, " ==============================")

    resM <- resolveResolutionM(sp, "habitat", resolutionConfig, sharedResolutionM)
    habitatOutputDir <- file.path(processedRoot, scaleLabel(resM))
    outputDir <- file.path(outputRoot, scaleLabel(resM))
    dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

    resKey <- as.character(resM)
    # One year's stack at a time (see keepStacksInMemory): loading all years first needs ~50 GB at 200 m.
    getStack <- function(yr) {
      key <- paste(resKey, yr)
      if (keepStacksInMemory && !is.null(stackCache[[key]])) return(stackCache[[key]])
      # Preferred: the year's stack as a FILE (written once by the uncertainty workflow's covcache step,
      # <outputRoot>/uncertainty/covcache/<scale>/habitat_<year>.tif -- the same layers loadHabitatCovariates()
      # returns). Reading it needs a few GB; building it in memory with loadHabitatCovariates() peaks at ~11.6 GB
      # (measured) and killed the 20 GB cluster tasks.
      cacheFile <- file.path(outputRoot, "uncertainty", "covcache", scaleLabel(resM), paste0("habitat_", yr, ".tif"))
      s <- if (file.exists(cacheFile)) {
        message("Habitat covariates for ", yr, ": reading the cached stack ", basename(cacheFile))
        terra::rast(cacheFile)
      } else loadHabitatCovariates(yr, habitatOutputDir)
      if (keepStacksInMemory) stackCache[[key]] <<- s
      s
    }

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
    if (!is.null(evalResult$foldModels)) {
      saveRDS(evalResult$foldModels,
              file.path(outputDir, paste0(spClean, "_foldModels_habitat.rds")))
    }
    saveRDS(brtPerf, outPerf)

    message("Performance: AUC = ", round(brtPerf$AUC, 3), " | TSS = ", round(brtPerf$TSS, 3),
            " | D2 = ", round(brtPerf$D2, 3))
    if (brtPerf$AUC < 0.7) warning("AUC < 0.7 for ", sp, " -- interpret with caution.")
    message("D2 at occurrence locations: ", round(evalResult$explDev, 3))
    saveRDS(data.frame(species = sp, D2 = evalResult$explDev),
            file.path(outputDir, paste0(spClean, "_expl_dev_habitat.rds")))

    yearsToDo <- predictionYears
    if (!is.null(yearChunk)) {
      nYears <- length(predictionYears); chunkSize <- ceiling(nYears / yearChunk[2])
      first <- (yearChunk[1] - 1) * chunkSize + 1
      yearsToDo <- if (first > nYears) integer(0) else predictionYears[first:min(nYears, yearChunk[1] * chunkSize)]
      message("Chunk ", yearChunk[1], " of ", yearChunk[2], ": years ", if (length(yearsToDo)) paste(range(yearsToDo), collapse = "-") else "(none)")
    }
    message("Predicting onto German habitat grid for ", length(yearsToDo), " years...")
    spPredFiles <- list()

    for (yr in yearsToDo) {
      outTif <- file.path(outputDir, paste0(spClean, "_pred_habitat_", yr, ".tif"))

      covStack <- getStack(yr)
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
      if (!keepStacksInMemory) { rm(covStack, predRaster); invisible(gc()) }
    }

    result[[sp]] <- list(modelPath = outModel, perfPath = outPerf,
                          predictions = spPredFiles, perf = brtPerf)
    message("  == Done: ", sp, " ==============================")
  }

  result
}
