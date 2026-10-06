#' Glue between the SpaDES events of models_Monitor and the uncertainty functions (unc*.R)
#'
#' The uncertainty workflow (option B: spatial-block bootstrap of the BRTs) is split into steps so the cluster can run
#' each as its own job; every step is a models_Monitor EVENT selected with the `runScale` parameter (cluster) or run in
#' order after metaModel when `uncertaintyReps` is set (single session). Method and outputs: UNCERTAINTY.md.

#' Uncertainty configuration from the module's parameters
uncCfgFromSim <- function(sim) {
  p <- P(sim)
  species <- p$uncertaintySpecies
  if (is.null(species) || length(species) == 0 || all(is.na(species))) species <- names(sim$inputsData$gerHabitat)
  if (is.null(species) || length(species) == 0)
    stop("uncertainty: no species -- set the uncertaintySpecies parameter (the full roster, in a stable order).")
  uncCfgFromParams(
    inputRoot = inputPath(sim), outputRoot = outputPath(sim), species = species,
    predictionYears = p$predictionYears, habitatYears = p$habitatYears, resolutionConfig = p$resolutionConfig,
    climateResolutionM = p$climateResolutionM, habitatResolutionM = p$habitatResolutionM,
    landscapeResolutionM = p$landscapeResolutionM, climateWindowLength = p$climateWindowLength,
    reps = p$uncertaintyReps, outYears = p$uncertaintyYears, nBands = p$uncertaintyBands,
    repBatch = p$uncertaintyRepBatch, cores = p$uncertaintyCores, blockMult = p$uncertaintyBlockMult,
    probs = p$uncertaintyProbs, tag = p$uncertaintyTag, baselineYear = p$uncertaintyBaselineYear,
    currentYear = p$uncertaintyCurrentYear, codeRoot = getwd())
}

#' Handle one uncertainty event
#'
#' With `runSpecies` / `runBand` set (cluster task) the event does that one species / band; with them unset it loops over
#' all species / all bands.
uncHandleEvent <- function(sim, eventType) {
  p <- P(sim)
  if (is.null(p$uncertaintyReps))
    stop("Uncertainty event '", eventType, "' needs the uncertaintyReps parameter (e.g. 0:50).")
  cfg <- uncCfgFromSim(sim)
  spRun <- if (!is.na(p$runSpecies)) p$runSpecies else cfg$species
  bandsRun <- if (!is.na(p$runBand)) as.integer(p$runBand) else seq_len(cfg$nBands)
  message("=== uncertainty: ", eventType, " | replicates ", cfg$repLabel, " | species: ",
          if (length(spRun) > 3) paste(length(spRun), "species") else paste(spRun, collapse = ", "),
          if (eventType %in% c("uncertaintyBand", "uncertaintySummarize", "uncertaintyCommunity"))
            paste0(" | band(s): ", if (length(bandsRun) > 3) paste0(min(bandsRun), "-", max(bandsRun)) else paste(bandsRun, collapse = ",")),
          " | cores ", cfg$cores, " ===")

  switch(eventType,
    uncertaintyCovcache = uncCovcache(cfg),
    uncertaintyFit = for (sp in spRun) uncFitSpecies(cfg, sp),
    uncertaintyCoarse = for (sp in spRun) uncCoarseSpecies(cfg, sp),
    uncertaintyOof = for (sp in spRun) uncOofSpecies(cfg, sp),
    uncertaintyRidge = for (sp in spRun) uncRidgeSpecies(cfg, sp),
    uncertaintyBand = for (sp in spRun) {
      ctx <- uncContext(cfg, sp)
      ridge <- readRDS(file.path(uncSpDir(cfg, sp, "ridge"), paste0("ridge_", cfg$repLabel, ".rds")))
      bands <- uncBands(cfg, sp)$bands
      for (k in bandsRun) for (yr in cfg$outYears) {
        t1 <- Sys.time()
        status <- uncPredictBandYear(cfg, sp, ctx, ridge, yr, bands[[k]])
        message(sp, " band ", k, " year ", yr, ": ", status, " (", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), " min)")
      }
    },
    uncertaintySummarize = for (sp in spRun) {
      bands <- uncBands(cfg, sp)$bands
      for (k in bandsRun) uncSummarizeSpeciesBand(cfg, sp, bands[[k]])
    },
    uncertaintyAssemble = for (sp in spRun) uncAssembleSpecies(cfg, sp),
    uncertaintyCommunity = {
      bands <- uncBands(cfg, cfg$species[1])$bands
      for (k in bandsRun) uncCommunityBand(cfg, bands[[k]])
    },
    stop("Unknown uncertainty event: ", eventType)
  )
  invisible(sim)
}
