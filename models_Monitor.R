## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "models_Monitor",
  description = paste("Fits BRT species distribution models per scale (European climate,",
                       "German landscape, German habitat) and combines them via a ridge",
                       "regression meta-model. Consumes inputs_Monitor's sim$inputsData",
                       "directly as training data."),
  keywords = c("bird monitor", "BRT", "species distribution model", "ridge regression", "meta-model"),
  authors = structure(list(list(given = "Tati", family = "Micheletti", role = c("aut", "cre"),
                                 email = "tati.micheletti@gmail.com", comment = NULL),
                           list(given = "Lisa", family = "Hildebrand", role = "aut",
                                email = "lisa.hildebrand@ufz.de", comment = NULL)), class = "person"),
  childModules = character(0),
  version = list(models_Monitor = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "models_Monitor.Rmd"),
  reqdPkgs = list("PredictiveEcology/SpaDES.core@development (>= 3.2.0)",
                   "terra", "dismo", "gbm", "glmnet", "PresenceAbsence", "geodata"),
  parameters = bindrows(
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used - e.g., a hash of the study",
                          "area obtained using `reproducible::studyAreaName()`"),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?"),

    ## Climate prediction range (must match dataPrep_Monitor's values) ------------------
    defineParameter("climateTargetYears", "numeric", 2005:2025, NA, NA,
                    "Target years to predict the European climate BRT onto (one bioclim",
                    "rolling-window raster per year). Must match dataPrep_Monitor's",
                    "climateTargetYears."),
    defineParameter("climateWindowLength", "numeric", 6, NA, NA,
                    "Rolling window length (years) used to compute the bioclim climatology.",
                    "Must match dataPrep_Monitor's climateWindowLength."),

    ## German prediction/training ranges -------------------------------------------------
    defineParameter("landscapeYears", "numeric", 2005:2025, NA, NA,
                    "Years to predict the German habitat and landscape BRTs (and the",
                    "meta-model) onto. Habitat/meta training data covers only habitatYears;",
                    "years outside that range are extrapolations/hindcasts. Must match",
                    "dataPrep_Monitor's/inputs_Monitor's landscapeYears."),
    defineParameter("habitatYears", "numeric", 2022:2025, NA, NA,
                    "Years with real habitat occurrence data -- the meta-model's training",
                    "years. Must match dataPrep_Monitor's/inputs_Monitor's habitatYears."),

    ## Scale resolutions (must match dataPrep_Monitor's/inputs_Monitor's copies) -------
    defineParameter("climateResolutionM", "numeric", 50000, NA, NA,
                    "Resolution (m) of the climate scale -- used, via scaleLabel(), to",
                    "locate that scale's processed covariates and name its inputs/outputs",
                    "subfolders. Must match dataPrep_Monitor's climateResolutionM."),
    defineParameter("habitatResolutionM", "numeric", 200, NA, NA,
                    "Resolution (m) of the habitat scale -- used, via scaleLabel(), to",
                    "locate that scale's processed covariates and name its inputs/outputs",
                    "subfolders. Must match dataPrep_Monitor's habitatResolutionM."),
    defineParameter("landscapeResolutionM", "numeric", 1000, NA, NA,
                    "Resolution (m) of the landscape scale -- used, via scaleLabel(), to",
                    "locate that scale's processed covariates and name its inputs/outputs",
                    "subfolders. Must match dataPrep_Monitor's landscapeResolutionM."),
    defineParameter("landscapeResolutionOverrides", "list", list(), NA, NA,
                    "Named list, species -> resolution (m), e.g.",
                    "list(\"Milvus milvus\" = 10000). Only takes effect for that species'",
                    "single-species cluster task (runSpecies set) -- the batched full run",
                    "still uses one shared landscapeResolutionM for every species (per-",
                    "species resolution in a batched multi-species fit isn't supported;",
                    "see resolveSpeciesResolution()). Requires the override resolution's",
                    "covariates to already exist (dataPrep_Monitor doesn't yet generate",
                    "more than one landscape resolution per run -- see improvements.md",
                    "item 4)."),

    ## BRT learning-rate search starting points (per Wiedenroth et al. tuning notes) -----
    defineParameter("europeInitialLR", "numeric", 0.01, NA, NA,
                    "Starting learning rate for the European climate BRT's optimizeBRT() search."),
    defineParameter("habitatInitialLR", "numeric", 0.08, NA, NA,
                    "Starting learning rate for the German habitat BRT's optimizeBRT() search."),
    defineParameter("landscapeInitialLR", "numeric", 0.08, NA, NA,
                    "Starting learning rate for the German landscape BRT's optimizeBRT() search."),

    ## Cluster-task restriction -- leave both NA for the normal full run; every
    ## default codepath is unchanged when they're NA. See tools/runClusterTask.R. -
    defineParameter("runScale", "character", NA_character_, NA, NA,
                    "NA (default): schedule and run all 4 stages for all species,",
                    "exactly as the original module did. One of \"europe\"/\"habitat\"/",
                    "\"landscape\"/\"meta\": schedule only that stage, for a single",
                    "cluster task. Requires runSpecies to also be set."),
    defineParameter("runSpecies", "character", NA_character_, NA, NA,
                    "NA (default): run every species, exactly as the original module",
                    "did. A single Latin species name: restrict this run to just that",
                    "species, for a single cluster task. Requires runScale to also be set."),

    ## Rerun control ------------------------------------------------------------------
    defineParameter("rerunModelEurope", "logical", FALSE, NA, NA,
                    "Should modelEurope be re-run even if sim$europeModels exists?"),
    defineParameter("rerunModelGerHabitat", "logical", FALSE, NA, NA,
                    "Should modelGerHabitat be re-run even if sim$habitatModels exists?"),
    defineParameter("rerunModelGerLandscape", "logical", FALSE, NA, NA,
                    "Should modelGerLandscape be re-run even if sim$landscapeModels exists?"),
    defineParameter("rerunMetaModel", "logical", FALSE, NA, NA,
                    "Should metaModel be re-run even if sim$metaModels exists?"),
    defineParameter("nBootTrend", "numeric", 0, NA, NA,
                    paste("If > 0, bootstrap the ridge meta-model this many times per",
                    "species to quantify model-fitting uncertainty in the per-year",
                    "area-mean trend (see bootstrapMetaModelTrend()), in addition to",
                    "spatial-averaging precision. 0 (default): off, no added cost."))
  ),
  inputObjects = bindrows(
    expectsInput("inputsData", "list",
                 "List with europe/gerHabitat/gerLandscape named-by-species lists, each with",
                 "`data` (occurrence/coordinates/foldID/resolved predictor columns) and",
                 "`predictors` -- produced by inputs_Monitor.")
  ),
  outputObjects = bindrows(
    createsOutput("europeModels", "list",
                  "Named list (by species) with modelPath/perfPath/predictions/perf for the",
                  "European climate BRT."),
    createsOutput("habitatModels", "list",
                  "Named list (by species) with modelPath/perfPath/predictions/perf for the",
                  "German habitat BRT."),
    createsOutput("landscapeModels", "list",
                  "Named list (by species) with modelPath/perfPath/predictions/perf for the",
                  "German landscape BRT."),
    createsOutput("metaModels", "list",
                  "Named list (by species) with modelPath/perfPath/varimpPath/predictions/",
                  "perf/varimp for the ridge regression meta-model -- the pipeline's final",
                  "per-species, per-year suitability output.")
  )
))

doEvent.models_Monitor = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      if (!is.na(P(sim)$runScale)) {
        if (is.na(P(sim)$runSpecies)) {
          stop("runScale is set but runSpecies is NA -- both must be set together ",
               "for a single cluster task (leave both NA for the normal full run).")
        }
        eventName <- switch(P(sim)$runScale,
                            europe = "modelEurope",
                            habitat = "modelGerHabitat",
                            landscape = "modelGerLandscape",
                            meta = "metaModel",
                            stop("runScale must be one of: europe, habitat, landscape, meta ",
                                 "(got: ", P(sim)$runScale, ")"))
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", eventName)
      } else {
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "modelEurope")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "modelGerHabitat")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "modelGerLandscape")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "metaModel")
      }
    },

    modelEurope = {
      # ! ----- EDIT BELOW ----- ! #
      inputsData <- sim$inputsData$europe
      if (!is.na(P(sim)$runSpecies)) inputsData <- inputsData[P(sim)$runSpecies]

      landscapeResM <- resolveSpeciesResolution(P(sim)$runSpecies, P(sim)$landscapeResolutionM,
                                                 P(sim)$landscapeResolutionOverrides)

      if (is.null(sim$europeModels) || P(sim)$rerunModelEurope) {
        sim$europeModels <- modelEurope(
          inputsData = inputsData,
          climateTargetYears = P(sim)$climateTargetYears,
          climateWindowLength = P(sim)$climateWindowLength,
          climateOutputDir = file.path(inputPath(sim), "predictors", "processed",
                                        scaleLabel(P(sim)$climateResolutionM)),
          outputDir = file.path(outputPath(sim), scaleLabel(P(sim)$climateResolutionM)),
          initialLR = P(sim)$europeInitialLR)
      }

      if (!is.na(P(sim)$runScale) && checkAllScalesReady(
            species = P(sim)$runSpecies, outputRoot = outputPath(sim),
            climateResolutionM = P(sim)$climateResolutionM,
            habitatResolutionM = P(sim)$habitatResolutionM,
            landscapeResolutionM = landscapeResM,
            climateTargetYears = P(sim)$climateTargetYears,
            landscapeYears = P(sim)$landscapeYears)) {
        message(P(sim)$runSpecies, ": all 3 scales ready -- also running metaModel in this task.")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "metaModel")
      }
      # ! ----- STOP EDITING ----- ! #
    },

    modelGerHabitat = {
      # ! ----- EDIT BELOW ----- ! #
      inputsData <- sim$inputsData$gerHabitat
      if (!is.na(P(sim)$runSpecies)) inputsData <- inputsData[P(sim)$runSpecies]

      landscapeResM <- resolveSpeciesResolution(P(sim)$runSpecies, P(sim)$landscapeResolutionM,
                                                 P(sim)$landscapeResolutionOverrides)

      if (is.null(sim$habitatModels) || P(sim)$rerunModelGerHabitat) {
        sim$habitatModels <- modelGerHabitat(
          inputsData = inputsData,
          predictionYears = P(sim)$landscapeYears,
          habitatOutputDir = file.path(inputPath(sim), "predictors", "processed",
                                        scaleLabel(P(sim)$habitatResolutionM)),
          outputDir = file.path(outputPath(sim), scaleLabel(P(sim)$habitatResolutionM)),
          initialLR = P(sim)$habitatInitialLR)
      }

      if (!is.na(P(sim)$runScale) && checkAllScalesReady(
            species = P(sim)$runSpecies, outputRoot = outputPath(sim),
            climateResolutionM = P(sim)$climateResolutionM,
            habitatResolutionM = P(sim)$habitatResolutionM,
            landscapeResolutionM = landscapeResM,
            climateTargetYears = P(sim)$climateTargetYears,
            landscapeYears = P(sim)$landscapeYears)) {
        message(P(sim)$runSpecies, ": all 3 scales ready -- also running metaModel in this task.")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "metaModel")
      }
      # ! ----- STOP EDITING ----- ! #
    },

    modelGerLandscape = {
      # ! ----- EDIT BELOW ----- ! #
      inputsData <- sim$inputsData$gerLandscape
      if (!is.na(P(sim)$runSpecies)) inputsData <- inputsData[P(sim)$runSpecies]

      landscapeResM <- resolveSpeciesResolution(P(sim)$runSpecies, P(sim)$landscapeResolutionM,
                                                 P(sim)$landscapeResolutionOverrides)
      if (landscapeResM != P(sim)$landscapeResolutionM) {
        message(P(sim)$runSpecies, ": using overridden landscape resolution ",
                landscapeResM, "m (default is ", P(sim)$landscapeResolutionM, "m).")
      }

      if (is.null(sim$landscapeModels) || P(sim)$rerunModelGerLandscape) {
        sim$landscapeModels <- modelGerLandscape(
          inputsData = inputsData,
          predictionYears = P(sim)$landscapeYears,
          landscapeOutputDir = file.path(inputPath(sim), "predictors", "processed",
                                          scaleLabel(landscapeResM)),
          outputDir = file.path(outputPath(sim), scaleLabel(landscapeResM)),
          initialLR = P(sim)$landscapeInitialLR)
      }

      if (!is.na(P(sim)$runScale) && checkAllScalesReady(
            species = P(sim)$runSpecies, outputRoot = outputPath(sim),
            climateResolutionM = P(sim)$climateResolutionM,
            habitatResolutionM = P(sim)$habitatResolutionM,
            landscapeResolutionM = landscapeResM,
            climateTargetYears = P(sim)$climateTargetYears,
            landscapeYears = P(sim)$landscapeYears)) {
        message(P(sim)$runSpecies, ": all 3 scales ready -- also running metaModel in this task.")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "metaModel")
      }
      # ! ----- STOP EDITING ----- ! #
    },

    metaModel = {
      # ! ----- EDIT BELOW ----- ! #
      clusterMode <- !is.na(P(sim)$runScale)

      # Which landscape resolution THIS species' own fitted model actually
      # lives under (may be overridden -- see landscapeResolutionOverrides).
      # The metamodel OUTPUT folder name below deliberately does NOT use
      # this -- it stays keyed to the shared default resolutions, so every
      # species' meta output lands in the same folder regardless of any
      # individual species' scale override (downstream reporting reads one
      # shared metaDir for all species).
      landscapeResM <- resolveSpeciesResolution(P(sim)$runSpecies, P(sim)$landscapeResolutionM,
                                                 P(sim)$landscapeResolutionOverrides)

      if (clusterMode) {
        # Reached either self-triggered (the scale event that just finished
        # already confirmed readiness) or via a directly-submitted
        # `--scale meta` task (e.g. a manual re-run for a species some
        # earlier task failed on) -- re-check either way, and no-op rather
        # than error if the other scales genuinely aren't ready yet.
        if (!checkAllScalesReady(
              species = P(sim)$runSpecies, outputRoot = outputPath(sim),
              climateResolutionM = P(sim)$climateResolutionM,
              habitatResolutionM = P(sim)$habitatResolutionM,
              landscapeResolutionM = landscapeResM,
              climateTargetYears = P(sim)$climateTargetYears,
              landscapeYears = P(sim)$landscapeYears)) {
          message(P(sim)$runSpecies, ": not all 3 scales ready yet -- skipping metaModel for now.")
          return(invisible(sim))
        }
      } else if (is.null(sim$europeModels) || is.null(sim$habitatModels) || is.null(sim$landscapeModels)) {
        stop("metaModel requires europeModels/habitatModels/landscapeModels -- check the ",
             "modelEurope/modelGerHabitat/modelGerLandscape events ran first.")
      }

      if (is.null(sim$metaModels) || P(sim)$rerunMetaModel) {
        inputsDataGerHabitat <- sim$inputsData$gerHabitat
        if (!is.na(P(sim)$runSpecies)) inputsDataGerHabitat <- inputsDataGerHabitat[P(sim)$runSpecies]

        # Static DEM-derived reference grid -- deliberately not any species'
        # habitat prediction, so this never depends on model output ordering.
        refRaster <- terra::rast(file.path(inputPath(sim), "predictors", "processed",
                                            scaleLabel(P(sim)$habitatResolutionM),
                                            "solar_radiation_habitat.tif"))

        resolutionsM <- c(europe = P(sim)$climateResolutionM,
                           habitat = P(sim)$habitatResolutionM,
                           landscape = P(sim)$landscapeResolutionM)

        sim$metaModels <- metaModel(
          inputsDataGerHabitat = inputsDataGerHabitat,
          habitatYears = P(sim)$habitatYears,
          predictionYears = P(sim)$landscapeYears,
          modelDirs = list(europe = file.path(outputPath(sim), scaleLabel(P(sim)$climateResolutionM)),
                            landscape = file.path(outputPath(sim), scaleLabel(landscapeResM)),
                            habitat = file.path(outputPath(sim), scaleLabel(P(sim)$habitatResolutionM))),
          refRaster = refRaster,
          outputDir = file.path(outputPath(sim), metamodelLabel(resolutionsM)),
          nBootTrend = P(sim)$nBootTrend,
          gadmCacheDir = file.path(inputPath(sim), "predictors", "raw", "gadm"))
      }
      # ! ----- STOP EDITING ----- ! #
    },

    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
