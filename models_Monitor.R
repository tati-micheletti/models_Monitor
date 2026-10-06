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
                   "PredictiveEcology/reproducible@development",
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

    ## Prediction years (must match dataPrep_Monitor's predictionYears) -----------------
    defineParameter("predictionYears", "numeric", NA_real_, NA, NA,
                    "Years to generate PREDICTION rasters for, across the climate/habitat/",
                    "landscape BRTs and the meta-model. No default -- must be supplied",
                    "explicitly (e.g. predictionYears from sharedConfig.R), so a caller can",
                    "never silently fall back to a stale range. Independent of each scale's",
                    "own FITTING years (dataPrep_Monitor's landscapeYears/habitatYears --",
                    "training data is NEVER filtered by this parameter, see",
                    "DECISIONS.md's 2026-09-28 \"Decouple fitting years from prediction",
                    "years\" entry). The European/habitat/landscape BRTs' own internal",
                    "predict loop auto-covers habitatYears on top of this (see",
                    "resolveScalePredictionYears()) so metaModel() always finds the",
                    "suitability rasters it needs to train on, regardless of how",
                    "restricted this is; metaModel()'s own output stays exactly this."),
    defineParameter("climateWindowLength", "numeric", 6, NA, NA,
                    "Rolling window length (years) used to compute the bioclim climatology.",
                    "Must match dataPrep_Monitor's climateWindowLength."),
    defineParameter("habitatYears", "list", NULL, NA, NA,
                    "Named list, species -> integer vector of years with real habitat",
                    "occurrence data for that species -- both habitat's own fitting-year",
                    "constraint (in dataPrep_Monitor) AND the meta-model's training years",
                    "here (the same underlying data-availability fact, not two separate",
                    "parameters -- see DECISIONS.md's 2026-09-28 entry). Per-species since",
                    "2026-10-01 (e.g. Buteo buteo/Sturnus vulgaris's real MhB point-count",
                    "data is negligible before ~2020, while other species genuinely span a",
                    "wider range) -- see resolveYearsPerSpecies() (sharedSpeciesConfig.R).",
                    "Must match dataPrep_Monitor's/inputs_Monitor's habitatYears. No",
                    "default -- runMe.R must always supply a fully-resolved list."),

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
    defineParameter("resolutionConfig", "list", NULL, NA, NA,
                    "NULL (default): every species uses the shared *ResolutionM parameters",
                    "above. Otherwise a named list, species -> scale -> resolution (m),",
                    "produced by extractResolutionConfig() (sharedSpeciesConfig.R, repo",
                    "root) from speciesConfig_general.csv's resolution_m column -- e.g. a",
                    "species with a coarser habitat range than the shared default (Milvus",
                    "milvus at a 15km landscape window instead of the shared 1km -- see",
                    "DECISIONS.md's 2026-09-28 entry). Resolved once by the orchestrating",
                    "script and passed in as a plain value, same pattern as sharedConfig.R's",
                    "other shared values. See resolveResolutionM()."),

    defineParameter("scalesToRun", "list", NULL, NA, NA,
                    "NULL (default): every species runs all 4 stages (climate/habitat/",
                    "landscape/meta), exactly today's behavior. Otherwise a named list,",
                    "species -> character vector naming a subset of \"climate\"/\"habitat\"/",
                    "\"landscape\"/\"meta\" -- e.g. list(\"Emberiza citrinella\" =",
                    "c(\"habitat\", \"landscape\")) to compare covariate variants without",
                    "the climate niche or the final meta-model step. A testing/debugging",
                    "control, not a persistent species fact -- unlike resolutionConfig, it",
                    "has no CSV-backed source and doesn't need setting in runMe.R for a",
                    "normal full run. See resolveScalesToRun()."),

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
    defineParameter("habitatYearChunk", "numeric", NULL, NA, NA,
                    "NULL (default): the habitat step predicts every year in one process. c(i, n): predict only the i-th",
                    "of n equal chunks of the prediction years (cluster tasks run n processes in a row, so memory resets)."),
    defineParameter("keepHabitatStacks", "logical", FALSE, NA, NA,
                    "FALSE (default): the 200 m habitat covariate stack is loaded one year at a time (~2.5 GB each).",
                    "TRUE: all years stay in memory and are shared across species (~50 GB; only for a big local machine)."),
    defineParameter("nBootTrend", "numeric", 0, NA, NA,
                    paste("If > 0, bootstrap the ridge meta-model this many times per",
                    "species to quantify model-fitting uncertainty in the per-year",
                    "area-mean trend (see bootstrapMetaModelTrend()), in addition to",
                    "spatial-averaging precision. 0 (default): off, no added cost.")),

    ## Uncertainty (option B: spatial-block bootstrap of the BRTs) -- see UNCERTAINTY.md. Off unless
    ## uncertaintyReps is set. Results go to <outputPath>/uncertainty[_<tag>]/. ----------------------
    defineParameter("uncertaintyReps", "numeric", NULL, NA, NA,
                    "NULL (default): no uncertainty analysis. Otherwise the replicate ids of THIS run, e.g. 0:50.",
                    "Replicate 0 is the main models run through the same machinery (a built-in consistency",
                    "check against the baseline maps; left out of every interval). Adding replicates later: a",
                    "new run with disjoint ids (51:100); nothing already computed changes. Single session: the",
                    "uncertainty events then run after metaModel. Cluster: select ONE step per task with",
                    "runScale = uncertaintyCovcache/Fit/Coarse/Ridge/Band/Summarize/Assemble/Community."),
    defineParameter("uncertaintySpecies", "character", NULL, NA, NA,
                    "The full species roster, in the order the cluster arrays index it (NULL: the species of",
                    "sim$inputsData). Needed by the steps that look at all species (community maps)."),
    defineParameter("uncertaintyYears", "numeric", NULL, NA, NA,
                    "Years to map (NULL: all predictionYears). The habitat training years are always computed."),
    defineParameter("uncertaintyBands", "numeric", 16, NA, NA,
                    "Number of horizontal bands the country is cut into (one cluster task per species x band)."),
    defineParameter("uncertaintyRepBatch", "numeric", 10, NA, NA,
                    "Replicates held in memory at once while predicting a band."),
    defineParameter("uncertaintyCores", "numeric", 1, NA, NA,
                    "Cores for forked prediction/fitting inside one task (ignored on Windows)."),
    defineParameter("uncertaintyBlockMult", "numeric", 1, NA, NA,
                    "Multiplier on the resampling block size (the cross-validation block size): sensitivity test."),
    defineParameter("uncertaintyProbs", "numeric", c(0.05, 0.95), NA, NA,
                    "Percentiles of the interval (default 5th-95th = a 90% interval)."),
    defineParameter("uncertaintyTag", "character", "", NA, NA,
                    "Non-empty: write to uncertainty_<tag>/ instead of uncertainty/ (timing and test runs)."),
    defineParameter("uncertaintyBaselineYear", "numeric", 2005, NA, NA,
                    "Baseline year of the change maps (must equal runIndex_Monitor's baselineYear)."),
    defineParameter("uncertaintyCurrentYear", "numeric", NA_real_, NA, NA,
                    "Report year of the change maps (NA: max habitat year, as runMe.R's currentYear)."),
    defineParameter("runBand", "numeric", NA_real_, NA, NA,
                    "Cluster task of the band-wise uncertainty steps: the one band (1..uncertaintyBands) to do.",
                    "NA: all bands.")
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
      if (identical(P(sim)$predictionYears, NA_real_)) {
        stop("models_Monitor's predictionYears parameter must be supplied explicitly ",
             "(e.g. predictionYears from sharedConfig.R) -- no default.")
      }
      uncertaintyEvents <- c("uncertaintyCovcache", "uncertaintyFit", "uncertaintyCoarse", "uncertaintyRidge",
                             "uncertaintyBand", "uncertaintySummarize", "uncertaintyAssemble", "uncertaintyCommunity")
      if (!is.na(P(sim)$runScale)) {
        if (is.na(P(sim)$runSpecies) && !(P(sim)$runScale %in% uncertaintyEvents)) {
          stop("runScale is set but runSpecies is NA -- both must be set together ",
               "for a single cluster task (leave both NA for the normal full run).")
        }
        eventName <- if (P(sim)$runScale %in% uncertaintyEvents) P(sim)$runScale else
          switch(P(sim)$runScale,
                 europe = "modelEurope",
                 habitat = "modelGerHabitat",
                 landscape = "modelGerLandscape",
                 meta = "metaModel",
                 stop("runScale must be one of: europe, habitat, landscape, meta, ",
                      paste(uncertaintyEvents, collapse = ", "), " (got: ", P(sim)$runScale, ")"))
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", eventName)
      } else {
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "modelEurope")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "modelGerHabitat")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "modelGerLandscape")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "metaModel")
        # Uncertainty (option B) only when asked for: runs after metaModel, in this order.
        if (!is.null(P(sim)$uncertaintyReps)) {
          for (ev in uncertaintyEvents) sim <- scheduleEvent(sim, time(sim), "models_Monitor", ev)
        }
      }
    },

    modelEurope = {
      # ! ----- EDIT BELOW ----- ! #
      inputsData <- sim$inputsData$europe
      if (!is.na(P(sim)$runSpecies)) inputsData <- inputsData[P(sim)$runSpecies]
      inputsData <- inputsData[vapply(names(inputsData), function(sp)
        "climate" %in% resolveScalesToRun(sp, P(sim)$scalesToRun), logical(1))]

      # Scale-level BRTs' own predict loop auto-covers habitatYears on top
      # of the user-requested predictionYears, so metaModel() always finds
      # the suitability rasters it needs to train on (see
      # resolveScalePredictionYears(); DECISIONS.md's 2026-09-28 entry).
      scalePredictionYears <- resolveScalePredictionYears(P(sim)$predictionYears, P(sim)$habitatYears)

      if (is.null(sim$europeModels) || P(sim)$rerunModelEurope) {
        sim$europeModels <- modelEurope(
          inputsData = inputsData,
          climateTargetYears = scalePredictionYears,
          climateWindowLength = P(sim)$climateWindowLength,
          climateOutputDir = file.path(inputPath(sim), "predictors", "processed",
                                        scaleLabel(P(sim)$climateResolutionM)),
          outputDir = file.path(outputPath(sim), scaleLabel(P(sim)$climateResolutionM)),
          initialLR = P(sim)$europeInitialLR,
          cachePath = cachePath(sim))
      }

      if (!is.na(P(sim)$runScale) && checkAllScalesReady(
            species = P(sim)$runSpecies, outputRoot = outputPath(sim),
            climateResolutionM = P(sim)$climateResolutionM,
            habitatResolutionM = resolveResolutionM(P(sim)$runSpecies, "habitat",
                                                     P(sim)$resolutionConfig, P(sim)$habitatResolutionM),
            landscapeResolutionM = resolveResolutionM(P(sim)$runSpecies, "landscape",
                                                       P(sim)$resolutionConfig, P(sim)$landscapeResolutionM),
            predictionYears = P(sim)$predictionYears,
            habitatYears = if (is.list(P(sim)$habitatYears)) P(sim)$habitatYears[[P(sim)$runSpecies]] else P(sim)$habitatYears)) {
        message(P(sim)$runSpecies, ": all 3 scales ready -- also running metaModel in this task.")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "metaModel")
      }
      # ! ----- STOP EDITING ----- ! #
    },

    modelGerHabitat = {
      # ! ----- EDIT BELOW ----- ! #
      inputsData <- sim$inputsData$gerHabitat
      if (!is.na(P(sim)$runSpecies)) inputsData <- inputsData[P(sim)$runSpecies]
      inputsData <- inputsData[vapply(names(inputsData), function(sp)
        "habitat" %in% resolveScalesToRun(sp, P(sim)$scalesToRun), logical(1))]

      scalePredictionYears <- resolveScalePredictionYears(P(sim)$predictionYears, P(sim)$habitatYears)

      if (is.null(sim$habitatModels) || P(sim)$rerunModelGerHabitat) {
        sim$habitatModels <- modelGerHabitat(
          inputsData = inputsData,
          predictionYears = scalePredictionYears,
          processedRoot = file.path(inputPath(sim), "predictors", "processed"),
          outputRoot = outputPath(sim),
          resolutionConfig = P(sim)$resolutionConfig,
          sharedResolutionM = P(sim)$habitatResolutionM,
          initialLR = P(sim)$habitatInitialLR,
          cachePath = cachePath(sim),
          keepStacksInMemory = P(sim)$keepHabitatStacks,
          yearChunk = P(sim)$habitatYearChunk)
      }

      if (!is.na(P(sim)$runScale) && checkAllScalesReady(
            species = P(sim)$runSpecies, outputRoot = outputPath(sim),
            climateResolutionM = P(sim)$climateResolutionM,
            habitatResolutionM = resolveResolutionM(P(sim)$runSpecies, "habitat",
                                                     P(sim)$resolutionConfig, P(sim)$habitatResolutionM),
            landscapeResolutionM = resolveResolutionM(P(sim)$runSpecies, "landscape",
                                                       P(sim)$resolutionConfig, P(sim)$landscapeResolutionM),
            predictionYears = P(sim)$predictionYears,
            habitatYears = if (is.list(P(sim)$habitatYears)) P(sim)$habitatYears[[P(sim)$runSpecies]] else P(sim)$habitatYears)) {
        message(P(sim)$runSpecies, ": all 3 scales ready -- also running metaModel in this task.")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "metaModel")
      }
      # ! ----- STOP EDITING ----- ! #
    },

    modelGerLandscape = {
      # ! ----- EDIT BELOW ----- ! #
      inputsData <- sim$inputsData$gerLandscape
      if (!is.na(P(sim)$runSpecies)) inputsData <- inputsData[P(sim)$runSpecies]
      inputsData <- inputsData[vapply(names(inputsData), function(sp)
        "landscape" %in% resolveScalesToRun(sp, P(sim)$scalesToRun), logical(1))]

      scalePredictionYears <- resolveScalePredictionYears(P(sim)$predictionYears, P(sim)$habitatYears)

      if (is.null(sim$landscapeModels) || P(sim)$rerunModelGerLandscape) {
        sim$landscapeModels <- modelGerLandscape(
          inputsData = inputsData,
          predictionYears = scalePredictionYears,
          processedRoot = file.path(inputPath(sim), "predictors", "processed"),
          outputRoot = outputPath(sim),
          resolutionConfig = P(sim)$resolutionConfig,
          sharedResolutionM = P(sim)$landscapeResolutionM,
          initialLR = P(sim)$landscapeInitialLR,
          cachePath = cachePath(sim))
      }

      if (!is.na(P(sim)$runScale) && checkAllScalesReady(
            species = P(sim)$runSpecies, outputRoot = outputPath(sim),
            climateResolutionM = P(sim)$climateResolutionM,
            habitatResolutionM = resolveResolutionM(P(sim)$runSpecies, "habitat",
                                                     P(sim)$resolutionConfig, P(sim)$habitatResolutionM),
            landscapeResolutionM = resolveResolutionM(P(sim)$runSpecies, "landscape",
                                                       P(sim)$resolutionConfig, P(sim)$landscapeResolutionM),
            predictionYears = P(sim)$predictionYears,
            habitatYears = if (is.list(P(sim)$habitatYears)) P(sim)$habitatYears[[P(sim)$runSpecies]] else P(sim)$habitatYears)) {
        message(P(sim)$runSpecies, ": all 3 scales ready -- also running metaModel in this task.")
        sim <- scheduleEvent(sim, time(sim), "models_Monitor", "metaModel")
      }
      # ! ----- STOP EDITING ----- ! #
    },

    metaModel = {
      # ! ----- EDIT BELOW ----- ! #
      clusterMode <- !is.na(P(sim)$runScale)

      # Which habitat/landscape resolution THIS species' own fitted models
      # actually live under (may be overridden -- see resolutionConfig). The
      # metamodel OUTPUT folder name below deliberately does NOT use these --
      # it stays keyed to the shared default resolutions, so every species'
      # meta output lands in the same folder regardless of any individual
      # species' scale override (downstream reporting reads one shared
      # metaDir for all species).
      landscapeResM <- resolveResolutionM(P(sim)$runSpecies, "landscape",
                                           P(sim)$resolutionConfig, P(sim)$landscapeResolutionM)
      habitatResM <- resolveResolutionM(P(sim)$runSpecies, "habitat",
                                         P(sim)$resolutionConfig, P(sim)$habitatResolutionM)

      if (clusterMode) {
        # Reached either self-triggered (the scale event that just finished
        # already confirmed readiness) or via a directly-submitted
        # `--scale meta` task (e.g. a manual re-run for a species some
        # earlier task failed on) -- re-check either way, and no-op rather
        # than error if the other scales genuinely aren't ready yet.
        if (!checkAllScalesReady(
              species = P(sim)$runSpecies, outputRoot = outputPath(sim),
              climateResolutionM = P(sim)$climateResolutionM,
              habitatResolutionM = habitatResM,
              landscapeResolutionM = landscapeResM,
              predictionYears = P(sim)$predictionYears,
              habitatYears = if (is.list(P(sim)$habitatYears)) P(sim)$habitatYears[[P(sim)$runSpecies]] else P(sim)$habitatYears)) {
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
        inputsDataGerHabitat <- inputsDataGerHabitat[vapply(names(inputsDataGerHabitat), function(sp)
          "meta" %in% resolveScalesToRun(sp, P(sim)$scalesToRun), logical(1))]

        # Static DEM-derived reference grid -- deliberately not any species'
        # habitat prediction, so this never depends on model output ordering.
        # Resolution appended to the filename (second safety layer beyond
        # the containing scaleLabel()-named folder).
        habitatResLabel <- scaleLabel(P(sim)$habitatResolutionM)
        refRaster <- terra::rast(file.path(inputPath(sim), "predictors", "processed",
                                            habitatResLabel,
                                            paste0("solar_radiation_habitat_", habitatResLabel, ".tif")))

        resolutionsM <- c(europe = P(sim)$climateResolutionM,
                           habitat = P(sim)$habitatResolutionM,
                           landscape = P(sim)$landscapeResolutionM)

        sim$metaModels <- metaModel(
          inputsDataGerHabitat = inputsDataGerHabitat,
          habitatYears = P(sim)$habitatYears,
          predictionYears = P(sim)$predictionYears,
          modelDirs = list(europe = file.path(outputPath(sim), scaleLabel(P(sim)$climateResolutionM)),
                            landscape = file.path(outputPath(sim), scaleLabel(landscapeResM)),
                            habitat = file.path(outputPath(sim), scaleLabel(habitatResM))),
          refRaster = refRaster,
          outputDir = file.path(outputPath(sim), metamodelLabel(resolutionsM)),
          nBootTrend = P(sim)$nBootTrend,
          gadmCacheDir = file.path(inputPath(sim), "predictors", "raw", "gadm"),
          cachePath = cachePath(sim))
      }
      # ! ----- STOP EDITING ----- ! #
    },

    uncertaintyCovcache = , uncertaintyFit = , uncertaintyCoarse = , uncertaintyRidge = ,
    uncertaintyBand = , uncertaintySummarize = , uncertaintyAssemble = , uncertaintyCommunity = {
      # ! ----- EDIT BELOW ----- ! #
      sim <- uncHandleEvent(sim, eventType)   # see R/uncSim.R and UNCERTAINTY.md
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
