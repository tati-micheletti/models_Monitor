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
                   "terra", "dismo", "gbm", "glmnet", "PresenceAbsence"),
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

    ## BRT learning-rate search starting points (per Wiedenroth et al. tuning notes) -----
    defineParameter("europeInitialLR", "numeric", 0.01, NA, NA,
                    "Starting learning rate for the European climate BRT's optimizeBRT() search."),
    defineParameter("habitatInitialLR", "numeric", 0.08, NA, NA,
                    "Starting learning rate for the German habitat BRT's optimizeBRT() search."),
    defineParameter("landscapeInitialLR", "numeric", 0.08, NA, NA,
                    "Starting learning rate for the German landscape BRT's optimizeBRT() search."),

    ## Rerun control ------------------------------------------------------------------
    defineParameter("rerunModelEurope", "logical", FALSE, NA, NA,
                    "Should modelEurope be re-run even if sim$europeModels exists?"),
    defineParameter("rerunModelGerHabitat", "logical", FALSE, NA, NA,
                    "Should modelGerHabitat be re-run even if sim$habitatModels exists?"),
    defineParameter("rerunModelGerLandscape", "logical", FALSE, NA, NA,
                    "Should modelGerLandscape be re-run even if sim$landscapeModels exists?"),
    defineParameter("rerunMetaModel", "logical", FALSE, NA, NA,
                    "Should metaModel be re-run even if sim$metaModels exists?")
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
      sim <- scheduleEvent(sim, time(sim), "models_Monitor", "modelEurope")
      sim <- scheduleEvent(sim, time(sim), "models_Monitor", "modelGerHabitat")
      sim <- scheduleEvent(sim, time(sim), "models_Monitor", "modelGerLandscape")
      sim <- scheduleEvent(sim, time(sim), "models_Monitor", "metaModel")
    },

    modelEurope = {
      # ! ----- EDIT BELOW ----- ! #
      if (is.null(sim$europeModels) || P(sim)$rerunModelEurope) {
        sim$europeModels <- modelEurope(
          inputsData = sim$inputsData$europe,
          climateTargetYears = P(sim)$climateTargetYears,
          climateWindowLength = P(sim)$climateWindowLength,
          climateOutputDir = file.path(outputPath(sim), "climate"),
          outputDir = file.path(outputPath(sim), "models_Monitor", "europe"),
          initialLR = P(sim)$europeInitialLR)
      }
      # ! ----- STOP EDITING ----- ! #
    },

    modelGerHabitat = {
      # ! ----- EDIT BELOW ----- ! #
      if (is.null(sim$habitatModels) || P(sim)$rerunModelGerHabitat) {
        sim$habitatModels <- modelGerHabitat(
          inputsData = sim$inputsData$gerHabitat,
          predictionYears = P(sim)$landscapeYears,
          habitatOutputDir = file.path(outputPath(sim), "habitat"),
          outputDir = file.path(outputPath(sim), "models_Monitor", "habitat"),
          initialLR = P(sim)$habitatInitialLR)
      }
      # ! ----- STOP EDITING ----- ! #
    },

    modelGerLandscape = {
      # ! ----- EDIT BELOW ----- ! #
      if (is.null(sim$landscapeModels) || P(sim)$rerunModelGerLandscape) {
        sim$landscapeModels <- modelGerLandscape(
          inputsData = sim$inputsData$gerLandscape,
          predictionYears = P(sim)$landscapeYears,
          landscapeOutputDir = file.path(outputPath(sim), "landscape"),
          outputDir = file.path(outputPath(sim), "models_Monitor", "landscape"),
          initialLR = P(sim)$landscapeInitialLR)
      }
      # ! ----- STOP EDITING ----- ! #
    },

    metaModel = {
      # ! ----- EDIT BELOW ----- ! #
      if (is.null(sim$metaModels) || P(sim)$rerunMetaModel) {
        if (is.null(sim$europeModels) || is.null(sim$habitatModels) || is.null(sim$landscapeModels)) {
          stop("metaModel requires europeModels/habitatModels/landscapeModels -- check the ",
               "modelEurope/modelGerHabitat/modelGerLandscape events ran first.")
        }

        # Static DEM-derived reference grid -- deliberately not any species'
        # habitat prediction, so this never depends on model output ordering.
        refRaster <- terra::rast(file.path(outputPath(sim), "habitat", "solar_radiation_habitat.tif"))

        sim$metaModels <- metaModel(
          inputsDataGerHabitat = sim$inputsData$gerHabitat,
          habitatYears = P(sim)$habitatYears,
          predictionYears = P(sim)$landscapeYears,
          modelDirs = list(europe = file.path(outputPath(sim), "models_Monitor", "europe"),
                            landscape = file.path(outputPath(sim), "models_Monitor", "landscape"),
                            habitat = file.path(outputPath(sim), "models_Monitor", "habitat")),
          refRaster = refRaster,
          outputDir = file.path(outputPath(sim), "models_Monitor", "meta"))
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
