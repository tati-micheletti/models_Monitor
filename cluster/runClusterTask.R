#' Standalone entry point for one species x scale model-fitting task
#'
#' Runs OUTSIDE SpaDES entirely -- no simList, no `sim$...`. It loads one
#' species' already-persisted `inputs_Monitor` output straight from disk and
#' calls the exact same modelEurope()/modelGerHabitat()/modelGerLandscape()
#' functions the full pipeline uses, writing to the same shared output
#' directories. Because those functions already cache per-species
#' (`isValidCachedRDS`) and per-species-per-year (`isValidPredictionRaster`),
#' any number of these tasks can run concurrently (e.g. as a SLURM job
#' array) and the results are identical to running the sequential loop --
#' just faster, since species/scale combinations don't depend on each other.
#'
#' Requires `inputs_Monitor`'s `collinearityCheck` event to have already run
#' for this `runName` (that's what writes the `_inputs.rds`/`_predictors.rds`
#' files this script reads).
#'
#' Usage (local test, one task):
#'   Rscript runClusterTask.R --scale habitat --index 3 --run-name test1
#' Usage (SLURM array task -- index comes from $SLURM_ARRAY_TASK_ID if
#' --index is omitted):
#'   Rscript runClusterTask.R --scale habitat --run-name test1
#'
#' Run from the birdMonitor repo root (paths below are relative to it).

## Avoid oversubscription if the allocated node has more cores than this
## task is given -- gbm/dismo don't multithread internally, but the raster
## (GDAL) calls in prediction can, so pin this explicitly regardless.
Sys.setenv(OMP_NUM_THREADS = "1", GDAL_NUM_THREADS = "1")

## ---- Canonical species list -------------------------------------------
## Must match runMe.R's `sharedSpecies`, same order -- `--index` is a
## 1-based position in this vector. Change one, change the other.
sharedSpecies <- c("Vanellus vanellus", "Milvus milvus", "Lanius collurio",
                    "Lullula arborea", "Alauda arvensis", "Saxicola rubetra",
                    "Emberiza calandra", "Emberiza citrinella", "Buteo buteo",
                    "Sturnus vulgaris", "Perdix perdix")

## Must match runMe.R's shared*Years values. (climateTargetYears and
## landscapeYears currently hold the same range, 2005:2025, but are kept as
## separate variables -- as in runMe.R -- since they parameterize different
## things and aren't guaranteed to stay equal.)
sharedClimateWindowLength <- 6
sharedClimateTargetYears <- 2005:2025
sharedLandscapeYears <- 2005:2025  # habitat + landscape + meta prediction years

## Must match runMe.R's shared*ResolutionM values -- these drive scaleLabel()
## folder names throughout inputs/ and outputs/.
sharedClimateResolutionM <- 50000
sharedHabitatResolutionM <- 200
sharedLandscapeResolutionM <- 1000

## ---- Required packages (normally loaded via models_Monitor's reqdPkgs;
## this script bypasses SpaDES, so they're loaded explicitly here) --------
for (pkg in c("terra", "dismo", "gbm", "glmnet", "PresenceAbsence")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed. Install it before submitting ",
         "cluster jobs (see the WKDV ticket's software question).")
  }
}

## ---- CLI argument parsing (no extra package dependency) ----------------
parseArgs <- function(args) {
  getArg <- function(flag, default = NULL) {
    i <- which(args == flag)
    if (length(i) == 0) return(default)
    args[i + 1]
  }
  slurmIdx <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
  list(
    scale    = getArg("--scale"),
    index    = as.integer(getArg("--index", slurmIdx)),
    runName  = getArg("--run-name", "test1"),
    repoRoot = getArg("--repo-root", getwd())
  )
}

opt <- parseArgs(commandArgs(trailingOnly = TRUE))

if (is.null(opt$scale) || !opt$scale %in% c("europe", "habitat", "landscape")) {
  stop("--scale must be one of: europe, habitat, landscape (got: ", opt$scale, ")")
}
if (is.na(opt$index) || opt$index < 1 || opt$index > length(sharedSpecies)) {
  stop("--index (or $SLURM_ARRAY_TASK_ID) must be an integer between 1 and ",
       length(sharedSpecies), " -- got: ", opt$index)
}

species <- sharedSpecies[opt$index]
spClean <- gsub(" ", "_", species)
inputRoot <- file.path(opt$repoRoot, "inputs")
outputRoot <- file.path(opt$repoRoot, "outputs", opt$runName)

message("=== Cluster task: scale = ", opt$scale, " | species = ", species,
        " (index ", opt$index, "/", length(sharedSpecies), ") ===")

## ---- Load this module's functions (identical code the full SpaDES run uses) ----
## Must happen before scaleLabel() is called below.
moduleRDir <- file.path(opt$repoRoot, "modules", "models_Monitor", "R")
rFiles <- list.files(moduleRDir, pattern = "\\.R$", full.names = TRUE)
rFiles <- rFiles[basename(rFiles) != "runClusterTask.R"]
invisible(lapply(rFiles, source))

scaleLabelFor <- c(europe = scaleLabel(sharedClimateResolutionM),
                    habitat = scaleLabel(sharedHabitatResolutionM),
                    landscape = scaleLabel(sharedLandscapeResolutionM))[[opt$scale]]

## ---- Load just this one species' persisted inputs_Monitor output -------
inputsDir <- file.path(inputRoot, "model_ready", scaleLabelFor)

dataFile <- file.path(inputsDir, paste0(spClean, "_inputs.rds"))
predFile <- file.path(inputsDir, paste0(spClean, "_predictors.rds"))
if (!file.exists(dataFile) || !file.exists(predFile)) {
  stop("Missing inputs_Monitor output for ", species, " at ", inputsDir,
       " -- has inputs_Monitor's collinearityCheck event run for runName '",
       opt$runName, "'?")
}

inputsData <- setNames(
  list(list(data = readRDS(dataFile), predictors = readRDS(predFile))),
  species
)

## ---- Dispatch to the right scale's fitting function ---------------------
predictorsProcessedDir <- file.path(inputRoot, "predictors", "processed", scaleLabelFor)

result <- switch(
  opt$scale,
  europe = modelEurope(
    inputsData = inputsData,
    climateTargetYears = sharedClimateTargetYears,
    climateWindowLength = sharedClimateWindowLength,
    climateOutputDir = predictorsProcessedDir,
    outputDir = file.path(outputRoot, scaleLabelFor),
    initialLR = 0.01),
  habitat = modelGerHabitat(
    inputsData = inputsData,
    predictionYears = sharedLandscapeYears,
    habitatOutputDir = predictorsProcessedDir,
    outputDir = file.path(outputRoot, scaleLabelFor),
    initialLR = 0.08),
  landscape = modelGerLandscape(
    inputsData = inputsData,
    predictionYears = sharedLandscapeYears,
    landscapeOutputDir = predictorsProcessedDir,
    outputDir = file.path(outputRoot, scaleLabelFor),
    initialLR = 0.08)
)

message("=== Done: ", opt$scale, " / ", species, " -- AUC = ",
        round(result[[species]]$perf$AUC, 3), " ===")
