# Test helper: sources the module's functions and builds an uncertainty cfg from the shared config files and the
# BIRDMONITOR_* environment variables (the same values tools/runUncertaintyTask.R passes to the module).
for (f in sort(list.files("modules/models_Monitor/R", pattern = "[.]R$", full.names = TRUE))) source(f)
testCfg <- function(repoRoot = getwd()) {
  setwd(repoRoot)
  env <- new.env(parent = globalenv())
  sys.source("tools/sharedConfig.R", envir = env); sys.source("tools/sharedSpeciesConfig.R", envir = env)
  gen <- if (file.exists("data/speciesConfig_general.csv")) env$loadSpeciesGeneralConfig("data/speciesConfig_general.csv") else NULL
  ev <- function(n, d = "") { v <- Sys.getenv(n, ""); if (nzchar(v)) v else d }
  rng <- function(x) { p <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]]); if (length(p) == 1) p else p[1]:p[2] }
  yrs <- if (nzchar(ev("BIRDMONITOR_UNC_YEARS"))) sort(unique(unlist(lapply(trimws(strsplit(ev("BIRDMONITOR_UNC_YEARS"), ",")[[1]]), rng)))) else env$predictionYears
  uncCfgFromParams(
    inputRoot = file.path(repoRoot, "inputs"), outputRoot = file.path(repoRoot, "outputs", ev("BIRDMONITOR_RUNNAME", "test4")),
    species = env$sharedSpecies, predictionYears = env$predictionYears,
    habitatYears = env$resolveYearsPerSpecies(env$sharedSpecies, "habitat", env$extractYearsConfig(gen), env$sharedHabitatYears),
    resolutionConfig = env$extractResolutionConfig(gen), climateResolutionM = env$sharedClimateResolutionM,
    habitatResolutionM = env$sharedHabitatResolutionM, landscapeResolutionM = env$sharedLandscapeResolutionM,
    climateWindowLength = env$sharedClimateWindowLength, reps = rng(ev("BIRDMONITOR_UNC_REPS", "0:50")), outYears = yrs,
    nBands = as.integer(ev("BIRDMONITOR_UNC_BANDS", "16")), repBatch = as.integer(ev("BIRDMONITOR_UNC_REPBATCH", "10")),
    blockMult = as.numeric(ev("BIRDMONITOR_UNC_BLOCKMULT", "1")), tag = ev("BIRDMONITOR_UNC_TAG", ""),
    baselineYear = as.integer(ev("BIRDMONITOR_UNC_BASELINE", "2005")), codeRoot = repoRoot,
    members = if (nzchar(ev("BIRDMONITOR_UNC_MEMBERS"))) strsplit(ev("BIRDMONITOR_UNC_MEMBERS"), ",")[[1]] else NULL)
}
