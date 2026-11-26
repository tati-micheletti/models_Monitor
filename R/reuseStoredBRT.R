#' The stored BRT model and performance of a species/scale, instead of a refit (2026-11-26)
#'
#' Used to ADD a year to a finished run without changing it: with `BIRDMONITOR_REUSE_MODELS=1` the landscape and habitat model functions read `<species>_BRT_<scale>.rds` and `<species>_perf_<scale>.rds` from the run
#' folder and predict the requested years with exactly that model. Without the variable (default) nothing changes: the BRT is fitted (or taken from `reproducible::Cache()`) as before.
#' Why it is needed: the cache of a cluster task lives in the job's own temporary folder, so every job refits, and the landscape/habitat BRT is NOT reproducible on a refit (the optimisation picks another number of
#' trees: Alauda landscape 4750 -> 4600, correlation of the maps 0.98, single cells differ by up to 0.75; checked 2026-11-26). A new year predicted with a refitted model would not be comparable with the other years.
#' @param outModel,outPerf Character. Paths of the stored model and performance files.
#' @return NULL (reuse not requested, or one of the files is missing: then the model is fitted as usual), or list(model, perf).
reuseStoredBRT <- function(outModel, outPerf) {
  if (!identical(Sys.getenv("BIRDMONITOR_REUSE_MODELS", ""), "1")) return(NULL)
  if (!file.exists(outModel) || !file.exists(outPerf)) {
    message("BIRDMONITOR_REUSE_MODELS=1 but no stored model/performance (", basename(outModel), "): this species/scale is fitted as usual")
    return(NULL)
  }
  m <- readRDS(outModel); p <- readRDS(outPerf)
  message("REUSE (BIRDMONITOR_REUSE_MODELS=1): stored ", basename(outModel), " and its performance (AUC ", round(p$AUC, 3), "), no refit")
  list(model = m, perf = p)
}
