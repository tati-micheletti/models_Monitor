#' Check whether a species' final-year prediction already exists for one scale
#'
#' Used by the cluster-task self-triggering logic (see models_Monitor.R's
#' modelEurope/modelGerHabitat/modelGerLandscape/metaModel events): after
#' finishing its own scale, a cluster task checks whether the other scales
#' are also done for its species, and if so, runs metaModel() for that
#' species too, in the same job -- no separate SLURM array or manual job
#' dependency needed just for the meta-model step.
#'
#' @param species Character. Latin species name.
#' @param scaleDir Character. That scale's output directory.
#' @param filePattern Character. `sprintf()` pattern with one `%d` for the
#'   year, e.g. `"_pred_EU_%d.tif"`.
#' @param lastYear Integer. The final year that scale predicts onto.
#' @return Logical.
isSpeciesScaleReady <- function(species, scaleDir, filePattern, lastYear) {
  spClean <- gsub(" ", "_", species)
  f <- file.path(scaleDir, paste0(spClean, sprintf(filePattern, lastYear)))
  isValidPredictionRaster(f)
}
