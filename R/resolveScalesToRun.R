#' Resolve which of the 4 model stages to run for one species
#'
#' `scalesToRun` (see `models_Monitor`'s parameter of the same name) is a
#' named list, species -> character vector naming a subset of
#' "climate"/"habitat"/"landscape"/"meta". It's a testing/debugging control
#' (e.g. fitting only habitat+landscape to compare covariate variants without
#' the climate niche or the final meta-model step), not a persistent
#' ecological fact about a species, so unlike `resolutionConfig` it has no
#' CSV-backed source and defaults to `NULL` -- every species runs all 4
#' stages, exactly today's behavior, unless explicitly overridden.
#'
#' @param species Character. Latin species name.
#' @param scalesToRun The named list, or NULL.
#' @return Character vector, the stages to run for this species.
resolveScalesToRun <- function(species, scalesToRun) {
  allScales <- c("climate", "habitat", "landscape", "meta")
  if (is.null(scalesToRun) || is.null(scalesToRun[[species]])) return(allScales)
  scalesToRun[[species]]
}
