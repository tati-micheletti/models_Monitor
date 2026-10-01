#' Resolve which years a scale-level BRT's own predict loop must cover
#'
#' `modelEurope()`/`modelGerHabitat()`/`modelGerLandscape()`'s own predict
#' loop must ALWAYS include `habitatYears` on top of whatever
#' `predictionYears` the user requested -- `metaModel()` trains against
#' those years' suitability rasters regardless of how restricted
#' `predictionYears` is (see DECISIONS.md's 2026-09-28 "Decouple fitting
#' years from prediction years" entry). `metaModel()`'s OWN predict loop
#' does NOT use this -- it stays exactly `predictionYears`, since its
#' final output should be exactly what was asked for; it can safely
#' assume the scale-level rasters it needs already exist, as they're a
#' superset.
#'
#' @param predictionYears Integer vector. User-requested prediction years.
#' @param habitatYears Integer vector, or named list (species -> integer
#'   vector). Years with real habitat occurrence data -- `metaModel()`'s
#'   training years. A list is flattened to the union across every
#'   species first -- this function covers the scale-level BRTs' own
#'   predict loop, which processes potentially many species in one batch,
#'   so it must cover every species' own training years regardless of
#'   which one requested them.
#' @return Integer vector, the union, sorted.
resolveScalePredictionYears <- function(predictionYears, habitatYears) {
  habitatYearsFlat <- if (is.list(habitatYears)) unlist(habitatYears, use.names = FALSE) else habitatYears
  sort(unique(c(predictionYears, habitatYearsFlat)))
}
