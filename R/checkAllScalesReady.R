#' Check whether a species' europe/habitat/landscape predictions all exist
#'
#' Called from all three scale-fitting events (whichever one actually
#' finishes last for a given species triggers the meta-model for it, in
#' the same job) and from the metaModel event itself (as a guard, so a
#' directly-submitted `--scale meta` task -- e.g. a manual re-run for a
#' species some earlier task failed on -- safely no-ops instead of
#' erroring if the other scales aren't ready yet).
#'
#' Checks readiness at `max(predictionYears)` (the user's own requested
#' output) AND every year in `habitatYears` (what `metaModel()` actually
#' trains on) -- not just the former, since restricting `predictionYears`
#' can otherwise make this report "ready" while the habitat-year files
#' `metaModel()` needs are still missing (see DECISIONS.md's 2026-09-28
#' "Decouple fitting years from prediction years" entry).
#'
#' @param species Character. Latin species name.
#' @param outputRoot Character. `outputPath(sim)`.
#' @param climateResolutionM,habitatResolutionM,landscapeResolutionM Numeric.
#'   Passed to `scaleLabel()` to locate each scale's output folder.
#' @param predictionYears Integer vector. User-requested prediction years.
#' @param habitatYears Integer vector. Years with real habitat occurrence
#'   data -- `metaModel()`'s training years.
#' @return Logical.
checkAllScalesReady <- function(species, outputRoot, climateResolutionM,
                                 habitatResolutionM, landscapeResolutionM,
                                 predictionYears, habitatYears) {
  requiredYears <- sort(unique(c(max(predictionYears), habitatYears)))
  all(vapply(requiredYears, function(yr) {
    isSpeciesScaleReady(species, file.path(outputRoot, scaleLabel(climateResolutionM)),
                        "_pred_EU_%d.tif", yr) &&
      isSpeciesScaleReady(species, file.path(outputRoot, scaleLabel(habitatResolutionM)),
                          "_pred_habitat_%d.tif", yr) &&
      isSpeciesScaleReady(species, file.path(outputRoot, scaleLabel(landscapeResolutionM)),
                          "_pred_landscape_%d.tif", yr)
  }, logical(1)))
}
