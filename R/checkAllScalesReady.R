#' Check whether a species' europe/habitat/landscape predictions all exist
#'
#' Called from all three scale-fitting events (whichever one actually
#' finishes last for a given species triggers the meta-model for it, in
#' the same job) and from the metaModel event itself (as a guard, so a
#' directly-submitted `--scale meta` task -- e.g. a manual re-run for a
#' species some earlier task failed on -- safely no-ops instead of
#' erroring if the other scales aren't ready yet).
#'
#' @param species Character. Latin species name.
#' @param outputRoot Character. `outputPath(sim)`.
#' @param climateResolutionM,habitatResolutionM,landscapeResolutionM Numeric.
#'   Passed to `scaleLabel()` to locate each scale's output folder.
#' @param climateTargetYears,landscapeYears Integer vectors. Used only for
#'   their last element (the final year each scale predicts onto).
#' @return Logical.
checkAllScalesReady <- function(species, outputRoot, climateResolutionM,
                                 habitatResolutionM, landscapeResolutionM,
                                 climateTargetYears, landscapeYears) {
  isSpeciesScaleReady(species, file.path(outputRoot, scaleLabel(climateResolutionM)),
                      "_pred_EU_%d.tif", max(climateTargetYears)) &&
    isSpeciesScaleReady(species, file.path(outputRoot, scaleLabel(habitatResolutionM)),
                        "_pred_habitat_%d.tif", max(landscapeYears)) &&
    isSpeciesScaleReady(species, file.path(outputRoot, scaleLabel(landscapeResolutionM)),
                        "_pred_landscape_%d.tif", max(landscapeYears))
}
