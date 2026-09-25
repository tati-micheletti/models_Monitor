#' Pull one scale's per-species starting learning rates out of the
#' per-species, per-scale general config
#'
#' `perSpeciesGeneralConfig` (see `models_Monitor`'s parameter of the same
#' name, sourced from `loadSpeciesGeneralConfig()`) is keyed
#' species -> scale -> named list of settings, of which only `brt_start_lr`
#' is relevant here. Each scale's `modelX()` function just needs a plain
#' species -> LR vector for its own scale, not the whole nested structure.
#'
#' @param perSpeciesGeneralConfig The nested list, or NULL.
#' @param scale Character, one of "climate"/"landscape"/"habitat".
#' @return Named numeric vector (species -> brt_start_lr), or NULL if
#'   `perSpeciesGeneralConfig` is NULL. A species missing this scale, or
#'   with a blank/NA `brt_start_lr`, is simply absent from the result (so
#'   `resolveStartingLR()` falls through to its own default/persisted value).
extractScaleStartingLR <- function(perSpeciesGeneralConfig, scale) {
  if (is.null(perSpeciesGeneralConfig)) return(NULL)
  lrs <- sapply(perSpeciesGeneralConfig, function(sp) {
    val <- sp[[scale]]$brt_start_lr
    if (is.null(val)) NA_real_ else val
  })
  lrs[!is.na(lrs)]
}
