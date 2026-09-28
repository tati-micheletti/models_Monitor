#' Resolve a species' resolution (m) for one scale, from `resolutionConfig`
#'
#' `resolutionConfig` (see `models_Monitor`'s parameter of the same name,
#' produced by `extractResolutionConfig()` in `sharedSpeciesConfig.R` from
#' `speciesConfig_general.csv`'s `resolution_m` column) is keyed
#' species -> scale -> resolution (m). A species with no entry, an NA entry,
#' or `species` itself being NA (e.g. a batched run with no single current
#' species) all fall through to `sharedDefault`.
#'
#' @param species Character. Latin species name, or `NA_character_`.
#' @param scale Character, one of "climate"/"habitat"/"landscape".
#' @param resolutionConfig The nested list, or NULL.
#' @param sharedDefault Numeric. The scale's shared default resolution (m).
#' @return Numeric. The resolution (m) to use for this species.
resolveResolutionM <- function(species, scale, resolutionConfig, sharedDefault) {
  if (is.null(resolutionConfig) || is.na(species)) return(sharedDefault)
  v <- resolutionConfig[[species]][[scale]]
  if (is.null(v) || is.na(v)) sharedDefault else v
}
