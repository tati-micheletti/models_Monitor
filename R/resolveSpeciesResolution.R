#' Resolve a species' resolution, honoring a per-species override
#'
#' Some species need a different landscape-scale resolution than the
#' shared default -- e.g. Milvus milvus's real home range (~4.4-6.0km
#' radius) is roughly an order of magnitude larger than the pipeline's
#' default 1km landscape window (see `improvements.md` item 4), while
#' Buteo buteo's is already close to it. `overrides` lets specific
#' species use a different resolution without changing the shared
#' default for everyone else.
#'
#' @param species Character. Latin species name.
#' @param defaultResolutionM Numeric. The shared default resolution (m).
#' @param overrides Named list or NULL, e.g. `list("Milvus milvus" = 10000)`.
#'   `NA`/`NULL`/empty is treated as "no overrides".
#' @return Numeric. The resolution (m) to use for this species.
resolveSpeciesResolution <- function(species, defaultResolutionM, overrides = NULL) {
  if (is.null(overrides) || length(overrides) == 0) return(defaultResolutionM)
  if (species %in% names(overrides)) overrides[[species]] else defaultResolutionM
}
