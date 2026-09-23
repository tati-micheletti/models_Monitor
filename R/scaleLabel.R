#' Turn a raster resolution (meters) into a compact, generalizable scale label
#'
#' Derives the folder-name label used throughout the pipeline's `inputs/`
#' and `outputs/` trees from the resolution parameter itself, rather than a
#' hardcoded "habitat"/"landscape"/"climate" string -- so adding, removing,
#' or renaming a scale never requires touching path-construction code.
#' Duplicated identically in dataPrep_Monitor/inputs_Monitor/models_Monitor
#' (see tools/check_duplicated_functions.R at the repo root).
#'
#' @param resolutionM Numeric. Raster resolution in meters.
#' @return Character, e.g. `scaleLabel(200)` -> "scale_02",
#'   `scaleLabel(1000)` -> "scale_1", `scaleLabel(50000)` -> "scale_50".
scaleLabel <- function(resolutionM) {
  paste0("scale_", sub("\\.", "", as.character(resolutionM / 1000)))
}
