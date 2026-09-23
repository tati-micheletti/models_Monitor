#' Build the metamodel output folder name from the scale labels it combines
#'
#' @param resolutionsM Named numeric vector of scale resolutions (meters)
#'   being combined, e.g. `c(habitat = 200, landscape = 1000, climate = 50000)`.
#' @return Character, e.g. "metamodel_02_1_50" for the vector above --
#'   ordered finest-to-coarsest resolution, independent of input order/names,
#'   so it stays meaningful if a scale is ever added, dropped, or resized.
metamodelLabel <- function(resolutionsM) {
  ordered <- sort(resolutionsM)
  labels <- vapply(ordered, scaleLabel, character(1))
  paste0("metamodel_", paste(sub("^scale_", "", labels), collapse = "_"))
}
