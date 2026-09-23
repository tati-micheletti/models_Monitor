#' Normalize a species' raw occurrence-probability series to a baseline-100 index
#'
#' `Iₜ = 100 * rawₜ / raw_baseline`, matching how Zbinden et al. (2005)
#' normalize per-species population indices to their chosen baseline year
#' (there: 2000; here: 2005, since that's the first year with genuinely
#' real landscape-scale training data in this pipeline -- see
#' `computeCombinedIndex.R`'s docstring for the reasoning).
#'
#' @param rawSeries Named numeric vector from `extractMetaProbSeries()`
#'   (names = years as character).
#' @param baselineYear Integer. The year whose index value is fixed at 100.
#' @return Named numeric vector, same names as `rawSeries`, rescaled so
#'   `result[as.character(baselineYear)] == 100`. `NA` in -> `NA` out.
computeSpeciesIndex <- function(rawSeries, baselineYear) {
  baseKey <- as.character(baselineYear)
  if (!baseKey %in% names(rawSeries) || is.na(rawSeries[[baseKey]])) {
    stop("No valid baseline-year (", baselineYear, ") value in rawSeries -- ",
         "cannot normalize.")
  }
  100 * rawSeries / rawSeries[[baseKey]]
}
