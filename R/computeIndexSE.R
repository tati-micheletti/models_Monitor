#' Standard error of a baseline-normalized species index, per year
#'
#' Companion to `computeSpeciesIndex()`: since `Index_t = 100 * raw_t /
#' raw_baseline` is a linear rescaling of `raw_t` (treating `raw_baseline`
#' as fixed, per `computeCombinedIndexAnalytical()`'s documented
#' simplification), its standard error scales the same way:
#' `SE(Index_t) = (100 / raw_baseline) * SE(raw_t)`.
#'
#' @param rawSE Named numeric vector from `extractMetaProbSD()` (names =
#'   years as character).
#' @param rawBaseline Numeric. The species' raw value at the baseline year
#'   (i.e. `rawSeries[[as.character(baselineYear)]]`, the same denominator
#'   `computeSpeciesIndex()` used).
#' @return Named numeric vector, same names as `rawSE`, of index SE values.
computeIndexSE <- function(rawSE, rawBaseline) {
  (100 / rawBaseline) * rawSE
}
