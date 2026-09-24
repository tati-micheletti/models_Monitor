#' Percent change between two years, with a CI that propagates BOTH years'
#' uncertainty -- not the index value's own CI
#'
#' `computeIndexSE()`/`computeCombinedIndexAnalytical()` treat the baseline
#' year as fixed/known, which is fine for the index *level* but means their
#' CI understates uncertainty in the *change* since baseline (it ignores
#' that the baseline year's own spatial mean has uncertainty too) --
#' and, reported as a CI on the absolute index value (e.g. "94.24, CI
#' [94.16, 94.32]"), it doesn't directly answer "how confident are we that
#' this declined by 5.8%", which is the actual question. This function
#' answers that question directly: `pctChange = 100 * (raw_target/raw_base
#' - 1)`, with its SE from the delta-method formula for a ratio of two
#' independent quantities (independent is a reasonable assumption here --
#' the target-year and baseline-year rasters come from separate model
#' predictions, not shared measurements):
#'   Var(A/B) ~= (A/B)^2 * ( Var(A)/A^2 + Var(B)/B^2 )
#' Both years' spatial-averaging uncertainty (`extractMetaProbSD()`) go in,
#' so this CI is wider (more honest) than treating the baseline as exactly
#' known, though still only spatial-averaging uncertainty, not
#' model-fitting or hindcast uncertainty (same caveat as elsewhere).
#'
#' @param rawTarget,rawBase Numeric. Raw (unnormalized) area-mean values at
#'   the target and baseline years (from `extractMetaProbSeries()`).
#' @param seTarget,seBase Numeric. Their standard errors (from
#'   `extractMetaProbSD()`).
#' @return Named numeric vector: `pctChange`, `se`, `lwr95`, `upr95`.
computeChangeCI <- function(rawTarget, rawBase, seTarget, seBase) {
  ratio <- rawTarget / rawBase
  varRatio <- ratio^2 * ((seTarget^2 / rawTarget^2) + (seBase^2 / rawBase^2))
  seRatio <- sqrt(varRatio)

  pctChange <- 100 * (ratio - 1)
  sePct <- 100 * seRatio

  c(pctChange = pctChange, se = sePct,
    lwr95 = pctChange - 1.96 * sePct, upr95 = pctChange + 1.96 * sePct)
}
