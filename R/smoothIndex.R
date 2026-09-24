#' Smooth an index time series for trend reporting
#'
#' PECBMS/DDA and the modern Living Planet Index both publish two curves
#' per indicator: the raw (noisy) year-to-year index, and a smoothed trend
#' line used for the headline "up/down/stable" read -- LPI has moved to a
#' per-population GAM for this; PECBMS/DDA apply smoothing at the
#' multi-species level. This does the same at whichever level it's called
#' on (a single combined-index series), using `stats::loess()` rather than
#' a GAM so it needs no package beyond base R (`stats::loess` and
#' `mgcv::gam` produce visually similar smooths for a single 21-point
#' series like ours; swap in `mgcv::gam(y ~ s(x))` here if a
#' publication-standard GAM is specifically wanted later).
#'
#' @param years Numeric vector of years (x-axis), same length as `values`.
#' @param values Numeric vector, the index series to smooth. `NA`s allowed.
#' @param span Numeric. `loess()` span (smoothing bandwidth) -- larger is
#'   smoother. Default 0.75 is `loess()`'s own default.
#' @return Numeric vector, same length as `values`: the smoothed series
#'   (`NA` wherever `values` was `NA`, or everywhere if fewer than 4
#'   non-NA points -- too few to fit a local regression).
smoothIndex <- function(years, values, span = 0.75) {
  valid <- !is.na(values)
  smoothed <- rep(NA_real_, length(values))

  if (sum(valid) < 4) return(smoothed)

  fit <- stats::loess(values[valid] ~ years[valid], span = span, degree = 2)
  smoothed[valid] <- as.numeric(stats::predict(fit))
  smoothed
}
