#' Combine per-species indices via the current MSI-tool Monte Carlo method
#'
#' The DDA/PECBMS/MSI-tool combined index has TWO published CI methods, from
#' two different papers, and this repo now implements BOTH side by side
#' (deliberately -- see `computeCombinedIndexAnalytical()`'s docstring for
#' the older one):
#' - **Old** (Gregory et al. 2005, Appendix A): analytical delta-method
#'   variance formula. `computeCombinedIndexAnalytical()`.
#' - **Current standard** (Soldaat, Pannekoek, Verweij, van Turnhout & van
#'   Strien 2017, "A Monte Carlo method to account for sampling error in
#'   multi-species indicators", Ecological Indicators 81:340-347 -- the
#'   MSI-tool distributed by Statistics Netherlands/CBS and used across
#'   PECBMS national partners): simulate each species' log-index
#'   `nSim` times by drawing from a normal distribution (mean = log(index),
#'   SD = that species' own log-scale standard error), take the geometric
#'   mean across species for each simulated draw (i.e. the mean of the
#'   simulated log-indices, exponentiated), then summarize the resulting
#'   `nSim` combined-index draws (percentile CI, matching this repo's other
#'   two combined-index functions' bootstrap CIs).
#'
#' Same point estimate as `computeCombinedIndex()`/`computeCombinedIndexAnalytical()`
#' (all three are the geometric mean of the same per-species indices) --
#' only the uncertainty differs across all three.
#'
#' @param indexMat Numeric matrix, species (rows) x years (columns), of
#'   per-species indices already baselined via `computeSpeciesIndex()`.
#' @param seIndexMat Numeric matrix, same dimensions/dimnames as `indexMat`,
#'   of each index value's standard error on the natural scale (see
#'   `computeIndexSE()`) -- converted internally to log scale via the
#'   delta-method approximation `SE(log(I)) ~= SE(I) / I`.
#' @param nSim Integer. Monte Carlo simulations per year (Soldaat et al. use
#'   1000).
#' @param seed Integer. For reproducible simulation.
#' @return data.frame with columns `year`, `SBI`, `se`, `lwr95`, `upr95`.
computeCombinedIndexMSI <- function(indexMat, seIndexMat, nSim = 1000, seed = 42) {
  years <- colnames(indexMat)
  set.seed(seed)

  rows <- lapply(years, function(yr) {
    idx <- indexMat[, yr]
    se <- seIndexMat[, yr]
    valid <- !is.na(idx) & !is.na(se) & idx > 0
    idx <- idx[valid]
    se <- se[valid]

    if (length(idx) == 0) {
      return(c(SBI = NA_real_, se = NA_real_, lwr95 = NA_real_, upr95 = NA_real_))
    }

    logIdx <- log(idx)
    logSE <- se / idx  # delta-method: SE(log(I)) ~= SE(I) / I

    simMeans <- vapply(seq_len(nSim), function(s) {
      draws <- stats::rnorm(length(logIdx), mean = logIdx, sd = logSE)
      mean(draws)
    }, numeric(1))

    msiDraws <- exp(simMeans)

    c(SBI = exp(mean(logIdx)), se = stats::sd(msiDraws),
      lwr95 = as.numeric(stats::quantile(msiDraws, 0.025, na.rm = TRUE)),
      upr95 = as.numeric(stats::quantile(msiDraws, 0.975, na.rm = TRUE)))
  })

  out <- do.call(rbind, rows)
  data.frame(year = as.integer(years), SBI = out[, "SBI"], se = out[, "se"],
             lwr95 = out[, "lwr95"], upr95 = out[, "upr95"], row.names = NULL)
}
