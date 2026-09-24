#' Combine per-species indices via the DDA/PECBMS analytical delta-method CI
#'
#' Alternative to `computeCombinedIndex()`'s bootstrap CI, not a replacement
#' for it: same geometric-mean point estimate (confirmed identical
#' methodology in Gregory et al. 2005, "Developing indicators for European
#' birds", Phil. Trans. R. Soc. B 360:269-288 -- the paper describing how
#' PECBMS/DDA's national TRIM-based species indices are combined into
#' multi-species indicators), but the *uncertainty* is computed analytically
#' from each species' own year-specific variance via a delta-method
#' (Taylor-linearization) approximation, Appendix A of that paper (credited
#' there to J. Pannekoek of Statistics Netherlands, the TRIM author):
#'   I = exp( (1/T) * sum(log(I_t)) )                         (their eq. A1)
#'   var(I) ~= (I^2 / T^2) * sum( var(I_t) / I_t^2 )           (their eq. A3)
#' where T is the number of species contributing that year.
#'
#' Two real differences from PECBMS/DDA's own use of this formula, both
#' consequences of our data being a complete predicted raster surface for
#' every year rather than raw site counts with gaps (which is what TRIM's
#' own modelling solves): (1) `var(I_t)` here comes from the spatial
#' standard error of each year's mean occurrence probability
#' (`extractMetaProbSD()`), not from a TRIM count model's sampling variance
#' -- a genuinely different source of uncertainty (model-surface spatial
#' variability, not survey-design sampling error); (2) baseline-year
#' variance is treated as fixed/known (not itself propagated), a
#' simplification -- true uncertainty is very slightly understated as a
#' result.
#'
#' This CI is a symmetric normal-approximation (`I +/- 1.96*SE`), unlike
#' `computeCombinedIndex()`'s percentile-based bootstrap CI (which can be
#' asymmetric). Report both if there's any doubt about which is more
#' appropriate -- they answer subtly different questions (this one:
#' uncertainty from each species' own spatial variability; the bootstrap:
#' uncertainty from which particular species happen to be in the set).
#'
#' @param indexMat Numeric matrix, species (rows) x years (columns), of
#'   per-species indices already baselined via `computeSpeciesIndex()`.
#' @param seIndexMat Numeric matrix, same dimensions/dimnames as `indexMat`,
#'   of each index value's standard error (see `computeIndexSE()`).
#' @return data.frame with columns `year`, `SBI`, `se`, `lwr95`, `upr95`.
computeCombinedIndexAnalytical <- function(indexMat, seIndexMat) {
  years <- colnames(indexMat)

  rows <- lapply(years, function(yr) {
    idx <- indexMat[, yr]
    se <- seIndexMat[, yr]
    valid <- !is.na(idx) & !is.na(se) & idx > 0
    idx <- idx[valid]
    se <- se[valid]
    Tn <- length(idx)

    if (Tn == 0) {
      return(c(SBI = NA_real_, se = NA_real_, lwr95 = NA_real_, upr95 = NA_real_))
    }

    geomMean <- exp(mean(log(idx)))
    varI <- (geomMean^2 / Tn^2) * sum((se^2) / (idx^2))
    seI <- sqrt(varI)

    c(SBI = geomMean, se = seI, lwr95 = geomMean - 1.96 * seI, upr95 = geomMean + 1.96 * seI)
  })

  out <- do.call(rbind, rows)
  data.frame(year = as.integer(years), SBI = out[, "SBI"], se = out[, "se"],
             lwr95 = out[, "lwr95"], upr95 = out[, "upr95"], row.names = NULL)
}
