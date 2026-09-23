#' Combine per-species indices into one multi-species index (Zbinden et al. 2005 method)
#'
#' `SBI_t = ( prod_i I_i,t )^(1/n) = exp(mean(log(I_i,t)))` -- the geometric
#' mean across species at each year, computed only over species with a
#' non-NA index that year (so a species missing a year just drops out,
#' matching "combined indices can be calculated for any species set").
#' Geometric (not arithmetic) mean is deliberate: it's what makes a
#' doubling in one species exactly offset a halving in another, since
#' these are ratio/percent-of-baseline values, not counts.
#'
#' 95% confidence intervals follow the same species-resampling logic as
#' Freeman, Baillie & Gregory (2001, BTO Research Report 251), the
#' statistical basis Zbinden et al. cite for their own CIs: resample which
#' species are included (with replacement), recompute the whole combined
#' series, repeat `nBoot` times, take the 2.5th/97.5th percentile per year.
#'
#' @param indexMat Numeric matrix, species (rows) x years (columns), of
#'   per-species indices already baselined via `computeSpeciesIndex()`
#'   (column names = years as character).
#' @param nBoot Integer. Bootstrap replicates for the CI. Set to 0 to skip
#'   CI computation (faster, e.g. for interactive exploration).
#' @param seed Integer. For reproducible bootstrap resampling.
#' @return data.frame with columns `year`, `SBI`, `lwr95`, `upr95` (the
#'   last two are `NA` if `nBoot == 0`), one row per year in `indexMat`.
computeCombinedIndex <- function(indexMat, nBoot = 999, seed = 42) {
  years <- colnames(indexMat)
  nSp <- nrow(indexMat)

  geomMeanCol <- function(mat, col) {
    v <- mat[, col]
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NA_real_)
    exp(mean(log(v)))
  }

  sbi <- vapply(years, function(yr) geomMeanCol(indexMat, yr), numeric(1))

  if (nBoot > 0) {
    set.seed(seed)
    boot <- replicate(nBoot, {
      spSample <- sample.int(nSp, nSp, replace = TRUE)
      vapply(years, function(yr) geomMeanCol(indexMat[spSample, , drop = FALSE], yr), numeric(1))
    })
    # boot is years x nBoot
    lwr <- apply(boot, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
    upr <- apply(boot, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
  } else {
    lwr <- rep(NA_real_, length(years))
    upr <- rep(NA_real_, length(years))
  }

  data.frame(year = as.integer(years), SBI = sbi, lwr95 = lwr, upr95 = upr,
             row.names = NULL)
}
