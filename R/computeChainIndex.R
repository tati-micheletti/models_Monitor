#' Combine per-species indices via the Living Planet Index chain method
#'
#' A genuinely different combination strategy from `computeCombinedIndex()`
#' (Zbinden et al. 2005/SBI-style), not just a different CI: instead of
#' taking the geometric mean of each year's *baseline-normalized level*
#' across species, it takes the geometric mean of each year-to-year
#' *log-ratio* (`lambda_t = log(I_t / I_{t-1})`) across species, then
#' cumulates ("chains") those combined log-ratios into an index series
#' (`Index_1 = 100`, `Index_t = Index_{t-1} * exp(mean_lambda_{t-1})`) --
#' following the Living Planet Index's log-linear chain method (McRae,
#' Deinet & Freeman 2017, PLOS ONE 12(1):e0169156; Collen et al. 2009,
#' Conservation Biology 23(2):317-327).
#'
#' Because a log-ratio between two years is unaffected by whatever a
#' species' baseline year was (`log((100*r_t/r_base)/(100*r_{t-1}/r_base))
#' = log(r_t/r_{t-1})`, `r_base` cancels), this method is mathematically
#' independent of baseline-year choice -- it always starts at 100 wherever
#' the series begins, sidestepping the "why 2005 and not 2022" question
#' entirely for this variant (though `computeSpeciesIndex()`'s baseline
#' choice still matters for the *other* (Zbinden-style) combined index and
#' for every per-species index reported alongside it).
#'
#' Real LPI implementations additionally: (1) fall back to this chain
#' method only for short/gappy time series and use a per-population GAM
#' otherwise (we have complete, gap-free series for all species, so that
#' distinction doesn't apply here); (2) diversity-weight the aggregation
#' across taxonomic/biogeographic groups when combining very different
#' species sets (not obviously relevant with the current roster of
#' similarly-surveyed farmland species -- see `sharedSpecies` in
#' sharedConfig.R for the count in force for any given run -- so left
#' unweighted here, matching Zbinden et al.'s own unweighted approach).
#'
#' @param indexMat Numeric matrix, species (rows) x years (columns), of
#'   per-species indices already baselined via `computeSpeciesIndex()`
#'   (column names = years as character, in chronological order).
#' @param nBoot Integer. Bootstrap replicates for the CI (same
#'   species-resampling logic as `computeCombinedIndex()`). 0 to skip.
#' @param seed Integer. For reproducible bootstrap resampling.
#' @return data.frame with columns `year`, `chainIndex`, `lwr95`, `upr95`.
computeChainIndex <- function(indexMat, nBoot = 999, seed = 42) {
  years <- colnames(indexMat)
  nYears <- length(years)
  nSp <- nrow(indexMat)

  lambdaMat <- matrix(NA_real_, nrow = nSp, ncol = nYears - 1,
                       dimnames = list(rownames(indexMat), NULL))
  for (i in seq_len(nYears - 1)) {
    lambdaMat[, i] <- log(indexMat[, i + 1] / indexMat[, i])
  }

  chainFromLambda <- function(lam) {
    meanLambda <- apply(lam, 2, mean, na.rm = TRUE)
    idx <- numeric(nYears)
    idx[1] <- 100
    for (i in seq_len(nYears - 1)) {
      idx[i + 1] <- idx[i] * exp(meanLambda[i])
    }
    idx
  }

  chainIdx <- chainFromLambda(lambdaMat)

  if (nBoot > 0) {
    set.seed(seed)
    boot <- replicate(nBoot, {
      spSample <- sample.int(nSp, nSp, replace = TRUE)
      chainFromLambda(lambdaMat[spSample, , drop = FALSE])
    })
    # boot is years x nBoot
    lwr <- apply(boot, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
    upr <- apply(boot, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
  } else {
    lwr <- rep(NA_real_, nYears)
    upr <- rep(NA_real_, nYears)
  }

  data.frame(year = as.integer(years), chainIndex = chainIdx, lwr95 = lwr, upr95 = upr,
             row.names = NULL)
}
