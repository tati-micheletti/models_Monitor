#' Leave-one-out variable importance for the meta-model
#'
#' Importance of a predictor = drop in `dev.ratio` when that predictor
#' is removed (its unique contribution), floored at 0 and normalised
#' across the three predictors so they sum to 1 -- scale-invariant and
#' comparable across species, unlike raw glmnet coefficients. Follows
#' Wiedenroth et al. NOTE: unlike Wiedenroth's version, no
#' `lower.limits = 0` is applied in `fitDevRatio()`, matching how the
#' actual ridge meta-model is fit (unconstrained), so this stays
#' internally consistent with it.
#'
#' @param trainDf data.frame with `occurrence`, `climate_mean_prob`,
#'   `landscape_mean_prob`, and `habitat_mean_prob` columns.
#' @param suitCols Character vector `c("climate_mean_prob",
#'   "landscape_mean_prob", "habitat_mean_prob")`, in that order.
#' @param lambda Numeric. Ridge penalty (the meta-model's chosen lambda).
#' @return One-row data.frame with `dev_ratio_full`, `imp_climate`,
#'   `imp_landscape`, `imp_habitat` (all NA if total importance is zero).
computeVariableImportance <- function(trainDf, suitCols, lambda) {
  X <- as.matrix(trainDf[, suitCols])
  y <- trainDf$occurrence

  fullDev <- fitDevRatio(X, y, lambda)

  dropOne <- function(dropCol) {
    keep <- setdiff(suitCols, dropCol)
    fitDevRatio(X[, keep, drop = FALSE], y, lambda)
  }

  impClimate <- max(fullDev - dropOne("climate_mean_prob"), 0)
  impLandscape <- max(fullDev - dropOne("landscape_mean_prob"), 0)
  impHabitat <- max(fullDev - dropOne("habitat_mean_prob"), 0)

  total <- sum(impClimate, impLandscape, impHabitat)

  if (total == 0) {
    return(data.frame(dev_ratio_full = fullDev,
                       imp_climate = NA, imp_landscape = NA, imp_habitat = NA))
  }

  data.frame(dev_ratio_full = fullDev,
             imp_climate = impClimate / total,
             imp_landscape = impLandscape / total,
             imp_habitat = impHabitat / total)
}
