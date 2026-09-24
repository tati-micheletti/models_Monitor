#' Bootstrap the ridge meta-model's fit to quantify trend uncertainty
#'
#' Refits the ridge regression `nBoot` times on bootstrap resamples (with
#' replacement) of the SAME training rows (`X`, `y`) `metaModel()` uses for
#' its single point-estimate fit, then re-applies each bootstrap model's
#' coefficients to the SAME already-computed per-cell suitability values for
#' each prediction year (`newXByYear`), returning the bootstrap distribution
#' of the area-mean predicted suitability per year.
#'
#' This quantifies genuine model-fitting (sampling) uncertainty in HOW THE
#' THREE SCALES ARE COMBINED -- distinct from, and additional to, the
#' spatial-averaging SE already used for area-mean precision (which treats
#' the fitted model as exactly known and only asks "how precisely do we know
#' the mean of a fixed surface"). It does NOT capture uncertainty in the
#' three underlying BRT scale models themselves -- refitting each BRT
#' ensemble `nBoot` times would be orders of magnitude more expensive
#' (each BRT fit is itself tens of minutes to hours; see
#' `improvements.md`) and is left as a candidate follow-up once real EVE
#' cluster capacity makes that many parallel refits feasible.
#'
#' @param X Numeric matrix, ridge training predictors (same as `metaModel()`'s
#'   `X`: columns `climate_mean_prob`/`landscape_mean_prob`/`habitat_mean_prob`).
#' @param y Numeric vector, ridge training response (same as `metaModel()`'s `y`).
#' @param newXByYear Named list (names = years, as character) of numeric
#'   matrices -- full-grid suitability values (non-NA cells only) for that
#'   year, same column order as `X`. This is the same matrix
#'   `predictRidgeToRaster()` builds internally as `xPred`.
#' @param nBoot Integer. Number of bootstrap replicates. Default 200.
#' @param seed Integer, for reproducibility.
#' @return Data frame with one row per year: `year`, `bootMean`, `bootSE`,
#'   `lwr95`/`upr95` (2.5th/97.5th percentiles of the bootstrap distribution
#'   of the area-mean predicted suitability), and `nBootOK` (how many of the
#'   `nBoot` resamples produced a usable fit -- resamples that happen to draw
#'   only one class are skipped, not treated as zero).
bootstrapMetaModelTrend <- function(X, y, newXByYear, nBoot = 200, seed = 42) {
  set.seed(seed)
  n <- nrow(X)
  years <- names(newXByYear)

  bootAreaMeans <- matrix(NA_real_, nrow = nBoot, ncol = length(years),
                           dimnames = list(NULL, years))

  for (b in seq_len(nBoot)) {
    idx <- sample.int(n, n, replace = TRUE)
    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]

    if (length(unique(yb)) < 2) next

    fitB <- tryCatch(
      glmnet::cv.glmnet(x = Xb, y = yb, family = "binomial", alpha = 0,
                         nfolds = 5, standardize = TRUE),
      error = function(e) NULL)
    if (is.null(fitB)) next

    for (yr in years) {
      predB <- as.vector(stats::predict(fitB, newx = newXByYear[[yr]],
                                          s = fitB$lambda.1se, type = "response"))
      bootAreaMeans[b, yr] <- mean(predB, na.rm = TRUE)
    }
  }

  data.frame(
    year = as.integer(years),
    bootMean = apply(bootAreaMeans, 2, mean, na.rm = TRUE),
    bootSE = apply(bootAreaMeans, 2, stats::sd, na.rm = TRUE),
    lwr95 = apply(bootAreaMeans, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    upr95 = apply(bootAreaMeans, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
    nBootOK = colSums(!is.na(bootAreaMeans)),
    row.names = NULL
  )
}
