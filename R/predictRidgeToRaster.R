#' Predict a fitted ridge model onto a stacked suitability raster
#'
#' Same reconstruction approach as `predictBRTToRaster()`: assigns
#' predicted values into a template layer's cells directly, so the
#' output is guaranteed grid-identical to `suitStack`.
#'
#' Returns the raster rather than writing it to a fixed path (unlike its
#' pre-2026-09-30 version) so it can be wrapped in `reproducible::Cache()`
#' -- see `metaModel()`, which calls this then writes the (possibly-
#' cached) result to that run's own known output path. Mirrors
#' `predictBRTToRaster()`'s own 2026-09-28 conversion for the same reason.
#'
#' @param suitStack SpatRaster stack with the three `suitCols` layers.
#' @param suitCols Character vector of predictor column names, in the
#'   order the ridge model was trained on.
#' @param ridgeModel A `cv.glmnet()` model.
#' @param lambda Numeric. Penalty to predict at (e.g. `lambda.1se`).
#' @param thresh Numeric. Binary classification threshold.
#' @return SpatRaster, two layers (`meta_prob`, `binary`).
predictRidgeToRaster <- function(suitStack, suitCols, ridgeModel, lambda, thresh) {
  predDf <- as.data.frame(suitStack, xy = TRUE, na.rm = FALSE)
  completeIdx <- stats::complete.cases(predDf[, suitCols])
  predComplete <- predDf[completeIdx, ]

  xPred <- as.matrix(predComplete[, suitCols])
  predVals <- as.vector(stats::predict(ridgeModel, newx = xPred, s = lambda, type = "response"))

  predFull <- rep(NA_real_, nrow(predDf))
  predFull[completeIdx] <- predVals

  rMeta <- suitStack[[1]]
  terra::values(rMeta) <- predFull
  names(rMeta) <- "meta_prob"

  rBinary <- rMeta >= thresh
  names(rBinary) <- "binary"

  combineTwoLayerRaster(rMeta, rBinary, c("meta_prob", "binary"))
}
