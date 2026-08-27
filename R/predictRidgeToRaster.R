#' Predict a fitted ridge model onto a stacked suitability raster
#'
#' Same reconstruction approach as `predictBRTToRaster()`: assigns
#' predicted values into a template layer's cells directly, so the
#' output is guaranteed grid-identical to `suitStack`.
#'
#' @param suitStack SpatRaster stack with the three `suitCols` layers.
#' @param suitCols Character vector of predictor column names, in the
#'   order the ridge model was trained on.
#' @param ridgeModel A `cv.glmnet()` model.
#' @param lambda Numeric. Penalty to predict at (e.g. `lambda.1se`).
#' @param thresh Numeric. Binary classification threshold.
#' @param outPath Character. Output file path.
#' @return Invisibly, `outPath`.
predictRidgeToRaster <- function(suitStack, suitCols, ridgeModel, lambda, thresh, outPath) {
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

  writeTwoLayerRaster(rMeta, rBinary, c("meta_prob", "binary"), outPath)
}
