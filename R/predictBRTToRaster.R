#' Predict a fitted BRT model onto a covariate raster stack
#'
#' Builds a prediction data.frame from `covStack`, predicts with
#' `brtModel`, and rebuilds a full-extent raster (NA where covariates
#' were incomplete) with `mean_prob` and a `binary` layer thresholded at
#' `thresh`. Writes via `writeTwoLayerRaster()` to work around a terra
#' stacking bug.
#'
#' Reconstructs onto `covStack`'s own grid by directly assigning into a
#' template layer's cell values (rather than gridding predicted points
#' back via coordinates), so the output raster is guaranteed
#' cell-for-cell identical in extent/resolution to the input -- this
#' matters downstream since the meta-model resamples every scale's
#' prediction onto a shared reference grid.
#'
#' @param covStack SpatRaster stack containing at least the columns in `predictors`.
#' @param predictors Character vector of predictor column names.
#' @param brtModel A fitted `gbm.step()` model.
#' @param thresh Numeric. Binary classification threshold.
#' @param outPath Character. Output file path.
#' @return Invisibly, `outPath`.
predictBRTToRaster <- function(covStack, predictors, brtModel, thresh, outPath) {
  predRast <- covStack[[predictors]]
  predDf <- as.data.frame(predRast, xy = TRUE, na.rm = FALSE)

  completeIdx <- stats::complete.cases(predDf[, predictors])
  predDfComplete <- predDf[completeIdx, ]

  predVals <- gbm::predict.gbm(brtModel, predDfComplete[, predictors],
                                n.trees = brtModel$gbm.call$best.trees, type = "response")

  predFull <- rep(NA_real_, nrow(predDf))
  predFull[completeIdx] <- predVals

  rPred <- predRast[[1]]
  terra::values(rPred) <- predFull
  names(rPred) <- "mean_prob"

  rBinary <- rPred >= thresh
  names(rBinary) <- "binary"

  writeTwoLayerRaster(rPred, rBinary, c("mean_prob", "binary"), outPath)
}
