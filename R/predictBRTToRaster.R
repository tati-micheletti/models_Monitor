#' Predict a fitted BRT model onto a covariate raster stack
#'
#' Builds a prediction data.frame from `covStack`, predicts with
#' `brtModel`, and rebuilds a full-extent raster (NA where covariates
#' were incomplete) with `mean_prob` and a `binary` layer thresholded at
#' `thresh`. Combined via `combineTwoLayerRaster()` to work around a terra
#' stacking bug.
#'
#' Reconstructs onto `covStack`'s own grid by directly assigning into a
#' template layer's cell values (rather than gridding predicted points
#' back via coordinates), so the output raster is guaranteed
#' cell-for-cell identical in extent/resolution to the input -- this
#' matters downstream since the meta-model resamples every scale's
#' prediction onto a shared reference grid.
#'
#' Returns the raster rather than writing it to a fixed path (unlike its
#' pre-2026-09-28 version) so it can be wrapped in `reproducible::Cache()`
#' -- see `modelGerLandscape()`/`GerHabitat()`/`modelEurope()`, which each
#' call this then write the (possibly-cached) result to that run's own
#' known output path.
#'
#' @param covStack SpatRaster stack containing at least the non-coordinate
#'   columns in `predictors` -- `x`/`y` (if present in `predictors`) are
#'   never expected as actual layers; they're derived from `covStack`'s own
#'   cell coordinates instead (see DECISIONS.md, 2026-09-26).
#' @param predictors Character vector of predictor column names.
#' @param brtModel A fitted `gbm.step()` model.
#' @param thresh Numeric. Binary classification threshold.
#' @return SpatRaster, two layers (`mean_prob`, `binary`).
predictBRTToRaster <- function(covStack, predictors, brtModel, thresh) {
  rasterPredictors <- setdiff(predictors, c("x", "y"))
  predRast <- covStack[[rasterPredictors]]
  # xy = TRUE supplies "x"/"y" columns (this raster's own cell coordinates,
  # same CRS/grid as everything else) for free -- exactly matching how the
  # training data's own x/y columns were extracted (terra::extract() from
  # the same covariate stack), so no separate synthetic layer is needed.
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

  combineTwoLayerRaster(rPred, rBinary, c("mean_prob", "binary"))
}
