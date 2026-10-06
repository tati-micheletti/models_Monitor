#' Predict a fitted BRT model onto a covariate raster stack
#'
#' Predicts with `brtModel` for every cell of `covStack` whose predictors are
#' all available, and rebuilds a full-extent raster (NA where covariates
#' were incomplete) with `mean_prob` and a `binary` layer thresholded at
#' `thresh`. Combined via `combineTwoLayerRaster()` to work around a terra
#' stacking bug.
#'
#' The cells are processed in blocks of whole rows (`chunkCells` cells at a
#' time), NOT as one giant data.frame: at 200 m the whole country is 14.85M
#' cells x ~20 layers, and building that data.frame (plus its complete-case
#' copy and gbm's own matrix) needed ~10 GB on top of the 2.5 GB covariate
#' stack -- the EVE habitat array was killed for running out of memory.
#' Blockwise prediction gives IDENTICAL values (a prediction depends only on
#' its own row) with a small fraction of the memory.
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
#' @param chunkCells Integer. Approximate number of cells predicted at a time (whole rows).
#' @return SpatRaster, two layers (`mean_prob`, `binary`).
predictBRTToRaster <- function(covStack, predictors, brtModel, thresh, chunkCells = 2e6) {
  rasterPredictors <- setdiff(predictors, c("x", "y"))
  predRast <- covStack[[rasterPredictors]]
  nr <- terra::nrow(predRast); nc <- terra::ncol(predRast)
  rowsPerChunk <- max(1L, as.integer(floor(chunkCells / nc)))
  # A file-backed stack (e.g. the habitat covariate cache) must be opened for the blockwise reads.
  terra::readStart(predRast)
  on.exit(terra::readStop(predRast), add = TRUE)

  predFull <- rep(NA_real_, nr * nc)
  for (r1 in seq(1L, nr, by = rowsPerChunk)) {
    nRows <- min(rowsPerChunk, nr - r1 + 1L)
    cellIds <- ((r1 - 1L) * nc + 1L):((r1 - 1L + nRows) * nc)
    # "x"/"y" are this raster's own cell coordinates (same CRS/grid as everything else), exactly matching
    # how the training data's own x/y columns were extracted (terra::extract() from the same covariate
    # stack), so no separate synthetic layer is needed.
    blockDf <- cbind(as.data.frame(terra::xyFromCell(predRast, cellIds)),
                     terra::readValues(predRast, row = r1, nrows = nRows, dataframe = TRUE))
    completeIdx <- stats::complete.cases(blockDf[, predictors, drop = FALSE])
    if (any(completeIdx)) {
      predFull[cellIds[completeIdx]] <- gbm::predict.gbm(
        brtModel, blockDf[completeIdx, predictors, drop = FALSE],
        n.trees = brtModel$gbm.call$best.trees, type = "response")
    }
  }

  rPred <- predRast[[1]]
  terra::values(rPred) <- predFull
  names(rPred) <- "mean_prob"

  rBinary <- rPred >= thresh
  names(rBinary) <- "binary"

  combineTwoLayerRaster(rPred, rBinary, c("mean_prob", "binary"))
}
