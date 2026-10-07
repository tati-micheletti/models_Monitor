#' Predict ANY model onto a covariate raster stack (blockwise)
#'
#' The algorithm-agnostic twin of `predictBRTToRaster()` (which stays untouched, so the BRT workflow and its caches do not change):
#' same blockwise reading (whole rows, `chunkCells` cells at a time), same complete-case rule, same output (`mean_prob` + `binary`),
#' but the prediction itself is done by `predictFun`, so GLM, GAM and random-forest fits (`algoFit()`) use the identical machinery.
#'
#' @param covStack SpatRaster with (at least) the non-coordinate predictors; "x"/"y" predictors come from the cell coordinates.
#' @param predictors Character. Predictor column names the model needs.
#' @param predictFun Function(data.frame) -> numeric vector of probabilities, one per row.
#' @param thresh Numeric. Threshold of the `binary` layer.
#' @param chunkCells Integer. Approximate number of cells per block.
#' @return SpatRaster, two layers (`mean_prob`, `binary`).
predictModelToRaster <- function(covStack, predictors, predictFun, thresh, chunkCells = 2e6) {
  rasterPredictors <- setdiff(predictors, c("x", "y"))
  predRast <- covStack[[rasterPredictors]]
  nr <- terra::nrow(predRast); nc <- terra::ncol(predRast)
  rowsPerChunk <- max(1L, as.integer(floor(chunkCells / nc)))
  terra::readStart(predRast)
  on.exit(terra::readStop(predRast), add = TRUE)

  predFull <- rep(NA_real_, nr * nc)
  for (r1 in seq(1L, nr, by = rowsPerChunk)) {
    nRows <- min(rowsPerChunk, nr - r1 + 1L)
    cellIds <- ((r1 - 1L) * nc + 1L):((r1 - 1L + nRows) * nc)
    blockDf <- cbind(as.data.frame(terra::xyFromCell(predRast, cellIds)),
                     terra::readValues(predRast, row = r1, nrows = nRows, dataframe = TRUE))
    completeIdx <- stats::complete.cases(blockDf[, predictors, drop = FALSE])
    if (any(completeIdx)) predFull[cellIds[completeIdx]] <- predictFun(blockDf[completeIdx, predictors, drop = FALSE])
  }

  rPred <- predRast[[1]]
  terra::values(rPred) <- predFull
  names(rPred) <- "mean_prob"
  rBinary <- rPred >= thresh
  names(rBinary) <- "binary"
  combineTwoLayerRaster(rPred, rBinary, c("mean_prob", "binary"))
}
