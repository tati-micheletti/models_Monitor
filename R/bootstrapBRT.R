#' Building blocks for the spatial-block bootstrap of the BRTs (uncertainty option "B")
#'
#' Plan and decisions: improvements.md item 15. Each replicate resamples spatial BLOCKS of the training
#' data with replacement, refits the BRT with the FIXED hyperparameters of the main fit (tree count,
#' learning rate, bag fraction, tree complexity -- as the Boreal Avian Modelling Project does), and
#' predicts. Everything derived (maps, changes, trends, area means) must be computed INSIDE each
#' replicate, then summarised across replicates; that is what carries spatial AND temporal uncertainty
#' into per-pixel trends.

#' Resampling-block id for every training row
#'
#' A square grid of `blockSizeM`, anchored at 0 (blocks only need to be spatial clusters of roughly the
#' spatial-CV block size; they do not have to coincide with the CV blocks, which are not saved for the
#' cluster tasks).
#'
#' @param spPa data.frame with columns `x`, `y` (EPSG:3035).
#' @param blockSizeM Numeric. Block side length (m).
#' @return Character vector, one block id per row.
bootstrapBlockIds <- function(spPa, blockSizeM) {
  paste(floor(spPa$x / blockSizeM), floor(spPa$y / blockSizeM), sep = "_")
}

#' Row indices of one block-bootstrap replicate
#'
#' Draws as many blocks as there are blocks, with replacement, and returns all rows of the drawn blocks
#' (a block drawn twice contributes its rows twice). Redraws (deterministically) if a draw lacks
#' presences or absences.
#'
#' @param blockIds Character vector from `bootstrapBlockIds()`.
#' @param occurrence 0/1 vector, same length.
#' @param seed Integer. Replicate seed (fixed per replicate so replicates can be ADDED later).
#' @param minPerClass Integer. Minimum presences and absences in a usable draw.
#' @return Integer vector of row indices.
bootstrapRowIndex <- function(blockIds, occurrence, seed, minPerClass = 10) {
  rowsByBlock <- split(seq_along(blockIds), blockIds)
  nBlocks <- length(rowsByBlock)
  for (attempt in 0:49) {
    set.seed(seed + attempt * 1000003L)
    idx <- unlist(rowsByBlock[sample.int(nBlocks, nBlocks, replace = TRUE)], use.names = FALSE)
    if (sum(occurrence[idx] == 1) >= minPerClass && sum(occurrence[idx] == 0) >= minPerClass) return(idx)
  }
  stop("bootstrapRowIndex(): no usable draw after 50 attempts (too few blocks/records?).")
}

#' Refit a BRT on a bootstrap sample with the main fit's FIXED hyperparameters
#'
#' @param spPa data.frame (model-ready table) with `occurrence` and the predictors.
#' @param predSel Character. Predictor columns.
#' @param brtM The main `gbm.step()` model (hyperparameters read from `$gbm.call`).
#' @param idx Integer row indices from `bootstrapRowIndex()`.
#' @return A `gbm` model (fitted with `keep.data = FALSE`, so small).
fitBootstrapBRT <- function(spPa, predSel, brtM, idx) {
  d <- spPa[idx, c("occurrence", predSel), drop = FALSE]
  gbm::gbm(formula = occurrence ~ ., distribution = "bernoulli", data = d,
           n.trees = brtM$gbm.call$best.trees, shrinkage = brtM$gbm.call$learning.rate,
           bag.fraction = brtM$gbm.call$bag.fraction, interaction.depth = brtM$gbm.call$tree.complexity,
           weights = rep(1, nrow(d)), verbose = FALSE, keep.data = FALSE)
}

#' Covariates of one year as a data.frame of the cells that can be predicted
#'
#' Same cell selection as `predictBRTToRaster()` (complete cases only; x/y from the raster itself), built
#' ONCE per year and scale and then reused by every replicate.
#'
#' @param covStack SpatRaster of one year's covariates.
#' @param predictors Character. Predictor columns (may include "x","y").
#' @return List: `df` (complete cells, predictor columns), `idx` (positions of those cells in the full
#'   raster), `template` (single-layer SpatRaster to write values into).
prepareCovariateCells <- function(covStack, predictors) {
  rasterPredictors <- setdiff(predictors, c("x", "y"))
  predRast <- covStack[[rasterPredictors]]
  predDf <- as.data.frame(predRast, xy = TRUE, na.rm = FALSE)
  complete <- stats::complete.cases(predDf[, predictors])
  list(df = predDf[complete, predictors, drop = FALSE], idx = which(complete), template = predRast[[1]])
}

#' Predict a BRT on prepared cells, optionally in parallel chunks
#'
#' @param model A `gbm` model.
#' @param cells Output of `prepareCovariateCells()`.
#' @param nTrees Integer. Trees to use.
#' @param nCores Integer. Cores for chunked prediction (forked; 1 on Windows).
#' @return Numeric vector of probabilities, one per complete cell.
predictBRTCells <- function(model, cells, nTrees, nCores = 1) {
  n <- nrow(cells$df)
  if (nCores <= 1 || .Platform$OS.type == "windows") {
    return(gbm::predict.gbm(model, cells$df, n.trees = nTrees, type = "response"))
  }
  chunks <- split(seq_len(n), cut(seq_len(n), nCores, labels = FALSE))
  unlist(parallel::mclapply(chunks, function(i)
    gbm::predict.gbm(model, cells$df[i, , drop = FALSE], n.trees = nTrees, type = "response"),
    mc.cores = nCores), use.names = FALSE)
}

#' Summarise replicates across the columns of a cells x B matrix
#'
#' @param m Numeric matrix, rows = cells, columns = replicates.
#' @param probs Numeric length 2. Lower/upper percentile (default 5th/95th = the 90% interval).
#' @return Matrix with columns mean, sd, lwr, upr, width.
summarizeReplicates <- function(m, probs = c(0.05, 0.95)) {
  q <- matrixStats::rowQuantiles(m, probs = probs, na.rm = TRUE)
  cbind(mean = rowMeans(m, na.rm = TRUE), sd = matrixStats::rowSds(m, na.rm = TRUE),
        lwr = q[, 1], upr = q[, 2], width = q[, 2] - q[, 1])
}

#' Share of replicates in which a change is negative (decrease) / positive (increase)
#'
#' @param change Numeric matrix, rows = cells, columns = replicates (late minus early, per replicate).
#' @return Matrix with columns shareDecrease, shareIncrease.
shareChangeSign <- function(change) {
  cbind(shareDecrease = rowMeans(change < 0, na.rm = TRUE),
        shareIncrease = rowMeans(change > 0, na.rm = TRUE))
}
