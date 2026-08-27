#' Fit a boosted regression tree, searching for a learning rate that
#' yields 1000-5000 trees
#'
#' Following Wiedenroth et al.: halves the learning rate if the model
#' fails to fit or produces too few trees, doubles it if it produces too
#' many, until `dismo::gbm.step()` converges within the target tree
#' range.
#'
#' @param data data.frame with the response and predictor columns.
#' @param gbmX Character vector of predictor column names.
#' @param gbmY Character. Response column name (must be 0/1).
#' @param initialLR Numeric. Starting learning rate.
#' @param treeComplexity Integer. `gbm.step()` tree.complexity, default 2.
#' @param bagFraction Numeric. `gbm.step()` bag.fraction, default 0.75.
#' @param minTrees Integer. Minimum acceptable best.trees, default 1000.
#' @param maxTrees Integer. Maximum acceptable best.trees, default 5000.
#' @return A fitted `gbm.step()` model object.
optimizeBRT <- function(data, gbmX, gbmY, initialLR, treeComplexity = 2,
                         bagFraction = 0.75, minTrees = 1000, maxTrees = 5000) {

  LR <- initialLR
  optimizing <- TRUE
  brtM <- NULL

  while (optimizing) {
    brtM <- try(
      dismo::gbm.step(data = data, gbm.x = gbmX, gbm.y = gbmY, family = "bernoulli",
                       tree.complexity = treeComplexity, bag.fraction = bagFraction,
                       learning.rate = LR, verbose = FALSE, plot.main = FALSE),
      silent = TRUE)

    if (inherits(brtM, "try-error") || is.null(brtM)) {
      LR <- LR / 2
      message("BRT failed -- halving LR to ", LR)
    } else if (brtM$gbm.call$best.trees < minTrees) {
      LR <- LR / 2
      message("Too few trees (", brtM$gbm.call$best.trees, ") -- halving LR to ", LR)
    } else if (brtM$gbm.call$best.trees > maxTrees) {
      LR <- LR * 2
      message("Too many trees (", brtM$gbm.call$best.trees, ") -- doubling LR to ", LR)
    } else {
      optimizing <- FALSE
      message("BRT converged: ", brtM$gbm.call$best.trees, " trees, LR = ", LR)
    }
  }

  brtM
}
