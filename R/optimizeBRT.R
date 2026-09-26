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
#' @param minLR Numeric. If halving would take the learning rate below this,
#'   stop and error instead of continuing -- default 1e-6. Below this, the
#'   problem is essentially always that `gbm.step()` cannot fit ANYTHING
#'   (e.g. a response with zero variance -- no presences or no absences at
#'   all), not a learning-rate tuning problem: halving a fixed value toward
#'   zero eventually underflows to exactly `0.0` in floating point, and
#'   `0 / 2 == 0`, so without this floor the search would loop forever,
#'   repeatedly calling `gbm.step(learning.rate = 0)` and never converging
#'   or erroring on its own (this happened for real -- Anthus pratensis at
#'   Europe scale, 2026-09-26, "total mean deviance = NaN" repeating with no
#'   end, root cause was 0 presences from a stale EBBA2 extract).
#' @param maxLR Numeric. Same idea for doubling, default 1. Symmetric safety
#'   net, not yet observed to trigger in practice.
#' @param maxIter Integer. Absolute cap on optimization iterations regardless
#'   of `minLR`/`maxLR`, default 50 -- a second, independent safety net.
#' @return A fitted `gbm.step()` model object.
optimizeBRT <- function(data, gbmX, gbmY, initialLR, treeComplexity = 2,
                         bagFraction = 0.75, minTrees = 1000, maxTrees = 5000,
                         minLR = 1e-6, maxLR = 1, maxIter = 50) {

  LR <- initialLR
  optimizing <- TRUE
  brtM <- NULL
  iter <- 0L

  while (optimizing) {
    iter <- iter + 1L
    if (iter > maxIter) {
      stop("optimizeBRT() did not converge within ", maxIter, " iterations ",
           "(last LR = ", LR, ") -- gbm.step() is very likely failing for a ",
           "reason no amount of learning-rate adjustment will fix (e.g. the ",
           "response has zero variance -- check for 0 presences or 0 ",
           "absences in this species' data).")
    }

    brtM <- try(
      dismo::gbm.step(data = data, gbm.x = gbmX, gbm.y = gbmY, family = "bernoulli",
                       tree.complexity = treeComplexity, bag.fraction = bagFraction,
                       learning.rate = LR, verbose = FALSE, plot.main = FALSE),
      silent = TRUE)

    if (inherits(brtM, "try-error") || is.null(brtM)) {
      LR <- LR / 2
      if (LR < minLR) {
        stop("optimizeBRT() gave up: learning rate fell below minLR (", minLR, ") ",
             "after ", iter, " failed attempt(s) -- gbm.step() cannot fit this data ",
             "regardless of learning rate. Most likely cause: the response has zero ",
             "variance (0 presences or 0 absences). Check this species' occurrence ",
             "data before re-running.")
      }
      message("BRT failed -- halving LR to ", LR)
    } else if (brtM$gbm.call$best.trees < minTrees) {
      LR <- LR / 2
      if (LR < minLR) {
        stop("optimizeBRT() gave up: learning rate fell below minLR (", minLR, ") ",
             "after ", iter, " attempt(s), still producing too few trees. Check this ",
             "species' occurrence data (sample size, class balance) before re-running.")
      }
      message("Too few trees (", brtM$gbm.call$best.trees, ") -- halving LR to ", LR)
    } else if (brtM$gbm.call$best.trees > maxTrees) {
      LR <- LR * 2
      if (LR > maxLR) {
        stop("optimizeBRT() gave up: learning rate exceeded maxLR (", maxLR, ") ",
             "after ", iter, " attempt(s), still producing too many trees.")
      }
      message("Too many trees (", brtM$gbm.call$best.trees, ") -- doubling LR to ", LR)
    } else {
      optimizing <- FALSE
      message("BRT converged: ", brtM$gbm.call$best.trees, " trees, LR = ", LR)
    }
  }

  brtM
}
