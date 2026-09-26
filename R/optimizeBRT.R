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
#'   give up and return NULL (with a warning) instead of continuing --
#'   default 1e-6. Below this, the problem is essentially always that
#'   `gbm.step()` cannot fit ANYTHING (e.g. a response with zero variance --
#'   no presences or no absences at all), not a learning-rate tuning
#'   problem: halving a fixed value toward zero eventually underflows to
#'   exactly `0.0` in floating point, and `0 / 2 == 0`, so without this
#'   floor the search would loop forever, repeatedly calling
#'   `gbm.step(learning.rate = 0)` and never converging or erroring on its
#'   own (this happened for real -- Anthus pratensis at Europe scale,
#'   2026-09-26, "total mean deviance = NaN" repeating with no end, root
#'   cause was 0 presences from a stale EBBA2 extract). A NULL return here
#'   means one species failing shouldn't take the whole run down -- see
#'   `modelEurope()`/`modelGerHabitat()`/`modelGerLandscape()`, which all
#'   skip to the next species (with a warning) when this returns NULL,
#'   rather than erroring.
#' @param maxLR Numeric. Same idea for doubling, default 1. Symmetric safety
#'   net, not yet observed to trigger in practice.
#' @param maxIter Integer. Absolute cap on optimization iterations regardless
#'   of `minLR`/`maxLR`, default 50 -- a second, independent safety net.
#' @return A fitted `gbm.step()` model object, or NULL if it gave up (see
#'   `minLR`/`maxLR`/`maxIter` above) -- always check for NULL before using
#'   the result.
optimizeBRT <- function(data, gbmX, gbmY, initialLR, treeComplexity = 2,
                         bagFraction = 0.75, minTrees = 1000, maxTrees = 5000,
                         minLR = 1e-6, maxLR = 1, maxIter = 50) {

  LR <- initialLR
  optimizing <- TRUE
  brtM <- NULL
  iter <- 0L

  giveUp <- function(reason) {
    warning("optimizeBRT() gave up: ", reason, " Most likely cause: the response ",
            "has zero variance (0 presences or 0 absences), or too small a sample ",
            "size to fit at all -- check this species' occurrence data. Returning ",
            "NULL; the caller should skip this species rather than treat it as fit.")
    NULL
  }

  while (optimizing) {
    iter <- iter + 1L
    if (iter > maxIter) {
      return(giveUp(paste0("did not converge within ", maxIter, " iterations (last LR = ", LR, ").")))
    }

    brtM <- try(
      dismo::gbm.step(data = data, gbm.x = gbmX, gbm.y = gbmY, family = "bernoulli",
                       tree.complexity = treeComplexity, bag.fraction = bagFraction,
                       learning.rate = LR, verbose = FALSE, plot.main = FALSE),
      silent = TRUE)

    if (inherits(brtM, "try-error") || is.null(brtM)) {
      LR <- LR / 2
      if (LR < minLR) {
        return(giveUp(paste0("learning rate fell below minLR (", minLR, ") after ",
                              iter, " failed attempt(s).")))
      }
      message("BRT failed -- halving LR to ", LR)
    } else if (brtM$gbm.call$best.trees < minTrees) {
      LR <- LR / 2
      if (LR < minLR) {
        return(giveUp(paste0("learning rate fell below minLR (", minLR, ") after ",
                              iter, " attempt(s), still producing too few trees.")))
      }
      message("Too few trees (", brtM$gbm.call$best.trees, ") -- halving LR to ", LR)
    } else if (brtM$gbm.call$best.trees > maxTrees) {
      LR <- LR * 2
      if (LR > maxLR) {
        return(giveUp(paste0("learning rate exceeded maxLR (", maxLR, ") after ",
                              iter, " attempt(s), still producing too many trees.")))
      }
      message("Too many trees (", brtM$gbm.call$best.trees, ") -- doubling LR to ", LR)
    } else {
      optimizing <- FALSE
      message("BRT converged: ", brtM$gbm.call$best.trees, " trees, LR = ", LR)
    }
  }

  brtM
}
