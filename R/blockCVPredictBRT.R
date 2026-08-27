#' Block cross-validated predictions from a fitted BRT model's hyperparameters
#'
#' Refits a `gbm()` per fold (using the exact tree count, shrinkage, bag
#' fraction, and tree complexity `brtModel` converged on) on the
#' out-of-fold rows, and predicts on the held-out fold. Uses the
#' `foldID` vector directly -- sourced from inputs_Monitor's
#' `sim$inputsData$<scale>[[species]]$data$foldID` column, which is
#' populated identically whether spatial blocking was real or mimicked
#' (see inputs_Monitor's mimicSpatialBlocks()) -- so this works
#' correctly regardless of that upstream toggle, with no spatial-join
#' reconstruction needed.
#'
#' @param data data.frame with the response, predictor, and fold columns.
#' @param gbmX Character vector of predictor column names.
#' @param gbmY Character. Response column name (must be 0/1).
#' @param foldID Integer vector, one fold id per row of `data`.
#' @param brtModel A fitted `gbm.step()` model (see `optimizeBRT()`) whose
#'   hyperparameters are reused for every fold's refit.
#' @return Numeric vector of out-of-fold predicted probabilities, same
#'   length and order as `data`/`foldID`.
blockCVPredictBRT <- function(data, gbmX, gbmY, foldID, brtModel) {
  cvPred <- rep(NA_real_, nrow(data))

  for (k in sort(unique(foldID))) {
    if (is.na(k)) next
    trainIdx <- which(foldID != k)
    testIdx <- which(foldID == k)
    if (length(trainIdx) == 0 || length(testIdx) == 0) next

    cvTrain <- data[trainIdx, ]
    cvTest <- data[testIdx, ]

    cvBrt <- gbm::gbm(formula = stats::as.formula(paste(gbmY, "~ .")),
                       distribution = "bernoulli",
                       data = cvTrain[, c(gbmY, gbmX)],
                       n.trees = brtModel$gbm.call$best.trees,
                       shrinkage = brtModel$gbm.call$learning.rate,
                       bag.fraction = brtModel$gbm.call$bag.fraction,
                       interaction.depth = brtModel$gbm.call$tree.complexity,
                       weights = rep(1, nrow(cvTrain)),
                       verbose = FALSE)

    cvPred[testIdx] <- gbm::predict.gbm(cvBrt, cvTest, type = "response",
                                         n.trees = brtModel$gbm.call$best.trees)
  }

  cvPred
}
