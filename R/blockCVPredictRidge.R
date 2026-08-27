#' Block cross-validated predictions from a ridge regression
#'
#' Nested CV: for each outer fold (from `foldID` -- see `blockCVPredictBRT()`
#' for why this is sourced directly from inputs_Monitor's table rather than
#' reconstructed via a spatial join), fits a fresh `cv.glmnet()` (its own
#' inner 5-fold CV to pick lambda) on the out-of-fold rows and predicts on
#' the held-out fold.
#'
#' @param X Numeric matrix of predictors.
#' @param y Numeric vector of 0/1 responses.
#' @param foldID Integer vector, one fold id per row of `X`/`y`.
#' @param innerNfolds Integer. `cv.glmnet()`'s own inner CV folds, default 5.
#' @return Numeric vector of out-of-fold predicted probabilities.
blockCVPredictRidge <- function(X, y, foldID, innerNfolds = 5) {
  cvPred <- rep(NA_real_, nrow(X))

  for (k in sort(unique(foldID))) {
    if (is.na(k)) next
    trainIdx <- which(foldID != k)
    testIdx <- which(foldID == k)
    if (length(trainIdx) == 0 || length(testIdx) == 0) next

    cvFit <- glmnet::cv.glmnet(x = X[trainIdx, , drop = FALSE], y = y[trainIdx],
                                family = "binomial", alpha = 0, nfolds = innerNfolds,
                                standardize = TRUE)

    cvPred[testIdx] <- as.vector(stats::predict(cvFit, newx = X[testIdx, , drop = FALSE],
                                                 s = cvFit$lambda.1se, type = "response"))
  }

  cvPred
}
