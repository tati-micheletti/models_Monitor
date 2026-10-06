#' The ONE place the ridge meta-model is fitted (cross-validated glmnet, alpha = 0)
#'
#' Every ridge fit of the pipeline goes through here -- the main meta-model (`fitRidgeOneSpecies()`), its block-CV
#' evaluation (`blockCVPredictRidge()`), the trend bootstrap, the uncertainty replicates (`uncRidgeSpecies()`) and the
#' out-of-fold check -- so they cannot drift apart.
#'
#' `lower.limits = 0`: a scale can never get a NEGATIVE weight. A combiner of three predictions of the same species'
#' occurrence has no biological reason to subtract one of them; a negative weight is the ridge compensating for
#' correlation between the inputs, and it makes the maps hard to defend. Wiedenroth et al. do the same ("we did not allow
#' negative coefficients in the model"). A scale whose weight would be negative is simply switched off (weight 0).
#'
#' @param x Numeric matrix, the three scale predictions.
#' @param y Numeric vector, 0/1 occurrence.
#' @param ... Passed to `glmnet::cv.glmnet()` (`nfolds` or `foldid`).
#' @return A `cv.glmnet` fit.
fitRidgeCv <- function(x, y, ...) {
  glmnet::cv.glmnet(x = x, y = y, family = "binomial", alpha = 0, standardize = TRUE, lower.limits = 0, ...)
}
