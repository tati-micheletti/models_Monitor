#' Alternative model families for one scale: GLM, GAM and random forest (the BRT stays in optimizeBRT())
#'
#' Wiedenroth et al. fit four algorithms at every scale (04c_model-training_200m.R) and average them (arithmetic mean) into the
#' scale's ensemble prediction. These functions follow their specification:
#'   glm : binomial GLM with a linear and a quadratic term for every predictor, AIC stepwise selection (`step()`); if the
#'         selected model does not converge, terms are dropped from the end until it does.
#'   gam : `mgcv::gam`, binomial, one smooth `s(x, k = 4)` per predictor; if the fit errors (too few distinct values in a
#'         predictor), the predictor with the fewest non-zero values is dropped and the fit retried.
#'   rf  : `randomForest` on the 0/1 response (regression forest, i.e. the predicted probability is the forest mean),
#'         1000 trees.
#' Every fitted model is returned as a small list (`algo`, `model`, `formula`, `predSel` = the predictors actually used) so the
#' rest of the pipeline treats the algorithms alike: `algoPredict()`, `algoRefit()` (same specification, new data: used for
#' block cross-validation and bootstrap replicates) and `blockCVPredictAlgo()`.

#' Fit one algorithm on a training table
#'
#' @param algo "glm", "gam" or "rf".
#' @param data data.frame with `response` (0/1) and the predictors.
#' @param predSel Character. Candidate predictors.
#' @param response Character. Response column (default "occurrence").
#' @param ntree Integer. Trees of the random forest.
#' @param seed Integer or NULL. Set before fitting (random forest, bagging).
#' @return List: `algo`, `model`, `formula`, `predSel` (predictors the final model uses), `response`, `ntree`, `note`.
algoFit <- function(algo, data, predSel, response = "occurrence", ntree = 1000L, seed = NULL) {
  algo <- match.arg(algo, c("glm", "gam", "rf"))
  if (!is.null(seed)) set.seed(seed)
  switch(algo,
    glm = algoFitGLM(data, predSel, response),
    gam = algoFitGAM(data, predSel, response),
    rf  = algoFitRF(data, predSel, response, ntree))
}

algoFitGLM <- function(data, predSel, response) {
  terms <- unlist(lapply(predSel, function(k) c(k, paste0("I(", k, "^2)"))))
  f <- stats::as.formula(paste(response, "~", paste(terms, collapse = " + ")))
  m <- stats::step(stats::glm(f, family = "binomial", data = data), trace = 0)
  note <- ""
  # Did not converge: drop the last coefficient until it does (as in the reference code).
  n <- 0L
  while (!isTRUE(m$converged)) {
    keep <- names(m$coefficients)[-1]
    if (length(keep) <= 1L) stop("GLM did not converge even after removing terms.")
    keep <- keep[-length(keep)]
    note <- paste0(note, "dropped ", names(m$coefficients)[length(m$coefficients)], "; ")
    m <- stats::glm(stats::as.formula(paste(response, "~", paste(keep, collapse = " + "))), family = "binomial", data = data)
    n <- n + 1L
  }
  used <- predSel[vapply(predSel, function(k) any(grepl(paste0("(^|\\()", k, "(\\)|\\^|$)"), names(m$coefficients))), logical(1))]
  list(algo = "glm", model = m, formula = stats::formula(m), predSel = used, response = response, note = note)
}

algoFitGAM <- function(data, predSel, response) {
  used <- predSel; note <- ""
  repeat {
    f <- stats::as.formula(paste(response, "~", paste0("s(", used, ", k = 4)", collapse = " + ")))
    m <- tryCatch(mgcv::gam(f, family = "binomial", data = data), error = function(e) NULL)
    if (!is.null(m)) break
    if (length(used) <= 1L) stop("GAM could not be fitted even with a single predictor.")
    nz <- vapply(used, function(k) sum(data[[k]] != 0, na.rm = TRUE), numeric(1))
    drop <- used[which.min(nz)]
    note <- paste0(note, "dropped ", drop, "; ")
    used <- setdiff(used, drop)
  }
  if (!isTRUE(m$converged)) stop("GAM did not converge.")
  list(algo = "gam", model = m, formula = f, predSel = used, response = response, note = note)
}

algoFitRF <- function(data, predSel, response, ntree) {
  f <- stats::as.formula(paste(response, "~", paste(predSel, collapse = " + ")))
  m <- suppressWarnings(randomForest::randomForest(formula = f, data = data[, c(response, predSel), drop = FALSE], ntree = ntree))   # regression forest on 0/1, as in the reference code (randomForest warns about it)
  list(algo = "rf", model = m, formula = f, predSel = predSel, response = response, ntree = ntree, note = "")
}

#' Predicted probability of presence
#'
#' @param fit A list returned by `algoFit()`.
#' @param newdata data.frame with at least `fit$predSel`.
#' @return Numeric vector, one value per row of `newdata`.
algoPredict <- function(fit, newdata) {
  nd <- newdata[, fit$predSel, drop = FALSE]
  p <- switch(fit$algo,
    glm = stats::predict(fit$model, nd, type = "response"),
    gam = stats::predict(fit$model, nd, type = "response"),
    rf  = stats::predict(fit$model, nd, type = "response"))
  pmin(pmax(as.numeric(p), 0), 1)
}

#' Refit the SAME specification on other data (no new model selection)
#'
#' Used for the block cross-validation and for bootstrap replicates, exactly like the reference code's `update(model, data = ...)`:
#' the selected terms/predictors and the settings stay fixed.
algoRefit <- function(fit, data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- switch(fit$algo,
    glm = stats::glm(fit$formula, family = "binomial", data = data),
    gam = mgcv::gam(fit$formula, family = "binomial", data = data),
    rf  = suppressWarnings(randomForest::randomForest(formula = fit$formula, data = data[, c(fit$response, fit$predSel), drop = FALSE], ntree = fit$ntree)))
  fit$model <- m
  fit
}

#' Block cross-validated predictions of one algorithm
#'
#' Same scheme as `blockCVPredictBRT()`: per fold, refit on the other folds, predict the held-out fold.
#'
#' @param fit A list returned by `algoFit()`.
#' @param data The training table (response + predictors).
#' @param foldID Integer vector, one fold per row of `data`.
#' @return Numeric vector of out-of-fold probabilities (NA where a fold could not be fitted).
blockCVPredictAlgo <- function(fit, data, foldID) {
  cvPred <- rep(NA_real_, nrow(data))
  for (k in sort(unique(foldID))) {
    if (is.na(k)) next
    trainIdx <- which(foldID != k); testIdx <- which(foldID == k)
    if (!length(trainIdx) || !length(testIdx)) next
    fk <- tryCatch(algoRefit(fit, data[trainIdx, , drop = FALSE], seed = stableSeed(c("cvAlgo", fit$algo, k, nrow(data), sum(data[[fit$response]])))),
                   error = function(e) NULL)
    if (is.null(fk)) next
    cvPred[testIdx] <- algoPredict(fk, data[testIdx, , drop = FALSE])
  }
  cvPred
}

#' Arithmetic mean of several algorithms' predictions (the ensemble of one scale)
#'
#' @param preds List (or data.frame) of equal-length numeric vectors, one per algorithm.
#' @return Numeric vector. NA only where every algorithm is NA; otherwise the mean of the available ones is NOT used -- a cell
#'   needs all algorithms, so the ensemble never silently becomes a different mixture from one cell to the next.
algoEnsembleMean <- function(preds) {
  m <- do.call(cbind, preds)
  out <- rowMeans(m)
  out[!stats::complete.cases(m)] <- NA_real_
  out
}
