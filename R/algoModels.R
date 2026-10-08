#' Model families of the multi-algorithm ensemble: a REGISTRY, so a model is switched on or off by name and a new one is a single entry
#'
#' Wiedenroth et al. fit four algorithms at every scale (04c_model-training_200m.R) and average them (arithmetic mean) into the scale's
#' ensemble prediction. Here every family is a registry entry (`algoRegistry()`) with three functions -- `fit`, `refit` (same specification,
#' new data: block cross-validation, bootstrap) and `predict` -- and the whole workflow (`algoScaleRun()`, `ensembleScaleRun()`, the SLURM
#' scripts) only ever uses the names. The BRT is a member like the others (its outputs come from the existing BRT workflow; see
#' `brtMemberRun()`); the names, in the canonical order, are `ALGO_MEMBERS = c("brt", "glm", "gam", "rf", "nn")`.
#'
#' TO ADD A MODEL: write `fitXxx(data, predSel, response, ..., seed, threads)` returning `list(algo = "xxx", model, predSel, response, note, ...)`,
#' add a `refit` and a `predict` function, put the three in `algoRegistry()` and the name in `ALGO_MEMBERS`. Nothing else changes.
#'
#'   glm : binomial GLM with a linear and a quadratic term for every predictor, AIC stepwise selection (`step()`); if the selected model does
#'         not converge, terms are dropped from the end until it does.
#'   gam : `mgcv::gam`, binomial, one smooth `s(x, k = 4)` per predictor; if the fit errors (too few distinct values in a predictor), the
#'         predictor with the fewest non-zero values is dropped and the fit retried.
#'   rf  : random forest on the 0/1 response, grown as a REGRESSION forest (the predicted probability is the forest mean), 1000 trees.
#'         Wiedenroth et al. use `randomForest::randomForest()`; here `ranger::ranger()` does the same job several times faster and can use
#'         several cores (a documented deviation, decided 2026-10-07), with `randomForest`'s regression defaults (`mtry = max(floor(p/3), 1)`,
#'         `min.node.size = 5`, bootstrap sampling).
#'   nn  : feed-forward neural network (`nnet`, one hidden layer, logistic output, weight decay), predictors standardised, the average of
#'         `nRep` fits from different random starts (a single start is unstable). FIRST VERSION with fixed size/decay (not tuned); an addition
#'         beyond Wiedenroth et al.

ALGO_MEMBERS <- c("brt", "glm", "gam", "rf", "nn")

#' Fit one family on a training table (see the registry)
#'
#' @param algo One of `ALGO_MEMBERS` except "brt" (the BRT is fitted by the existing BRT workflow).
#' @param data data.frame with `response` (0/1) and the predictors.
#' @param predSel Character. Candidate predictors.
#' @param response Character. Response column (default "occurrence").
#' @param ntree Integer. Trees of the random forest.
#' @param threads Integer. Threads of the random forest (fit and prediction).
#' @param seed Integer or NULL. Set before fitting.
#' @return List: `algo`, `model`, `predSel` (predictors the final model uses), `response`, `note` (+ family-specific fields).
algoFit <- function(algo, data, predSel, response = "occurrence", ntree = 1000L, seed = NULL, threads = 1L) {
  reg <- algoRegistry()
  if (!algo %in% names(reg)) stop("algoFit(): unknown or non-fittable algorithm '", algo, "' (available: ", paste(names(reg), collapse = ", "), ")")
  if (!is.null(seed)) set.seed(seed)
  reg[[algo]]$fit(data = data, predSel = predSel, response = response, ntree = ntree, seed = seed, threads = threads)
}

#' Predicted probability of presence
#'
#' @param fit A list returned by `algoFit()`.
#' @param newdata data.frame with at least `fit$predSel`.
#' @return Numeric vector, one value per row of `newdata`.
algoPredict <- function(fit, newdata) {
  p <- algoRegistry()[[fit$algo]]$predict(fit, newdata[, fit$predSel, drop = FALSE])
  pmin(pmax(as.numeric(p), 0), 1)
}

#' Refit the SAME specification on other data (no new model selection)
#'
#' Used for the block cross-validation and for bootstrap replicates, like the reference code's `update(model, data = ...)`: the selected terms,
#' predictors and settings stay fixed.
algoRefit <- function(fit, data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  algoRegistry()[[fit$algo]]$refit(fit, data, seed)
}

#' The registry: name -> list(fit, refit, predict)
algoRegistry <- function() list(
  glm = list(fit = function(data, predSel, response, ...) algoFitGLM(data, predSel, response),
             refit = function(fit, data, seed) { fit$model <- stats::glm(fit$formula, family = "binomial", data = data); fit },
             predict = function(fit, nd) stats::predict(fit$model, nd, type = "response")),
  gam = list(fit = function(data, predSel, response, ...) algoFitGAM(data, predSel, response),
             refit = function(fit, data, seed) { fit$model <- mgcv::gam(fit$formula, family = "binomial", data = data); fit },
             predict = function(fit, nd) { loadNamespace("mgcv"); stats::predict(fit$model, nd, type = "response") }),   # S3 method needs the package loaded
  rf  = list(fit = function(data, predSel, response, ntree, seed, threads, ...) algoFitRF(data, predSel, response, ntree, seed, threads),
             refit = function(fit, data, seed) { fit$model <- algoFitRF(data, fit$predSel, fit$response, fit$ntree, seed, if (is.null(fit$threads)) 1L else fit$threads)$model; fit },
             predict = function(fit, nd) ranger::predictions(stats::predict(fit$model, data = nd, num.threads = if (is.null(fit$threads)) 1L else fit$threads))),
  nn  = list(fit = function(data, predSel, response, seed, ...) algoFitNN(data, predSel, response, seed),
             refit = function(fit, data, seed) { n <- algoFitNN(data, fit$predSel, fit$response, seed, size = fit$size, decay = fit$decay, nRep = fit$nRep); fit$model <- n$model; fit },
             predict = function(fit, nd) algoPredictNN(fit, nd))
)

algoFitGLM <- function(data, predSel, response) {
  terms <- unlist(lapply(predSel, function(k) c(k, paste0("I(", k, "^2)"))))
  f <- stats::as.formula(paste(response, "~", paste(terms, collapse = " + ")))
  m <- stats::step(stats::glm(f, family = "binomial", data = data), trace = 0)
  note <- ""
  # Did not converge: drop the last coefficient until it does (as in the reference code).
  while (!isTRUE(m$converged)) {
    keep <- names(m$coefficients)[-1]
    if (length(keep) <= 1L) stop("GLM did not converge even after removing terms.")
    keep <- keep[-length(keep)]
    note <- paste0(note, "dropped ", names(m$coefficients)[length(m$coefficients)], "; ")
    m <- stats::glm(stats::as.formula(paste(response, "~", paste(keep, collapse = " + "))), family = "binomial", data = data)
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

algoFitRF <- function(data, predSel, response, ntree, seed = NULL, threads = 1L) {
  f <- stats::as.formula(paste(response, "~", paste(predSel, collapse = " + ")))
  m <- ranger::ranger(formula = f, data = data[, c(response, predSel), drop = FALSE], num.trees = ntree,
                      mtry = max(floor(length(predSel) / 3), 1L), min.node.size = 5L, num.threads = threads, seed = seed)
  list(algo = "rf", model = m, formula = f, predSel = predSel, response = response, ntree = ntree, threads = threads, note = "")
}

# neural network: standardised predictors, the mean of `nRep` single-hidden-layer fits (logistic output, cross-entropy, weight decay)
# Defaults (2026-10-08, tools/diagnoseNN.R pilot, block-CV, 3 scales): 3 hidden units, decay 0.1, up to 1000 iterations. The earlier 10 units / 300 iterations
# did not converge in 15 of 15 landscape fits; 3 units / decay 0.1 matched or beat the best setting at every scale and converged.
algoFitNN <- function(data, predSel, response, seed = NULL, size = 3L, decay = 0.1, nRep = 5L, maxit = 1000L) {
  X <- as.matrix(data[, predSel, drop = FALSE])
  mu <- colMeans(X); sdv <- apply(X, 2, stats::sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
  Xs <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
  base <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed
  fits <- lapply(seq_len(nRep), function(i) {
    set.seed(base + i)
    nnet::nnet(x = Xs, y = data[[response]], size = size, decay = decay, entropy = TRUE, maxit = maxit, MaxNWts = 10000L, trace = FALSE)
  })
  list(algo = "nn", model = list(fits = fits, mu = mu, sdv = sdv), predSel = predSel, response = response, size = size, decay = decay,
       nRep = nRep, note = "")
}
algoPredictNN <- function(fit, nd) {
  loadNamespace("nnet")      # predict() on a saved fit needs the package's S3 method, also in a fresh R process
  Xs <- sweep(sweep(as.matrix(nd[, fit$predSel, drop = FALSE]), 2, fit$model$mu, "-"), 2, fit$model$sdv, "/")
  rowMeans(vapply(fit$model$fits, function(m) as.numeric(stats::predict(m, Xs, type = "raw")), numeric(nrow(Xs))))
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
#' @param preds List (or data.frame) of equal-length numeric vectors, one per member.
#' @return Numeric vector, NA where ANY member is NA: a cell needs all members, so the ensemble never silently becomes a different mixture
#'   from one cell to the next.
algoEnsembleMean <- function(preds) {
  m <- do.call(cbind, preds)
  out <- rowMeans(m)
  out[!stats::complete.cases(m)] <- NA_real_
  out
}

#' Name of an ensemble made of `members`
#'
#' The full default set (brt, glm, gam, rf) is called "ens" (as before); any other selection is "ens_<members in canonical order joined by
#' '-'>", e.g. "ens_brt-gam" or "ens_brt-glm-gam-rf-nn". The name is used in file names, in the meta-model folder (`metamodel_<label>_<name>`)
#' and as the index output tag, so different ensembles never overwrite each other.
ensTag <- function(members) {
  members <- ALGO_MEMBERS[ALGO_MEMBERS %in% members]
  if (!length(members)) stop("ensTag(): no valid members")
  if (identical(members, c("brt", "glm", "gam", "rf"))) "ens" else paste0("ens_", paste(members, collapse = "-"))
}
