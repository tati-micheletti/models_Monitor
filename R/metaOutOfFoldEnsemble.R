#' Honest meta-model inputs for the ENSEMBLE: out-of-fold ensemble predictions at the habitat records
#'
#' The ensemble twin of `metaOutOfFoldCheck()`: each scale's out-of-fold prediction at the habitat records is the arithmetic mean of the
#' out-of-fold predictions of its members (BRT, GLM, GAM, random forest), prepared by `ensembleScaleRun()` as `<sp>_oofhab_ens_<scale>.rds`
#' (every member predicted each record with a model that never saw the record's block). The ridge is then fitted and evaluated exactly like
#' `metaOutOfFoldCheck()` does (same function `fitRidgeCv()`, same spatial folds).
#'
#' @param cfg Configuration from `uncCfgFromParams()`.
#' @param sp Character. Species.
#' @param tag Name of the ensemble (`ensTag(members)`: "ens", "ens_brt-gam", ...).
#' @return List with the same fields `metaModel()` uses from `metaOutOfFoldCheck()`: `perfOutOfFold`, `oof` (`X`, `y`, `foldID`),
#'   `coefficients`, `n`; `perfInSample` is NULL (the in-sample comparison is computed by `metaModel()` from the maps).
metaOutOfFoldEnsemble <- function(cfg, sp, tag = "ens") {
  hy <- cfg$habitatYearsOf(sp)
  hab <- uncTrainingTable(cfg, sp, "habitat")
  rows <- which(hab$year %in% hy)
  y <- hab$occurrence[rows]; foldHab <- hab$foldID[rows]
  suitCols <- c("climate_mean_prob", "landscape_mean_prob", "habitat_mean_prob")
  X <- do.call(cbind, lapply(c("climate", "landscape", "habitat"), function(s) {
    f <- ensFile(cfg, sp, s, tag, "oofhab")
    if (!file.exists(f)) stop("Ensemble out-of-fold inputs missing: ", f, " (run ensembleScaleRun() first)")
    v <- readRDS(f)
    if (length(v) != length(rows)) stop(basename(f), " has ", length(v), " values for ", length(rows), " habitat records")
    v }))
  colnames(X) <- suitCols
  ok <- stats::complete.cases(X); Xk <- X[ok, , drop = FALSE]; yk <- y[ok]
  set.seed(42)
  cvFit <- fitRidgeCv(Xk, yk, nfolds = 10)
  cvPred <- blockCVPredictRidge(Xk, yk, foldHab[ok])
  v <- !is.na(cvPred)
  coefs <- matrix(as.vector(stats::coef(cvFit, s = cvFit$lambda.1se)), 1, dimnames = list("outOfFold", c("intercept", "climate", "landscape", "habitat")))
  list(species = sp, n = sum(ok), coefficients = coefs, perfOutOfFold = evalSDM(yk[v], cvPred[v]), perfInSample = NULL,
       oof = list(X = Xk, y = yk, foldID = foldHab[ok]))
}
