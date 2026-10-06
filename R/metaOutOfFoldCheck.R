#' Honest accuracy of the combined (meta-model) map: refit the meta-model on OUT-OF-FOLD inputs and compare
#'
#' Why: `metaModel()` fits its ridge regression on the three scale predictions at the habitat records. Those predictions come
#' from BRTs fitted on ALL records, so the combiner's own block cross-validation (AUC ~0.99) is flattered: by the time a block is
#' held out for the combiner, its inputs already contain what the BRTs learned from that very block (improvements.md item 16).
#'
#' What this does, per species, on the habitat records of the habitat years:
#'   in-sample inputs  = what the pipeline does: the MAIN model of each scale predicts the record.
#'   out-of-fold inputs = each record is predicted by a FOLD model that never saw it:
#'       habitat   : the fold model of the record's own fold (`foldID` of the habitat table);
#'       landscape, climate : the fold model of the fold of the NEAREST record of that scale's table (a record that sits in
#'                   a block is held out together with its block; the blocks themselves were not saved, so this is the
#'                   nearest-neighbour approximation).
#' Then the ridge is fitted on both versions exactly like `fitRidgeOneSpecies()` and evaluated with the habitat table's own
#' block folds (`blockCVPredictRidge()` + `evalSDM()`). Predictions of the coarse scales are made in the neighbourhood of the
#' records only and read with bilinear interpolation, like the pipeline's resampling to the 200 m grid.
#'
#' @param cfg Uncertainty/baseline configuration from `uncCfgFromParams()`.
#' @param sp Character. Species.
#' @return List: `performance` (data.frame, one row per variant), `coefficients` (matrix, incl. the pipeline's own as a
#'   sanity check against the in-sample refit), `n` (records), `perfOutOfFold` / `perfInSample` (the two meta-model
#'   `evalSDM()` rows; `perfOutOfFold` is what `metaModel()` reports).
metaOutOfFoldCheck <- function(cfg, sp) {
  spClean <- gsub(" ", "_", sp)
  hy <- cfg$habitatYearsOf(sp)
  tabs <- list(habitat = uncTrainingTable(cfg, sp, "habitat"), landscape = uncTrainingTable(cfg, sp, "landscape"),
               climate = uncTrainingTable(cfg, sp, "climate"))
  hab <- tabs$habitat
  rows <- which(hab$year %in% hy)
  y <- hab$occurrence[rows]; yrs <- hab$year[rows]; foldHab <- hab$foldID[rows]
  xy <- cbind(hab$x[rows], hab$y[rows])
  mains <- lapply(c(climate = "climate", landscape = "landscape", habitat = "habitat"), function(s) uncMainModel(cfg, sp, s))
  folds <- lapply(c(climate = "climate", landscape = "landscape", habitat = "habitat"), function(s) {
    f <- file.path(uncMainDir(cfg, sp, s), paste0(spClean, "_foldModels_", s, ".rds"))
    if (!file.exists(f)) stop("Fold models missing: ", f, " (they are saved by the model arrays)")
    readRDS(f)
  })
  nTrees <- function(m) if (!is.null(m$gbm.call$best.trees)) m$gbm.call$best.trees else m$n.trees

  # fold of the nearest record of a scale's own table (stands in for "the block this location belongs to")
  nearestFold <- function(tbl) {
    pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = 3035)
    ref <- sf::st_as_sf(data.frame(x = tbl$x, y = tbl$y), coords = c("x", "y"), crs = 3035)
    tbl$foldID[as.integer(sf::st_nearest_feature(pts, ref))]
  }
  foldOf <- list(climate = nearestFold(tabs$climate), landscape = nearestFold(tabs$landscape))

  # --- coarse scales: predict around the records, read bilinearly --------------------------------------------------------
  inSample <- matrix(NA_real_, length(rows), 3, dimnames = list(NULL, c("climate", "landscape", "habitat")))
  oof <- inSample
  for (scale in c("climate", "landscape")) {
    lab <- scaleLabel(cfg$resOf(sp)[[scale]])
    procDir <- file.path(cfg$inputRoot, "predictors", "processed", lab)
    predSel <- mains[[scale]]$var.names
    for (yr in unique(yrs)) {
      ry <- which(yrs == yr)
      cov <- if (scale == "climate") {
        terra::rast(file.path(procDir, paste0("bioclim_", yr - (cfg$climateWindowLength - 1), "-", yr, "_", lab, ".tif")))
      } else {
        cv <- loadCovariates(yr, procDir, NULL); names(cv) <- gsub("_[0-9]{4}$", "", names(cv)); cv
      }
      keep <- uncNeighbourCells(cov, terra::cellFromXY(cov, xy[ry, , drop = FALSE]), reach = 2L)
      cells <- uncCells(cov, predSel, keepCells = keep)
      readAt <- function(model, rr) {
        r <- terra::rast(cov[[1]]); v <- rep(NA_real_, terra::ncell(r))
        v[cells$idx] <- gbm::predict.gbm(model, cells$df, n.trees = nTrees(model), type = "response")
        terra::values(r) <- v
        terra::extract(r, xy[rr, , drop = FALSE], method = "bilinear")[, 1]
      }
      inSample[ry, scale] <- readAt(mains[[scale]], ry)
      for (k in unique(foldOf[[scale]][ry])) {
        rr <- ry[foldOf[[scale]][ry] == k]
        oof[rr, scale] <- readAt(folds[[scale]][[as.character(k)]], rr)
      }
    }
  }

  # --- habitat: predict at the record's own predictors --------------------------------------------------------------------
  predHab <- mains$habitat$var.names
  inSample[, "habitat"] <- gbm::predict.gbm(mains$habitat, hab[rows, predHab, drop = FALSE], n.trees = nTrees(mains$habitat), type = "response")
  for (k in unique(foldHab)) {
    rr <- which(foldHab == k)
    oof[rr, "habitat"] <- gbm::predict.gbm(folds$habitat[[as.character(k)]], hab[rows[rr], predHab, drop = FALSE],
                                           n.trees = nTrees(folds$habitat[[as.character(k)]]), type = "response")
  }

  # --- ridge exactly like fitRidgeOneSpecies(), evaluated like evalRidgeOneSpecies() ---------------------------------------
  suitCols <- c("climate_mean_prob", "landscape_mean_prob", "habitat_mean_prob")
  fitAndEval <- function(X) {
    ok <- stats::complete.cases(X); Xk <- X[ok, , drop = FALSE]; yk <- y[ok]
    colnames(Xk) <- suitCols
    set.seed(42)
    cvFit <- fitRidgeCv(Xk, yk, nfolds = 10)
    cvPred <- blockCVPredictRidge(Xk, yk, foldHab[ok])
    v <- !is.na(cvPred)
    list(perf = evalSDM(yk[v], cvPred[v]), coef = as.vector(stats::coef(cvFit, s = cvFit$lambda.1se)), n = sum(ok))
  }
  fIn <- fitAndEval(inSample); fOof <- fitAndEval(oof)
  single <- function(v, label) { ok <- !is.na(v); cbind(data.frame(variant = label), evalSDM(y[ok], v[ok])) }
  perf <- rbind(
    single(oof[, "habitat"], "habitat BRT alone, out-of-fold (honest)"),
    single(oof[, "landscape"], "landscape BRT alone, out-of-fold"),
    single(oof[, "climate"], "climate BRT alone, out-of-fold"),
    single(inSample[, "habitat"], "habitat BRT alone, in-sample (optimistic)"),
    cbind(data.frame(variant = "meta-model, in-sample inputs (what the pipeline does)"), fIn$perf),
    cbind(data.frame(variant = "meta-model, OUT-OF-FOLD inputs (honest)"), fOof$perf))
  coefs <- rbind(inSample = fIn$coef, outOfFold = fOof$coef)
  colnames(coefs) <- c("intercept", "climate", "landscape", "habitat")
  pipeFile <- list.files(cfg$outputRoot, pattern = paste0("^", spClean, "_ridge_meta[.]rds$"), recursive = TRUE, full.names = TRUE)
  pipe <- if (length(pipeFile)) { m <- readRDS(pipeFile[1]); as.vector(stats::coef(m$model, s = m$lambda)) } else rep(NA_real_, 4)
  coefs <- rbind(coefs, pipeline = pipe)
  list(species = sp, n = fIn$n, performance = cbind(species = sp, perf), coefficients = coefs,
       perfOutOfFold = fOof$perf, perfInSample = fIn$perf)
}
