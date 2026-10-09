#' Band-wise prediction: the expensive part (habitat BRT at 200 m) and the ridge meta-model, per replicate
#'
#' The country window on the 200 m reference grid is cut into horizontal bands (uncBands()). For one band and
#' one year, `uncBandSuitability()` rebuilds exactly what the baseline's metaModel() does: predict each scale,
#' resample it BILINEARLY onto the reference grid, stack climate/landscape/habitat. `uncApplyRidge()` then
#' applies the replicate's own ridge coefficients.

#' Models of one species and the run's replicate ids, loaded once per task
uncContext <- function(cfg, sp) {
  m <- lapply(c("climate", "landscape", "habitat"), function(s)
    readRDS(file.path(uncSpDir(cfg, sp, "models"), paste0(s, "_", cfg$repLabel, ".rds"))))
  names(m) <- c("climate", "landscape", "habitat")
  list(models = m, ids = m$habitat$ids, coarseDir = uncSpDir(cfg, sp, "coarse", cfg$repLabel))
}

#' Suitability of one scale-triple for a set of replicates on a band template
#'
#' @param ctx From `uncContext()`.
#' @param year Integer.
#' @param template SpatRaster: the band (or a sub-window of a band) on the 200 m reference grid.
#' @param pos Integer positions (into `ctx$ids`) of the replicates to compute.
#' @param ptsXY Optional 2-column matrix of point coordinates: if given, the habitat BRT is only evaluated in the
#'   +-2-cell neighbourhood of these points (all that bilinear resampling needs at them) -- used to extract
#'   meta-model training values cheaply.
#' @return List of three numeric matrices (cells of `template` x replicates): `clim`, `land`, `hab`; NULL if a coarse
#'   scale has no prediction for this year.
uncBandSuitability <- function(cfg, sp, ctx, year, template, pos, ptsXY = NULL) {
  nCell <- terra::ncell(template); nRep <- length(pos)
  layerNames <- paste0("rep_", ctx$ids[pos])
  out <- list()

  # coarse scales: tiny rasters, resample the whole stack at once
  coarseScales <- c(clim = "climate", land = "landscape")
  for (key in names(coarseScales)) {
    nm <- coarseScales[[key]]
    f <- file.path(ctx$coarseDir, sprintf("%s_%d.tif", nm, year))
    # like the baseline's metaModel(): a year without one scale's prediction is skipped, not filled with NA
    if (!file.exists(f)) return(NULL)
    out[[key]] <- terra::values(terra::resample(terra::rast(f)[[layerNames]], template, method = "bilinear"), mat = TRUE)
  }

  # habitat: predict in a window around the template (padding for the bilinear neighbours), then resample
  covFile <- uncCovcacheFile(cfg, sp, year)
  if (!file.exists(covFile)) stop("Habitat covariate cache missing for ", year, ": ", covFile, " (run the covcache step)")
  covAll <- terra::rast(covFile)
  pad <- 4 * terra::res(covAll)[1]
  e <- terra::ext(template)
  padExt <- terra::intersect(terra::ext(terra::xmin(e) - pad, terra::xmax(e) + pad, terra::ymin(e) - pad, terra::ymax(e) + pad),
                             terra::ext(covAll))
  if (is.null(padExt)) { out$hab <- matrix(NA_real_, nCell, nRep); return(out) }
  covWin <- terra::crop(covAll, padExt, snap = "out")
  mod <- ctx$models$habitat
  missingPreds <- setdiff(mod$predSel, c(names(covWin), "x", "y"))
  if (length(missingPreds)) stop("Year ", year, ": habitat predictors missing in the cache: ", paste(missingPreds, collapse = ", "))
  keep <- if (is.null(ptsXY)) NULL else uncNeighbourCells(covWin, terra::cellFromXY(covWin, ptsXY), reach = 2L)
  cells <- uncCells(covWin, mod$predSel, keepCells = keep)
  vals <- matrix(NA_real_, terra::ncell(covWin), nRep)
  if (nrow(cells$df) > 0) vals[cells$idx, ] <- uncPredictScale(cfg, sp, "habitat", mod, pos, cells$df)
  hr <- terra::rast(covWin[[1]], nlyrs = nRep); terra::values(hr) <- vals
  out$hab <- terra::values(terra::resample(hr, template, method = "bilinear"), mat = TRUE)
  out
}

#' Apply replicate ridge coefficients to the suitability matrices
#'
#' @param S Output of `uncBandSuitability()`.
#' @param coefs Matrix (replicates x 4): intercept, climate, landscape, habitat (same order as the baseline's `suitCols`).
#' @return Matrix cells x replicates of probabilities (NA where any scale is NA).
uncApplyRidge <- function(S, coefs) {
  P <- matrix(NA_real_, nrow(S$hab), ncol(S$hab))
  for (j in seq_len(ncol(P)))
    P[, j] <- stats::plogis(coefs[j, 1] + coefs[j, 2] * S$clim[, j] + coefs[j, 3] * S$land[, j] + coefs[j, 4] * S$hab[, j])
  P
}

# ---- int16 storage ---------------------------------------------------------------------------------------

.uncScale <- 30000L

#' Save valid cells x replicates of probabilities as packed 16-bit integers (resolution 3.3e-5)
uncWriteInt16 <- function(path, idx, P, ids) {
  v <- as.integer(round(P * .uncScale)); v[is.na(v)] <- -1L
  uncSaveRDS(list(idx = as.integer(idx), dim = dim(P), ids = ids, raw = writeBin(v, raw(), size = 2)), path)
}

#' saveRDS via a temporary name + rename, so a job killed mid-write never leaves a half file that looks valid
uncSaveRDS <- function(obj, path) {
  tmp <- paste0(path, ".part")
  saveRDS(obj, tmp)
  file.rename(tmp, path)
  invisible(path)
}

#' Read it back: list(idx, P (double matrix valid cells x replicates), ids)
uncReadInt16 <- function(path) {
  x <- readRDS(path)
  v <- readBin(x$raw, "integer", n = prod(x$dim), size = 2, signed = TRUE)
  v[v < 0L] <- NA
  list(idx = x$idx, P = matrix(v / .uncScale, x$dim[1], x$dim[2]), ids = x$ids)
}

#' Predict one band and one year for all replicates of the run, store int16 + the area-mean sums
uncPredictBandYear <- function(cfg, sp, ctx, ridge, year, band) {
  predDir <- uncSpDir(cfg, sp, "pred", cfg$repLabel)
  outF <- file.path(predDir, sprintf("%d_band%02d.rds", year, band$k))
  statF <- file.path(predDir, sprintf("stats_%d_band%02d.rds", year, band$k))
  if (file.exists(outF) && file.exists(statF)) return(invisible("cached"))

  usable <- which(!is.na(ridge$coef[, 1]))
  batches <- split(usable, ceiling(seq_along(usable) / cfg$repBatch))
  Pall <- NULL; valid <- NULL; sums <- numeric(length(usable)); col <- 0L
  for (bt in batches) {
    S <- uncBandSuitability(cfg, sp, ctx, year, band$template, bt)
    if (is.null(S)) return(invisible("skipped (a scale has no prediction for this year)"))
    P <- uncApplyRidge(S, ridge$coef[bt, , drop = FALSE]); rm(S)
    if (is.null(valid)) {
      valid <- which(!is.na(P[, 1]))
      Pall <- matrix(NA_real_, length(valid), length(usable))
    }
    Pv <- P[valid, , drop = FALSE]; rm(P)
    if (anyNA(Pv)) warning(sp, " ", year, " band ", band$k, ": replicate-specific NA cells (", sum(is.na(Pv)), ") -- stored as NA")
    cols <- col + seq_along(bt); Pall[, cols] <- Pv; col <- col + length(bt)
  }
  if (is.null(valid)) valid <- integer(0)
  sums <- if (length(valid)) colSums(Pall, na.rm = TRUE) else rep(0, length(usable))
  uncWriteInt16(outF, valid, Pall, ctx$ids[usable])
  uncSaveRDS(list(ids = ctx$ids[usable], sum = sums, n = length(valid)), statF)
  invisible("done")
}

# ---- ridge meta-model per replicate -----------------------------------------------------------------------

#' Fit the ridge meta-model of every replicate of this run
#'
#' Training values: the replicate's own climate/landscape/habitat suitability at the habitat records of the
#' habitat years, OUT-OF-FOLD (each record predicted by a refit of the replicate's BRT that never saw its block,
#' uncOofSpecies()), so the weights are honest. Training rows: the
#' replicate's habitat bootstrap draw (the same blocks that refit the habitat BRT), so one resample of the
#' habitat data feeds both. Replicates fit with spatial-block folds for lambda (the baseline's plain 10-fold CV
#' would put copies of a resampled block in training and test folds); replicate 0 mimics the baseline exactly.
uncRidgeSpecies <- function(cfg, sp) {
  outF <- file.path(uncSpDir(cfg, sp, "ridge"), paste0("ridge_", cfg$repLabel, ".rds"))
  if (file.exists(outF)) { message("Ridge coefficients exist -- skipping: ", outF); return(invisible(readRDS(outF))) }
  ctx <- uncContext(cfg, sp)
  spPa <- uncTrainingTable(cfg, sp, "habitat")
  hy <- cfg$habitatYearsOf(sp)
  rowsAll <- which(spPa$year %in% hy)
  xy <- cbind(spPa$x, spPa$y)
  nRep <- length(ctx$ids)
  # HONEST inputs: out-of-fold predictions of each replicate's own scale models at the habitat records (uncOofSpecies()).
  oofF <- file.path(uncSpDir(cfg, sp, "oof"), paste0("oof_", cfg$repLabel, ".rds"))
  if (!file.exists(oofF)) stop("Out-of-fold inputs missing: ", oofF, " -- run the 'oof' step before 'ridge' (the replicate weights are trained on them).")
  oof <- readRDS(oofF)
  if (!identical(as.integer(oof$ids), as.integer(ctx$ids)) || !identical(as.integer(oof$rows), as.integer(rowsAll)))
    stop("Out-of-fold inputs do not match this run's replicates/records: ", oofF)
  Sx <- array(NA_real_, c(nrow(spPa), 3, nRep), dimnames = list(NULL, c("clim", "land", "hab"), NULL))
  Sx[rowsAll, , ] <- oof$S

  blockIds <- uncBlockIds(spPa$x, spPa$y, ctx$models$habitat$blockSizeM)
  coef <- matrix(NA_real_, nRep, 4, dimnames = list(paste0("rep_", ctx$ids), c("intercept", "climate", "landscape", "habitat")))
  info <- data.frame(replicate = ctx$ids, nTrain = NA_integer_, nPres = NA_integer_, nAbs = NA_integer_, lambda = NA_real_,
                     seed = NA_integer_, note = "", stringsAsFactors = FALSE)
  for (j in seq_len(nRep)) {
    b <- ctx$ids[j]
    idx <- if (b == 0) seq_len(nrow(spPa)) else ctx$models$habitat$draws[[j]]
    keep <- idx[spPa$year[idx] %in% hy]
    X <- cbind(Sx[keep, "clim", j], Sx[keep, "land", j], Sx[keep, "hab", j]); y <- spPa$occurrence[keep]
    ok <- stats::complete.cases(X); X <- X[ok, , drop = FALSE]; y <- y[ok]; blk <- blockIds[keep][ok]
    info$nTrain[j] <- length(y); info$nPres[j] <- sum(y == 1); info$nAbs[j] <- sum(y == 0)
    if (info$nPres[j] < 10 || info$nAbs[j] < 10) { info$note[j] <- "too few records -- replicate dropped"; next }
    colnames(X) <- c("climate_mean_prob", "landscape_mean_prob", "habitat_mean_prob")
    seed <- if (b == 0) 42L else stableSeed(c("ridgeBoot", sp, b))
    info$seed[j] <- seed
    # The penalty of a replicate is chosen by THE SAME RULE as the baseline's (ten random folds, lambda.1se), because a bootstrap must repeat the same procedure on the
    # resampled data. (Until 2026-10-10 the replicates used spatial BLOCK folds: 2-3 x larger penalties, weights shrunk towards zero -- Perdix perdix to exactly zero --
    # and trends much weaker than the baseline's, so the baseline fell outside its own interval; that variant is kept as a sensitivity analysis:
    # BIRDMONITOR_UNC_RIDGEFOLDS=block.)
    ridgeFolds <- Sys.getenv("BIRDMONITOR_UNC_RIDGEFOLDS", "random")
    fit <- tryCatch(uncWithSeed(seed, {
      if (b == 0 || ridgeFolds != "block") fitRidgeCv(X, y, nfolds = 10)
      else {
        ub <- unique(blk); nf <- min(10L, length(ub))
        foldOf <- stats::setNames(sample(rep_len(seq_len(nf), length(ub))), ub)
        fitRidgeCv(X, y, foldid = unname(foldOf[blk]))
      }
    }), error = function(e) { info$note[j] <<- conditionMessage(e); NULL })
    if (is.null(fit)) next
    info$lambda[j] <- fit$lambda.1se
    coef[j, ] <- as.vector(stats::coef(fit, s = fit$lambda.1se))
  }
  if (all(is.na(coef[, 1]))) stop("No replicate ridge model could be fitted for ", sp)
  res <- list(ids = ctx$ids, coef = coef, info = info, suitCols = c("climate_mean_prob", "landscape_mean_prob", "habitat_mean_prob"))
  uncSaveRDS(res, outF)
  message("Ridge: ", sum(!is.na(coef[, 1])), " of ", nRep, " replicate meta-models fitted -> ", outF)
  invisible(res)
}
