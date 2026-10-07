#' Out-of-fold scale predictions at the habitat records, per replicate: the HONEST inputs of the replicate ridge models
#'
#' The meta-model weights must be trained on predictions the scale models made for records they never saw, otherwise the
#' scale that overfits most (habitat) gets too much weight (stacked generalization: Wolpert 1992, Breiman 1996). The baseline
#' does this with the saved fold models of the main BRTs (`metaOutOfFoldCheck()`). A replicate has no saved fold models, so
#' here they are refitted, with the replicate's own bootstrap draw:
#'   fold k model of replicate b, scale s = the replicate's BRT, refitted with the main model's fixed hyperparameters on the
#'   rows of the draw whose block-CV fold (`foldID` of the scale's own table) is not k;
#'   it predicts the habitat records that belong to fold k (habitat: the record's own fold; climate/landscape: the fold of
#'   the nearest record of that scale's table, the same approximation as the baseline).
#' Replicate 0 uses the saved fold models of the main BRTs, so it equals the baseline.

#' Everything that does not depend on the replicate: the habitat records, their folds, cached covariates/cells
uncOofContext <- function(cfg, sp) {
  hy <- cfg$habitatYearsOf(sp)
  tabs <- list(habitat = uncTrainingTable(cfg, sp, "habitat"), landscape = uncTrainingTable(cfg, sp, "landscape"),
               climate = uncTrainingTable(cfg, sp, "climate"))
  hab <- tabs$habitat
  rows <- which(hab$year %in% hy)
  xy <- cbind(hab$x[rows], hab$y[rows])
  nearestFold <- function(tbl) {
    tbl$foldID[uncNearestIdx(xy, cbind(tbl$x, tbl$y))]
  }
  list(cfg = cfg, sp = sp, tabs = tabs, hab = hab, rows = rows, xy = xy, y = hab$occurrence[rows], yrs = hab$year[rows],
       foldOf = list(climate = nearestFold(tabs$climate), landscape = nearestFold(tabs$landscape), habitat = hab$foldID[rows]),
       cache = new.env(parent = emptyenv()))
}

#' Number of trees of a main-model fold fit or of a replicate fit
uncOofNTrees <- function(m) if (!is.null(m$gbm.call$best.trees)) m$gbm.call$best.trees else m$n.trees

#' Covariates of a coarse scale and year (cached in the context), plus the cells near the records of fold `k`
uncOofCells <- function(ctx, scale, yr, k, predSel, pos) {
  key <- paste(scale, yr, k, sep = "|")
  if (!is.null(ctx$cache[[key]])) return(ctx$cache[[key]])
  cfg <- ctx$cfg
  covKey <- paste("cov", scale, yr, sep = "|")
  if (is.null(ctx$cache[[covKey]])) {
    lab <- scaleLabel(cfg$resOf(ctx$sp)[[scale]])
    procDir <- file.path(cfg$inputRoot, "predictors", "processed", lab)
    ctx$cache[[covKey]] <- if (scale == "climate") {
      terra::rast(file.path(procDir, paste0("bioclim_", yr - (cfg$climateWindowLength - 1), "-", yr, "_", lab, ".tif")))
    } else {
      cv <- loadCovariates(yr, procDir, NULL); names(cv) <- gsub("_[0-9]{4}$", "", names(cv)); cv
    }
  }
  cov <- ctx$cache[[covKey]]
  keep <- uncNeighbourCells(cov, terra::cellFromXY(cov, ctx$xy[pos, , drop = FALSE]), reach = 2L)
  out <- list(cov = cov, cells = uncCells(cov, predSel, keepCells = keep), pos = pos)
  ctx$cache[[key]] <- out
  out
}

#' Predict the habitat records of fold `k` (positions in `ctx$rows`) with one fold model of one scale
#'
#' @param predictFun NULL (the model is a `gbm`) or a function(data.frame) -> probabilities, which is used INSTEAD of the gbm
#'   prediction (GLM/GAM/random-forest fold models of the ensemble: `function(df) algoPredict(foldFit, df)`).
#' @return List: `pos` (positions in `ctx$rows`) and `value` (predicted probability).
uncOofPredictFold <- function(ctx, scale, model, predSel, k, predictFun = NULL) {
  pos <- which(ctx$foldOf[[scale]] == k)
  if (!length(pos)) return(list(pos = integer(0), value = numeric(0)))
  nt <- if (is.null(predictFun)) uncOofNTrees(model) else NA
  predictRows <- function(df) if (is.null(predictFun)) gbm::predict.gbm(model, df, n.trees = nt, type = "response") else predictFun(df)
  if (scale == "habitat") {
    return(list(pos = pos, value = predictRows(ctx$hab[ctx$rows[pos], predSel, drop = FALSE])))
  }
  value <- rep(NA_real_, length(pos))
  for (yr in unique(ctx$yrs[pos])) {
    pp <- pos[ctx$yrs[pos] == yr]
    cc <- uncOofCells(ctx, scale, yr, k, predSel, pp)
    r <- terra::rast(cc$cov[[1]]); v <- rep(NA_real_, terra::ncell(r))
    v[cc$cells$idx] <- predictRows(cc$cells$df)
    terra::values(r) <- v
    value[match(pp, pos)] <- terra::extract(r, ctx$xy[pp, , drop = FALSE], method = "bilinear")[, 1]
  }
  list(pos = pos, value = value)
}

#' Out-of-fold inputs of every replicate of this run, for one species
#'
#' Writes `oof/oof_<run>.rds`: `S` (habitat records of the habitat years x 3 scales x replicates), `ids`, `rows`.
uncOofSpecies <- function(cfg, sp, minPerClass = 10) {
  spClean <- gsub(" ", "_", sp)
  outF <- file.path(uncSpDir(cfg, sp, "oof"), paste0("oof_", cfg$repLabel, ".rds"))
  if (file.exists(outF)) { message("Out-of-fold inputs exist -- skipping: ", outF); return(invisible(readRDS(outF))) }
  ctx <- uncOofContext(cfg, sp)
  ids <- cfg$reps
  keys <- c(climate = "clim", landscape = "land", habitat = "hab")
  S <- array(NA_real_, c(length(ctx$rows), 3, length(ids)), dimnames = list(NULL, unname(keys), NULL))

  for (scale in names(keys)) {
    mod <- readRDS(file.path(uncSpDir(cfg, sp, "models"), paste0(scale, "_", cfg$repLabel, ".rds")))
    stopifnot(identical(as.integer(mod$ids), as.integer(ids)))
    brtM <- uncMainModel(cfg, sp, scale)
    tbl <- ctx$tabs[[scale]]; fold <- tbl$foldID
    folds <- sort(unique(ctx$foldOf[[scale]]))
    mainFolds <- if (0 %in% ids) {
      f <- file.path(uncMainDir(cfg, sp, scale), paste0(spClean, "_foldModels_", scale, ".rds"))
      if (!file.exists(f)) stop("Fold models of the main BRT missing: ", f)
      readRDS(f)
    }
    t0 <- Sys.time()
    oneRep <- function(j) {
      b <- ids[j]; v <- rep(NA_real_, length(ctx$rows))
      for (k in folds) {
        m <- if (b == 0) mainFolds[[as.character(k)]] else {
          tr <- mod$draws[[j]]; tr <- tr[fold[tr] != k]
          if (sum(tbl$occurrence[tr] == 1) < minPerClass || sum(tbl$occurrence[tr] == 0) < minPerClass) NULL
          else uncWithSeed(stableSeed(c("oofFit", sp, scale, b, k)), uncFitBRT(tbl, mod$predSel, brtM, tr))
        }
        if (is.null(m)) next
        p <- uncOofPredictFold(ctx, scale, m, mod$predSel, k)
        v[p$pos] <- p$value
      }
      v
    }
    # warm the (replicate-independent) covariate/cell cache once, so forked workers do not each rebuild it
    if (scale != "habitat") for (k in folds) { pos <- which(ctx$foldOf[[scale]] == k)
      for (yr in unique(ctx$yrs[pos])) invisible(uncOofCells(ctx, scale, yr, k, mod$predSel, pos[ctx$yrs[pos] == yr])) }
    res <- if (cfg$cores > 1L && length(ids) > 1L) parallel::mclapply(seq_along(ids), oneRep, mc.cores = min(cfg$cores, length(ids))) else lapply(seq_along(ids), oneRep)
    bad <- vapply(res, function(r) !is.numeric(r) || length(r) != length(ctx$rows), logical(1))
    if (any(bad)) stop(scale, ": ", sum(bad), " replicate(s) failed in the out-of-fold step, e.g. replicate ", ids[which(bad)[1]])
    for (j in seq_along(ids)) S[, keys[[scale]], j] <- res[[j]]
    message(scale, ": out-of-fold inputs of ", length(ids), " replicates (", length(folds), " folds each) in ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  }
  out <- list(ids = ids, rows = ctx$rows, S = S)
  uncSaveRDS(out, outF)
  message("Out-of-fold inputs -> ", outF)
  invisible(out)
}
