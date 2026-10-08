#' Alternative algorithms (GLM, GAM, random forest) at one scale of one species: fit, block cross-validation, predictions for all years
#'
#' One call = one (species, scale, algorithm) job of the ensemble workflow (Wiedenroth et al.: every scale is fitted with several algorithms,
#' whose predictions are averaged; improvements.md item 7). The BRT of each scale is the existing one (`modelEurope()` etc.) and is not
#' refitted here. Everything uses the SAME predictors and spatial folds as the BRT of that scale, so the algorithms are comparable.
#'
#' Files written next to the BRT outputs of that scale (`uncMainDir()`), `<sp>` = species with underscores, `<scale>` = climate | landscape |
#' habitat, `<algo>` = glm | gam | rf:
#'   `<sp>_<algo>_<scale>.rds`            the fitted model (list from `algoFit()`)
#'   `<sp>_perf_<algo>_<scale>.rds`       block-CV performance (`evalSDM()` row, with the threshold)
#'   `<sp>_cvpred_<algo>_<scale>.rds`     out-of-fold prediction of every training row
#'   `<sp>_oofhab_<algo>_<scale>.rds`     out-of-fold prediction at every habitat record (the honest input of the meta-model; see uncOof.R)
#'   `<sp>_pred_<algo>_<scale>_<year>.tif`  prediction map (`mean_prob`, `binary`) for every year
#' Existing complete files are reused (a prediction map only if it is not older than its climatology / model).

#' Path of one ensemble file (see above); `kind` = model | perf | cvpred | oofhab | pred (pred needs `yr`)
ensFile <- function(cfg, sp, scale, algo, kind, yr = NULL) {
  spClean <- gsub(" ", "_", sp); d <- uncMainDir(cfg, sp, scale)
  # the BRT's maps are those of the existing BRT workflow (its model, folds and performance files keep their own names)
  if (algo == "brt" && kind == "pred")
    return(file.path(d, sprintf("%s_pred_%s_%d.tif", spClean, c(climate = "EU", landscape = "landscape", habitat = "habitat")[[scale]], yr)))
  switch(kind,
    members = file.path(d, sprintf("%s_perf_members_%s_%s.rds", spClean, algo, scale)),
    sd     = file.path(d, sprintf("%s_sd_%s_%s_%d.tif", spClean, algo, scale, yr)),
    model  = file.path(d, sprintf("%s_%s_%s.rds", spClean, algo, scale)),
    perf   = file.path(d, sprintf("%s_perf_%s_%s.rds", spClean, algo, scale)),
    cvpred = file.path(d, sprintf("%s_cvpred_%s_%s.rds", spClean, algo, scale)),
    oofhab = file.path(d, sprintf("%s_oofhab_%s_%s.rds", spClean, algo, scale)),
    pred   = file.path(d, sprintf("%s_pred_%s_%s_%d.tif", spClean, algo, scale, yr)),
    stop("ensFile(): unknown kind ", kind))
}

#' Predictors of a scale = those of its BRT (same set for every algorithm)
ensPredictors <- function(cfg, sp, scale) uncMainModel(cfg, sp, scale)$var.names

#' Climatology file of a climate target year (NA for the other scales)
ensClimateFile <- function(cfg, sp, yr) {
  lab <- scaleLabel(cfg$resOf(sp)[["climate"]])
  file.path(cfg$inputRoot, "predictors", "processed", lab, paste0("bioclim_", yr - (cfg$climateWindowLength - 1), "-", yr, "_", lab, ".tif"))
}

#' Covariate stack of one scale and year (NULL if unavailable)
ensCovStack <- function(cfg, sp, scale, yr) {
  lab <- scaleLabel(cfg$resOf(sp)[[scale]]); dir <- file.path(cfg$inputRoot, "predictors", "processed", lab)
  switch(scale,
    climate = { f <- ensClimateFile(cfg, sp, yr); if (file.exists(f)) terra::rast(f) else NULL },
    landscape = { cv <- loadCovariates(yr, dir, NULL); if (!is.null(cv)) names(cv) <- gsub("_[0-9]{4}$", "", names(cv)); cv },
    habitat = { f <- uncCovcacheFile(cfg, sp, yr)
      if (file.exists(f)) terra::rast(f) else { warning("habitat covariate cache missing: ", f, " (run the covcache step first)"); NULL } })
}

#' The BRT as a member of the ensemble: its out-of-fold predictions, rebuilt from the fold models the BRT workflow saved (no refitting)
#'
#' Writes `<sp>_cvpred_brt_<scale>.rds` (every training row) and `<sp>_oofhab_brt_<scale>.rds` (every habitat record), the same files the other
#' members have, so the ensemble step treats all members alike. The BRT maps are the existing ones (`ensFile(..., "brt", "pred", yr)`).
brtMemberRun <- function(cfg, sp, scale) {
  spClean <- gsub(" ", "_", sp); d <- uncMainDir(cfg, sp, scale)
  fCv <- ensFile(cfg, sp, scale, "brt", "cvpred"); fOof <- ensFile(cfg, sp, scale, "brt", "oofhab"); ff <- file.path(d, paste0(spClean, "_foldModels_", scale, ".rds"))
  if (!file.exists(ff)) stop(sp, " | ", scale, " | brt: fold models of the BRT missing: ", ff)
  if (file.exists(fCv) && file.exists(fOof) && file.mtime(fCv) >= file.mtime(ff)) return(invisible(TRUE))
  tbl <- uncTrainingTable(cfg, sp, scale); predSel <- ensPredictors(cfg, sp, scale); fm <- readRDS(ff)
  nT <- function(m) if (!is.null(m$gbm.call$best.trees)) m$gbm.call$best.trees else m$n.trees
  cvB <- rep(NA_real_, nrow(tbl))
  for (k in names(fm)) { te <- which(tbl$foldID == as.integer(k)); if (length(te)) cvB[te] <- gbm::predict.gbm(fm[[k]], tbl[te, predSel, drop = FALSE], n.trees = nT(fm[[k]]), type = "response") }
  ctx <- uncOofContext(cfg, sp); oofB <- rep(NA_real_, length(ctx$rows))
  for (k in names(fm)) { p <- uncOofPredictFold(ctx, scale, fm[[k]], predSel, as.integer(k)); oofB[p$pos] <- p$value }
  ok <- !is.na(cvB)
  saveRDS(cvB, fCv); saveRDS(oofB, fOof); saveRDS(evalSDM(tbl$occurrence[ok], cvB[ok]), ensFile(cfg, sp, scale, "brt", "perf"))
  message(sp, " | ", scale, " | brt: out-of-fold predictions rebuilt from the fold models")
  invisible(TRUE)
}

algoScaleRun <- function(cfg, sp, scale, algo, threads = 1L, years = uncAllYears(cfg, sp)) {
  stopifnot(scale %in% c("climate", "landscape", "habitat"), algo %in% ALGO_MEMBERS)
  if (algo == "brt") return(brtMemberRun(cfg, sp, scale))
  tag <- paste0(sp, " | ", scale, " | ", algo)
  tbl <- uncTrainingTable(cfg, sp, scale); predSel <- ensPredictors(cfg, sp, scale)
  stopifnot(all(predSel %in% names(tbl)), "foldID" %in% names(tbl))
  fMod <- ensFile(cfg, sp, scale, algo, "model"); fPerf <- ensFile(cfg, sp, scale, algo, "perf")
  fCv <- ensFile(cfg, sp, scale, algo, "cvpred"); fOof <- ensFile(cfg, sp, scale, algo, "oofhab")

  # 1. final model
  if (file.exists(fMod)) { fit <- readRDS(fMod); fit$threads <- threads; message(tag, ": model exists -- reused") } else {
    t0 <- Sys.time()
    fit <- algoFit(algo, tbl, predSel, seed = stableSeed(c("algoFit", sp, scale, algo)), threads = threads)
    saveRDS(fit, fMod)
    message(tag, ": fitted in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min (", length(fit$predSel), " predictors", if (nzchar(fit$note)) paste0("; ", fit$note), ")")
  }

  # 2. block CV and the out-of-fold prediction at the habitat records (one pass; fold models are not kept)
  if (!(file.exists(fPerf) && file.exists(fCv) && file.exists(fOof))) {
    t0 <- Sys.time()
    ctx <- uncOofContext(cfg, sp)
    cvPred <- rep(NA_real_, nrow(tbl)); oofHab <- rep(NA_real_, length(ctx$rows))
    for (k in sort(unique(tbl$foldID))) {
      tr <- which(tbl$foldID != k); te <- which(tbl$foldID == k)
      if (!length(tr) || !length(te)) next
      fk <- tryCatch(algoRefit(fit, tbl[tr, , drop = FALSE], seed = stableSeed(c("cvAlgo", sp, scale, algo, k))), error = function(e) NULL)
      if (is.null(fk)) { warning(tag, ": fold ", k, " could not be fitted"); next }
      fk$threads <- threads
      cvPred[te] <- algoPredict(fk, tbl[te, , drop = FALSE])
      p <- uncOofPredictFold(ctx, scale, NULL, predSel, k, predictFun = function(df) algoPredict(fk, df))
      oofHab[p$pos] <- p$value
    }
    perf <- evalSDM(tbl$occurrence[!is.na(cvPred)], cvPred[!is.na(cvPred)])
    saveRDS(cvPred, fCv); saveRDS(oofHab, fOof); saveRDS(perf, fPerf)
    message(tag, ": block-CV AUC ", round(perf$AUC, 3), " | TSS ", round(perf$TSS, 3), " | D2 ", round(perf$D2, 3), " (",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
  }
  perf <- readRDS(fPerf)

  # 3. predictions
  for (yr in years) {
    out <- ensFile(cfg, sp, scale, algo, "pred", yr)
    inputs <- c(fMod, if (scale == "climate") ensClimateFile(cfg, sp, yr))
    if (isValidPredictionRaster(out) && all(file.mtime(out) >= file.mtime(inputs[file.exists(inputs)]))) next
    cov <- ensCovStack(cfg, sp, scale, yr)
    if (is.null(cov)) { message(tag, " ", yr, ": covariates unavailable -- skipped"); next }
    miss <- setdiff(predSel, c(names(cov), "x", "y"))
    if (length(miss)) { warning(tag, " ", yr, ": missing predictors ", paste(miss, collapse = ", ")); next }
    t0 <- Sys.time()
    r <- predictModelToRaster(cov, predSel, function(df) algoPredict(fit, df), thresh = perf$thresh)
    terra::writeRaster(r, out, overwrite = TRUE)
    message(tag, " ", yr, ": saved ", basename(out), " (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
    rm(cov, r); invisible(gc())
  }
  invisible(TRUE)
}
