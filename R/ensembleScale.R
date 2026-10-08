#' The ensemble of one scale of one species: arithmetic mean of ANY chosen set of members
#'
#' Wiedenroth et al. average the algorithms of a scale (arithmetic mean) into that scale's prediction, which the meta-model then combines with the
#' other two scales. Here the members are chosen by name (`members`, a subset of `ALGO_MEMBERS` = brt, glm, gam, rf, nn): all of them, a single
#' one (then the "ensemble" is just that model), or any combination, so "one out, two out" experiments only need another `members` vector.
#' Averaging maps is cheap, so many ensembles can be compared without refitting anything: the members are fitted ONCE by `algoScaleRun()`.
#'
#' Output (next to the members' files, `uncMainDir()`), `<sp>` = species with underscores, `<tag>` = `ensTag(members)` ("ens" for the default
#' brt+glm+gam+rf, otherwise e.g. "ens_brt-gam"), so ensembles of different members never overwrite each other:
#'   `<sp>_perf_<tag>_<scale>.rds`        block-CV performance of the ensemble (mean of the members' out-of-fold predictions)
#'   `<sp>_cvpred_<tag>_<scale>.rds`      that out-of-fold prediction
#'   `<sp>_oofhab_<tag>_<scale>.rds`      out-of-fold ensemble prediction at every habitat record (honest meta-model input)
#'   `<sp>_perf_members_<tag>_<scale>.rds` block-CV performance of every member and of the ensemble, one table
#'   `<sp>_pred_<tag>_<scale>_<year>.tif` ensemble map (`mean_prob`, `binary`); NA where ANY member is NA, so a cell never changes mixture
#'   `<sp>_sd_<tag>_<scale>_<year>.tif`   standard deviation ACROSS the members, cell by cell (disagreement between algorithms)
#' The BRT member needs only the fold models of the BRT workflow (rebuilt automatically by `brtMemberRun()`).
ensembleScaleRun <- function(cfg, sp, scale, members = c("brt", "glm", "gam", "rf"), years = uncAllYears(cfg, sp)) {
  members <- ALGO_MEMBERS[ALGO_MEMBERS %in% members]
  tag <- ensTag(members)
  d <- uncMainDir(cfg, sp, scale); spClean <- gsub(" ", "_", sp)
  msg <- paste0(sp, " | ", scale, " | ", tag)
  tbl <- uncTrainingTable(cfg, sp, scale)
  if ("brt" %in% members) brtMemberRun(cfg, sp, scale)
  fPerf <- ensFile(cfg, sp, scale, tag, "perf"); fCv <- ensFile(cfg, sp, scale, tag, "cvpred"); fOof <- ensFile(cfg, sp, scale, tag, "oofhab")
  need <- unlist(lapply(members, function(a) c(ensFile(cfg, sp, scale, a, "cvpred"), ensFile(cfg, sp, scale, a, "oofhab"))))
  if (!all(file.exists(need))) stop(msg, ": missing member results (run algoScaleRun for them): ", paste(basename(need[!file.exists(need)]), collapse = ", "))
  inputsStamp <- max(file.mtime(need))

  # 1. performance and out-of-fold inputs
  if (!(file.exists(fPerf) && file.exists(fCv) && file.exists(fOof) && file.mtime(fPerf) >= inputsStamp)) {
    cvAll <- setNames(lapply(members, function(a) readRDS(ensFile(cfg, sp, scale, a, "cvpred"))), members)
    oofAll <- setNames(lapply(members, function(a) readRDS(ensFile(cfg, sp, scale, a, "oofhab"))), members)
    cvEns <- algoEnsembleMean(cvAll); oofEns <- algoEnsembleMean(oofAll)
    ok <- !is.na(cvEns)
    perf <- evalSDM(tbl$occurrence[ok], cvEns[ok])
    perfAll <- do.call(rbind, lapply(names(cvAll), function(n) { o <- !is.na(cvAll[[n]]); cbind(member = n, evalSDM(tbl$occurrence[o], cvAll[[n]][o])) }))
    saveRDS(cvEns, fCv); saveRDS(oofEns, fOof); saveRDS(perf, fPerf)
    saveRDS(rbind(perfAll, cbind(member = tag, perf)), ensFile(cfg, sp, scale, tag, "members"))
    message(msg, ": ensemble block-CV AUC ", round(perf$AUC, 3), " | TSS ", round(perf$TSS, 3), " | D2 ", round(perf$D2, 3),
            " (members: ", paste(sprintf("%s %.3f", perfAll$member, perfAll$AUC), collapse = ", "), ")")
  }
  perf <- readRDS(fPerf)

  # 2. maps
  for (yr in years) {
    out <- ensFile(cfg, sp, scale, tag, "pred", yr); outSd <- ensFile(cfg, sp, scale, tag, "sd", yr)
    mf <- vapply(members, function(a) ensFile(cfg, sp, scale, a, "pred", yr), "")
    if (!all(file.exists(mf))) { message(msg, " ", yr, ": a member map is missing -- skipped"); next }
    if (isValidPredictionRaster(out) && file.exists(outSd) && file.mtime(out) >= max(file.mtime(mf)) && file.mtime(out) >= file.mtime(fPerf)) next
    layers <- lapply(mf, function(f) terra::rast(f)[["mean_prob"]])
    m <- terra::mean(terra::rast(layers))        # na.rm = FALSE: NA where any member is NA
    names(m) <- "mean_prob"
    b <- m >= perf$thresh; names(b) <- "binary"
    terra::writeRaster(combineTwoLayerRaster(m, b, c("mean_prob", "binary")), out, overwrite = TRUE)
    # disagreement between the members (standard deviation across them, cell by cell), as in Wiedenroth et al.'s 04e: structural uncertainty
    if (length(members) > 1L) { sdr <- terra::stdev(terra::rast(layers), pop = FALSE); names(sdr) <- "sd_members"; terra::writeRaster(sdr, outSd, overwrite = TRUE); rm(sdr)
    } else terra::writeRaster(terra::rast(layers)[[1]] * 0, outSd, overwrite = TRUE)
    message(msg, " ", yr, ": saved ", basename(out))
    rm(layers, m, b); invisible(gc())
  }
  invisible(TRUE)
}
