#' Steps that run once per species BEFORE any map is made: fit the replicate BRTs, cache the
#' habitat covariates, and predict the coarse (climate, landscape) scales.

#' Resampling block size per scale, from the SAME rule as the spatial cross-validation
#'
#' `cv_spatial_autocor()` range, capped and floored exactly as inputs_Monitor does (climate 200-1500 km;
#' habitat and landscape floor = 2 x resolution, cap 200 km). The autocorrelation step is run once per species
#' and cached, so later top-up runs resample with identical blocks.
#'
#' @return Named numeric vector (m): climate, landscape, habitat.
uncBlockSizes <- function(cfg, sp) {
  f <- file.path(uncSpDir(cfg, sp), "blocksize.rds")
  if (file.exists(f)) return(readRDS(f))
  res <- cfg$resOf(sp)
  limits <- list(climate = c(min = 200000, max = 1500000),
                 landscape = c(min = 2 * res[["landscape"]], max = 200000),
                 habitat = c(min = 2 * res[["habitat"]], max = 200000))
  out <- vapply(names(limits), function(s) {
    spPa <- uncTrainingTable(cfg, sp, s)
    sfOcc <- sf::st_as_sf(spPa, coords = c("x", "y"), crs = 3035)
    determineBlockSize(sfOcc, maxBlockSizeM = limits[[s]][["max"]], minBlockSizeM = limits[[s]][["min"]])
  }, numeric(1))
  uncSaveRDS(out, f)
  out
}

#' Fit the replicate BRTs of one species (all three scales) for the replicate ids of this run
#'
#' Per scale and replicate: draw spatial blocks with replacement (seeded by species/scale/replicate id), refit the
#' BRT with the main model's fixed hyperparameters. Replicate 0, if requested, is the main model itself
#' (consistency check). Writes `models/<scale>_<run>.rds` and appends one row per fit to `replicate_log.csv`.
uncFitSpecies <- function(cfg, sp) {
  spClean <- gsub(" ", "_", sp)
  outDir <- uncSpDir(cfg, sp, "models")
  blockSizes <- uncBlockSizes(cfg, sp)
  message("Resampling block sizes (km): ", paste(names(blockSizes), round(blockSizes / 1000, 1), collapse = ", "),
          if (cfg$blockMult != 1) paste0("  x multiplier ", cfg$blockMult) else "")
  git <- uncGitInfo(cfg$codeRoot)
  allLogs <- list()

  for (scaleKey in c("climate", "landscape", "habitat")) {
    f <- file.path(outDir, paste0(scaleKey, "_", cfg$repLabel, ".rds"))
    if (file.exists(f)) { message(scaleKey, ": ", basename(f), " exists -- skipping"); next }
    spPa <- uncTrainingTable(cfg, sp, scaleKey)
    brtM <- uncMainModel(cfg, sp, scaleKey)
    predSel <- brtM$var.names
    stopifnot(all(c("x", "y", "occurrence") %in% names(spPa)), all(predSel %in% names(spPa)))
    blockSizeM <- blockSizes[[scaleKey]] * cfg$blockMult
    blockIds <- uncBlockIds(spPa$x, spPa$y, blockSizeM)
    message(scaleKey, ": ", nrow(spPa), " records, ", length(unique(blockIds)), " resampling blocks, ",
            length(cfg$reps), " replicates, ", length(predSel), " predictors")

    fitOne <- function(b) {
      t0 <- Sys.time()
      seedDraw <- stableSeed(c("bootDraw", sp, scaleKey, b)); seedFit <- stableSeed(c("bootFit", sp, scaleKey, b))
      if (b == 0) {
        idx <- structure(seq_len(nrow(spPa)), attempts = 0L, nBlocks = length(unique(blockIds))); model <- brtM
      } else {
        idx <- uncRowIndex(blockIds, spPa$occurrence, seedDraw)
        model <- uncWithSeed(seedFit, uncFitBRT(spPa, predSel, brtM, idx))
      }
      list(model = model, idx = as.integer(idx),
           log = data.frame(species = sp, scale = scaleKey, replicate = b, runLabel = cfg$repLabel,
                            seedDraw = if (b == 0) NA else seedDraw, seedFit = if (b == 0) NA else seedFit,
                            attempts = attr(idx, "attempts"), blocksTotal = attr(idx, "nBlocks"),
                            blockSizeM = blockSizeM, rowsDrawn = length(idx),
                            presDrawn = sum(spPa$occurrence[idx] == 1), absDrawn = sum(spPa$occurrence[idx] == 0),
                            bestTrees = brtM$gbm.call$best.trees, learningRate = brtM$gbm.call$learning.rate,
                            bagFraction = brtM$gbm.call$bag.fraction, treeComplexity = brtM$gbm.call$tree.complexity,
                            fitSeconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2),
                            git = git, rVersion = paste0(R.version$major, ".", R.version$minor),
                            time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stringsAsFactors = FALSE))
    }
    fits <- if (cfg$cores > 1L) parallel::mclapply(cfg$reps, fitOne, mc.cores = min(cfg$cores, length(cfg$reps))) else lapply(cfg$reps, fitOne)
    bad <- vapply(fits, function(x) !is.list(x) || is.null(x$model), logical(1))
    if (any(bad)) stop(scaleKey, ": ", sum(bad), " replicate fit(s) failed, e.g. replicate ", cfg$reps[which(bad)[1]])
    uncSaveRDS(list(ids = cfg$reps, models = lapply(fits, `[[`, "model"), draws = lapply(fits, `[[`, "idx"),
                 predSel = predSel, nTrees = brtM$gbm.call$best.trees, blockSizeM = blockSizeM, scaleKey = scaleKey),
            f)
    allLogs <- c(allLogs, lapply(fits, `[[`, "log"))
    message(scaleKey, ": saved ", basename(f), " (", round(file.size(f) / 1e6), " MB)")
  }

  if (length(allLogs)) {
    logDf <- do.call(rbind, allLogs)
    logFile <- file.path(uncSpDir(cfg, sp), "replicate_log.csv")
    utils::write.table(logDf, logFile, sep = ",", row.names = FALSE, col.names = !file.exists(logFile),
                       append = file.exists(logFile), qmethod = "double")
    message("Replicate log -> ", logFile)
  }
  infoFile <- file.path(uncSpDir(cfg, sp), paste0("run_info_", cfg$repLabel, ".txt"))
  if (!file.exists(infoFile)) writeLines(c(paste("git:", git), capture.output(print(utils::sessionInfo()))), infoFile)
  invisible(TRUE)
}

#' Write the habitat covariate stacks of all needed years once (shared by all species of that resolution)
#'
#' The baseline's own loader builds each year's stack on the fly (about a minute and several GB per year); the
#' band tasks read windows of these cached files instead.
uncCovcache <- function(cfg) {
  labels <- vapply(cfg$species, function(s) scaleLabel(cfg$resOf(s)[["habitat"]]), character(1))
  for (lab in unique(labels)) {
    spOf <- cfg$species[labels == lab]
    years <- sort(unique(unlist(lapply(spOf, function(s) uncAllYears(cfg, s)))))
    habitatDir <- file.path(cfg$inputRoot, "predictors", "processed", lab)
    for (yr in years) {
      f <- uncCovcacheFile(cfg, spOf[1], yr)
      ok <- file.exists(f) && tryCatch(terra::nlyr(terra::rast(f)) > 0, error = function(e) FALSE)
      if (ok) { message(lab, " ", yr, ": cached"); next }
      dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
      cov <- loadHabitatCovariates(yr, habitatDir)
      if (is.null(cov)) { warning(lab, " ", yr, ": covariates unavailable -- not cached"); next }
      part <- sub("[.]tif$", ".part.tif", f)
      terra::writeRaster(cov, part, datatype = "FLT4S", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
      file.rename(part, f)
      message(lab, " ", yr, ": cached -> ", f, " (", round(file.size(f) / 1e6), " MB)")
      rm(cov); invisible(gc())
    }
  }
  invisible(TRUE)
}

#' Predict the replicate climate and landscape BRTs for every needed year
#'
#' Cheap scales (50 km and 1 km grids). One multi-layer GeoTIFF per scale and year, one layer per replicate
#' (`rep_<id>`), on the scale's own grid; the band tasks resample them onto the 200 m reference grid exactly
#' like the baseline's `loadSuitability()`.
uncCoarseSpecies <- function(cfg, sp) {
  years <- uncAllYears(cfg, sp)
  res <- cfg$resOf(sp)
  outDir <- uncSpDir(cfg, sp, "coarse", cfg$repLabel)
  for (scaleKey in c("climate", "landscape")) {
    mod <- readRDS(file.path(uncSpDir(cfg, sp, "models"), paste0(scaleKey, "_", cfg$repLabel, ".rds")))
    lab <- scaleLabel(res[[scaleKey]])
    procDir <- file.path(cfg$inputRoot, "predictors", "processed", lab)
    for (yr in years) {
      f <- file.path(outDir, sprintf("%s_%d.tif", scaleKey, yr))
      if (file.exists(f) && tryCatch(terra::nlyr(terra::rast(f)) == length(mod$ids), error = function(e) FALSE)) next
      cov <- if (scaleKey == "climate") {
        bf <- file.path(procDir, paste0("bioclim_", yr - (cfg$climateWindowLength - 1), "-", yr, "_", lab, ".tif"))
        if (!file.exists(bf)) { message(scaleKey, " ", yr, ": climatology missing -- skipping"); next }
        terra::rast(bf)
      } else {
        cv <- loadCovariates(yr, procDir, NULL)
        if (is.null(cv)) { message(scaleKey, " ", yr, ": covariates unavailable -- skipping"); next }
        names(cv) <- gsub("_\\d{4}$", "", names(cv)); cv
      }
      missingPreds <- setdiff(mod$predSel, c(names(cov), "x", "y"))
      if (length(missingPreds)) { warning(scaleKey, " ", yr, ": missing predictors ", paste(missingPreds, collapse = ", ")); next }
      cells <- uncCells(cov, mod$predSel)
      P <- uncPredictMany(mod$models, cells$df, mod$nTrees, cfg$cores)
      vals <- matrix(NA_real_, terra::ncell(cov), ncol(P)); vals[cells$idx, ] <- P
      r <- terra::rast(cov[[1]], nlyrs = ncol(P)); terra::values(r) <- vals
      names(r) <- paste0("rep_", mod$ids)
      terra::writeRaster(r, f, datatype = "FLT4S", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
      message(scaleKey, " ", yr, ": ", length(mod$ids), " replicates -> ", basename(f))
    }
  }
  invisible(TRUE)
}
