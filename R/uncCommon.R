#' Shared helpers for the uncertainty workflow (option B: spatial-block bootstrap of the BRTs)
#'
#' Plan, decisions and the exact method: improvements.md item 15, UNCERTAINTY.md (module root) and
#' DECISIONS.md (2026-10-06). Everything here is plain R (no SpaDES session): the tasks read the
#' finished baseline outputs (model_ready tables, main BRTs, covariates) and write ONLY under
#' outputs/<runName>/uncertainty/. Nothing in the baseline run is touched.

#' Label of a replicate run, e.g. reps_001-050 (a later top-up run is e.g. reps_051-100)
uncRepLabel <- function(ids) sprintf("reps_%03d-%03d", min(ids), max(ids))

#' Run `expr` under a fixed RNG seed and restore the caller's RNG state afterwards
#'
#' So that seeding a replicate never changes the random numbers of whatever code runs next.
uncWithSeed <- function(seed, expr) {
  hadSeed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (hadSeed) oldSeed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (hadSeed) assign(".Random.seed", oldSeed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

#' Seed for the attempt-th redraw of a replicate (kept in the int range)
uncAttemptSeed <- function(seed, attempt) as.integer((as.numeric(seed) + attempt * 1000003) %% 2147483629 + 1)

#' Short git hashes of the repo and its submodules (for the replicate log)
uncGitInfo <- function(repoRoot) {   # repoRoot = cfg$codeRoot
  one <- function(dir) tryCatch(system2("git", c("-C", shQuote(dir), "rev-parse", "--short", "HEAD"),
                                        stdout = TRUE, stderr = FALSE)[1], error = function(e) NA_character_,
                                warning = function(w) NA_character_)
  mods <- c("dataPrep_Monitor", "inputs_Monitor", "models_Monitor", "runIndex_Monitor")
  paste(c(paste0("root=", one(repoRoot)),
          paste0(mods, "=", vapply(mods, function(m) one(file.path(repoRoot, "modules", m)), character(1)))),
        collapse = ";")
}

#' Configuration of one uncertainty task, from plain values (the module passes its parameters; tests and the
#' cluster pre-flight pass values read from the shared config files)
#'
#' @param inputRoot Character. The `inputs/` directory (inputPath(sim)).
#' @param outputRoot Character. The run's output directory (outputPath(sim), e.g. outputs/test4).
#' @param species Character. Species Latin names, in the order the cluster arrays index them.
#' @param predictionYears,outYears Integer. Prediction years of the baseline / years to map (default: all prediction years).
#' @param habitatYears Named list (species -> years) or a flat vector: the habitat-scale training years.
#' @param resolutionConfig Named list species -> scale -> resolution (m), or NULL.
#' @param climateResolutionM,habitatResolutionM,landscapeResolutionM Shared default resolutions (m).
#' @param climateWindowLength Integer. Years in the climatology window.
#' @param reps Integer. Replicate ids of THIS run (0 = the main models, no resampling: a built-in consistency check,
#'   left out of every interval; later top-up runs use disjoint ids, e.g. 51:100).
#' @param nBands Integer. Number of horizontal bands the country is cut into.
#' @param repBatch Integer. Replicates held in memory at once.
#' @param cores Integer. Cores for forked prediction.
#' @param blockMult Numeric. Multiplier on the resampling block size (sensitivity test).
#' @param probs Numeric length 2. Interval percentiles (default 5th-95th = 90%).
#' @param tag Character. Non-empty: write to `uncertainty_<tag>` (timing/test runs; never mixes with real replicates).
#' @param baselineYear,currentYear Integer. Index baseline (2005) and the report year (default max habitat year).
#' @param codeRoot Character. Repo root, only used to log git commits.
uncCfgFromParams <- function(inputRoot, outputRoot, species, predictionYears, habitatYears, resolutionConfig = NULL,
                             climateResolutionM, habitatResolutionM, landscapeResolutionM, climateWindowLength,
                             reps, outYears = NULL, nBands = 16L, repBatch = 10L, cores = 1L, blockMult = 1,
                             probs = c(0.05, 0.95), tag = "", baselineYear = 2005L, currentYear = NULL, codeRoot = getwd(), members = NULL) {
  # members: NULL = the BRT-only workflow (default, unchanged). Otherwise the names of the models averaged in every replicate's scale
  # prediction (subset of ALGO_MEMBERS, e.g. c("brt","glm","gam","rf") = the ensemble "ens"); the run then defaults to tag ensTag(members).
  if (!is.null(members) && length(members) && !all(is.na(members))) {
    members <- ALGO_MEMBERS[ALGO_MEMBERS %in% members]
    if (is.null(tag) || is.na(tag) || !nzchar(tag)) tag <- ensTag(members)
  } else members <- NULL
  hyOf <- function(sp) if (is.list(habitatYears)) habitatYears[[sp]] else habitatYears
  allHy <- unlist(if (is.list(habitatYears)) habitatYears else list(habitatYears), use.names = FALSE)
  if (.Platform$OS.type == "windows") cores <- 1L
  list(
    inputRoot = normalizePath(inputRoot, winslash = "/", mustWork = FALSE),
    outputRoot = normalizePath(outputRoot, winslash = "/", mustWork = FALSE),
    codeRoot = codeRoot, tag = if (is.null(tag) || is.na(tag)) "" else tag,
    members = members, fitMembers = setdiff(members, "brt"),
    species = species,
    outYears = if (is.null(outYears)) predictionYears else outYears,
    predictionYears = predictionYears,
    baselineYear = as.integer(baselineYear),
    currentYear = if (is.null(currentYear) || is.na(currentYear)) max(allHy) else as.integer(currentYear),
    habitatYearsOf = hyOf,
    resOf = function(sp) c(climate = resolveResolutionM(sp, "climate", resolutionConfig, climateResolutionM),
                           landscape = resolveResolutionM(sp, "landscape", resolutionConfig, landscapeResolutionM),
                           habitat = resolveResolutionM(sp, "habitat", resolutionConfig, habitatResolutionM)),
    climateWindowLength = climateWindowLength,
    sharedRes = c(europe = climateResolutionM, habitat = habitatResolutionM, landscape = landscapeResolutionM),
    reps = as.integer(reps), repLabel = uncRepLabel(reps),
    nBands = as.integer(nBands), repBatch = as.integer(repBatch), cores = max(1L, as.integer(cores)),
    blockMult = as.numeric(blockMult), probs = as.numeric(probs)
  )
}

# ---- paths ----------------------------------------------------------------------------------------------

.uncModelSuffix <- c(climate = "EU", landscape = "landscape", habitat = "habitat")

#' Root of all uncertainty outputs of this run (BIRDMONITOR_UNC_TAG=timing gives a separate "uncertainty_timing" folder
#' for test runs, so they can never mix with the real replicates)
uncRoot <- function(cfg) file.path(cfg$outputRoot, paste0("uncertainty", if (nzchar(cfg$tag)) paste0("_", cfg$tag) else ""))

#' Per-species output folder (created)
uncSpDir <- function(cfg, sp, ...) {
  d <- file.path(uncRoot(cfg), gsub(" ", "_", sp), ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Folder with the baseline outputs of one scale (main BRT, predictions)
uncMainDir <- function(cfg, sp, scaleKey) file.path(cfg$outputRoot, scaleLabel(cfg$resOf(sp)[[scaleKey]]))

#' The baseline's main BRT of one species and scale
uncMainModel <- function(cfg, sp, scaleKey) {
  f <- file.path(uncMainDir(cfg, sp, scaleKey), paste0(gsub(" ", "_", sp), "_BRT_", .uncModelSuffix[[scaleKey]], ".rds"))
  if (!file.exists(f)) stop("Main BRT missing: ", f, " (has the baseline model array run for this species/scale?)")
  readRDS(f)
}

#' Folder of the model-ready training tables of one scale
#'
#' inputs_Monitor writes every species' table under the SHARED resolution of that scale (e.g. `scale_1` for landscape), even
#' for a species with its own landscape resolution (Lanius collurio 700 m, Buteo buteo 5 km); only model OUTPUTS and
#' covariates live in the per-species resolution folders.
uncModelReadyDir <- function(cfg, scaleKey) {
  shared <- c(climate = "europe", landscape = "landscape", habitat = "habitat")[[scaleKey]]
  file.path(cfg$inputRoot, "model_ready", scaleLabel(cfg$sharedRes[[shared]]))
}

#' The model-ready training table of one species and scale (inputs_Monitor output)
uncTrainingTable <- function(cfg, sp, scaleKey) {
  f <- file.path(uncModelReadyDir(cfg, scaleKey), paste0(gsub(" ", "_", sp), "_inputs.rds"))
  if (!file.exists(f)) stop("Training table missing: ", f)
  readRDS(f)
}

#' Habitat covariate cache file (all layers of one year, written once by uncCovcache())
uncCovcacheFile <- function(cfg, sp, year) {
  lab <- scaleLabel(cfg$resOf(sp)[["habitat"]])
  file.path(cfg$outputRoot, "uncertainty", "covcache", lab, paste0("habitat_", year, ".tif"))   # shared by all tags
}

#' Years for which coarse (climate/landscape) predictions and the habitat cache are needed
uncAllYears <- function(cfg, sp) sort(unique(c(cfg$outYears, cfg$habitatYearsOf(sp))))

# ---- spatial windows -----------------------------------------------------------------------------------

#' The country window on the 200 m reference grid, cut into `nBands` horizontal bands
#'
#' The reference grid is the static solar-radiation raster the baseline meta-model resamples every scale onto;
#' the bands partition ITS cells, so every output cell belongs to exactly one band.
#'
#' @return List: `window` (template SpatRaster of the whole window), `bands` (list of per-band
#'   templates, each with `$ext`), `res`.
uncBands <- function(cfg, sp) {
  firstYear <- uncAllYears(cfg, sp)[1]
  cov1 <- terra::rast(uncCovcacheFile(cfg, sp, firstYear))
  habLab <- scaleLabel(cfg$resOf(sp)[["habitat"]])
  ref <- terra::rast(file.path(cfg$inputRoot, "predictors", "processed", habLab,
                               paste0("solar_radiation_habitat_", habLab, ".tif")))
  win <- terra::crop(ref, terra::ext(cov1), snap = "out")
  win <- terra::rast(win)                       # geometry only
  nr <- terra::nrow(win); r <- terra::res(win)
  edges <- unique(round(seq(0, nr, length.out = cfg$nBands + 1)))
  bands <- lapply(seq_len(length(edges) - 1), function(k) {
    r1 <- edges[k] + 1; r2 <- edges[k + 1]
    ymaxB <- terra::ymax(win) - (r1 - 1) * r[2]; yminB <- terra::ymax(win) - r2 * r[2]
    tmpl <- terra::rast(xmin = terra::xmin(win), xmax = terra::xmax(win), ymin = yminB, ymax = ymaxB,
                        resolution = r, crs = terra::crs(win))
    list(k = k, rows = c(r1, r2), template = tmpl)
  })
  list(window = win, bands = bands, res = r)
}

# ---- bootstrap pieces -----------------------------------------------------------------------------------

#' Resampling-block id of every training row (square grid of `blockSizeM`, anchored at 0)
uncBlockIds <- function(x, y, blockSizeM) paste(floor(x / blockSizeM), floor(y / blockSizeM), sep = "_")

#' Row indices of one block-bootstrap draw
#'
#' Draws as many blocks as there are blocks, with replacement; all rows of a drawn block enter (a block drawn
#' twice contributes its rows twice). A draw lacking presences or absences is redrawn (deterministically).
#'
#' @return Integer row indices; attributes `attempts` and `nBlocks`.
uncRowIndex <- function(blockIds, occurrence, seed, minPerClass = 10, maxAttempts = 50) {
  rowsByBlock <- split(seq_along(blockIds), blockIds)
  nBlocks <- length(rowsByBlock)
  for (attempt in 0:(maxAttempts - 1)) {
    idx <- uncWithSeed(uncAttemptSeed(seed, attempt),
                       unlist(rowsByBlock[sample.int(nBlocks, nBlocks, replace = TRUE)], use.names = FALSE))
    if (sum(occurrence[idx] == 1) >= minPerClass && sum(occurrence[idx] == 0) >= minPerClass)
      return(structure(idx, attempts = attempt + 1L, nBlocks = nBlocks))
  }
  stop("uncRowIndex(): no usable draw after ", maxAttempts, " attempts (too few blocks/records?).")
}

#' Refit a BRT on a bootstrap sample with the main fit's FIXED hyperparameters
#'
#' @param spPa data.frame with `occurrence` and the predictors.
#' @param predSel Character. Predictor columns.
#' @param brtM The main `gbm.step()` model (tree count, learning rate, bag fraction and tree complexity are read from it).
#' @param idx Integer row indices of the draw.
#' @return A `gbm` fitted with `keep.data = FALSE` (small).
uncFitBRT <- function(spPa, predSel, brtM, idx) {
  d <- spPa[idx, c("occurrence", predSel), drop = FALSE]
  gbm::gbm(formula = occurrence ~ ., distribution = "bernoulli", data = d,
           n.trees = brtM$gbm.call$best.trees, shrinkage = brtM$gbm.call$learning.rate,
           bag.fraction = brtM$gbm.call$bag.fraction, interaction.depth = brtM$gbm.call$tree.complexity,
           weights = rep(1, nrow(d)), verbose = FALSE, keep.data = FALSE)
}

#' Predict many gbm models on one data.frame (forked over models) -> matrix cells x models
uncPredictMany <- function(models, df, nTrees, cores = 1L) {
  one <- function(m) gbm::predict.gbm(m, df, n.trees = nTrees, type = "response")
  res <- if (cores > 1L && length(models) > 1L) parallel::mclapply(models, one, mc.cores = min(cores, length(models))) else lapply(models, one)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) stop("uncPredictMany(): prediction failed in ", sum(bad), " of ", length(models), " models.")
  do.call(cbind, res)
}

#' Complete-case cells of a covariate stack as a data.frame (same selection as predictBRTToRaster())
#'
#' @param covStack SpatRaster.
#' @param predictors Character. Predictor columns (may include "x","y", which come from the cell coordinates).
#' @param keepCells Optional integer vector of cell numbers: only these cells are considered.
#' @return List: `df` (complete cells), `idx` (their cell numbers in `covStack`).
uncCells <- function(covStack, predictors, keepCells = NULL) {
  rasterPredictors <- setdiff(predictors, c("x", "y"))
  predRast <- covStack[[rasterPredictors]]
  if (is.null(keepCells)) {
    predDf <- as.data.frame(predRast, xy = TRUE, na.rm = FALSE)
    cellIds <- seq_len(nrow(predDf))
  } else {
    cellIds <- sort(unique(keepCells[!is.na(keepCells)]))
    predDf <- cbind(as.data.frame(terra::xyFromCell(predRast, cellIds)), terra::extract(predRast, cellIds))
  }
  complete <- stats::complete.cases(predDf[, predictors, drop = FALSE])
  list(df = predDf[complete, predictors, drop = FALSE], idx = cellIds[complete])
}

#' Row/column-neighbourhood (+-`reach` cells) of a set of cells, as cell numbers
uncNeighbourCells <- function(r, cells, reach = 2L) {
  cells <- cells[!is.na(cells)]
  if (!length(cells)) return(integer(0))
  rc <- terra::rowColFromCell(r, cells)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  out <- lapply(-reach:reach, function(dr) lapply(-reach:reach, function(dc) {
    rr <- rc[, 1] + dr; cc <- rc[, 2] + dc
    ok <- rr >= 1 & rr <= nr & cc >= 1 & cc <= nc
    terra::cellFromRowCol(r, rr[ok], cc[ok])
  }))
  sort(unique(unlist(out)))
}

#' Index of the nearest reference point for every query point (planar distance, plain R: no sf needed)
#'
#' Same result as `sf::st_nearest_feature()` for projected coordinates (EPSG:3035 metres); ties go to the first reference point.
#' Chunked, so memory stays at a few hundred MB however many points there are.
#'
#' @param q Numeric matrix (n x 2), query x/y.
#' @param ref Numeric matrix (m x 2), reference x/y.
#' @return Integer vector of length n.
uncNearestIdx <- function(q, ref, chunk = 500L) {
  q <- as.matrix(q); ref <- as.matrix(ref)
  out <- integer(nrow(q))
  for (s in seq(1L, nrow(q), by = chunk)) {
    i <- s:min(s + chunk - 1L, nrow(q))
    d <- outer(q[i, 1], ref[, 1], "-")^2 + outer(q[i, 2], ref[, 2], "-")^2
    out[i] <- max.col(-d, ties.method = "first")
  }
  out
}

# ---- ensemble members inside the replicates --------------------------------------------------------------------------

.uncMemberCache <- new.env(parent = emptyenv())

#' The replicate fits of one non-BRT member at one scale (models/<scale>_<member>_<run>.rds), cached per process
uncMemberFits <- function(cfg, sp, scale, member) {
  f <- file.path(uncSpDir(cfg, sp, "models"), paste0(scale, "_", member, "_", cfg$repLabel, ".rds"))
  if (is.null(.uncMemberCache[[f]])) {
    if (!file.exists(f)) stop("Replicate fits of member '", member, "' missing: ", f, " (run the 'fit' step with the members set)")
    .uncMemberCache[[f]] <- readRDS(f)
  }
  .uncMemberCache[[f]]
}

#' Predict a list of algorithm fits on one data.frame -> matrix rows x fits
uncPredictFits <- function(fits, df, cores = 1L) {
  one <- function(ft) algoPredict(ft, df)
  res <- if (cores > 1L && length(fits) > 1L) parallel::mclapply(fits, one, mc.cores = min(cores, length(fits))) else lapply(fits, one)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) stop("uncPredictFits(): prediction failed in ", sum(bad), " of ", length(fits), " fits.")
  do.call(cbind, res)
}

#' Scale prediction of the replicates `pos`: the BRT alone (default) or the MEAN of the run's members
#'
#' @param mod The BRT bundle of the scale (`models/<scale>_<run>.rds`).
#' @return Matrix rows of `df` x replicates `pos`. With `cfg$members` NULL this is exactly `uncPredictMany()` of the BRTs.
uncPredictScale <- function(cfg, sp, scale, mod, pos, df) {
  useBrt <- is.null(cfg$members) || "brt" %in% cfg$members
  parts <- list()
  if (useBrt) parts$brt <- uncPredictMany(mod$models[pos], df, mod$nTrees, cfg$cores)
  for (m in cfg$fitMembers) {
    mf <- uncMemberFits(cfg, sp, scale, m)
    if (!identical(as.integer(mf$ids), as.integer(mod$ids))) stop("Replicate ids of member '", m, "' differ from the BRT's at scale ", scale)
    parts[[m]] <- uncPredictFits(mf$fits[pos], df, cfg$cores)
  }
  if (length(parts) == 1L) parts[[1]] else Reduce(`+`, parts) / length(parts)
}
