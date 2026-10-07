#' Summaries across replicates: per-year maps, change maps, per-pixel trend, area-mean series, indices
#'
#' All summaries are computed INSIDE the replicate first (change = late minus early, trend = slope, area mean,
#' index) and only then summarised across replicates, so spatial AND temporal uncertainty carry into the derived
#' layers. Intervals are percentile intervals of the fitted probability surface (default 5th-95th = 90%).

#' Mean, SD, lower, upper percentile and width across the columns of a (cells x replicates) matrix without NAs
uncSummarizeMatrix <- function(m, probs = c(0.05, 0.95)) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    q <- matrixStats::rowQuantiles(m, probs = probs, na.rm = FALSE, useNames = FALSE)
    cbind(mean = rowMeans(m), sd = matrixStats::rowSds(m), lwr = q[, 1], upr = q[, 2], width = q[, 2] - q[, 1])
  } else {                                                # slow fallback (tests, tiny inputs)
    q <- t(apply(m, 1, stats::quantile, probs = probs, names = FALSE))
    cbind(mean = rowMeans(m), sd = apply(m, 1, stats::sd), lwr = q[, 1], upr = q[, 2], width = q[, 2] - q[, 1])
  }
}

#' Summary of a change (or trend) matrix: mean, SD, lower, upper, width, share of replicates decreasing / increasing
#' Layer names: `<prefix>Mean/Sd/Lwr/Upr/Width` (prefix `delta` for a change, `slope` for the per-decade trend) and the two shares.
uncSummarizeChange <- function(d, probs = c(0.05, 0.95), prefix = "delta") {
  s <- uncSummarizeMatrix(d, probs)
  out <- cbind(s, shareDecrease = rowMeans(d < 0), shareIncrease = rowMeans(d > 0))
  colnames(out) <- c(paste0(prefix, c("Mean", "Sd", "Lwr", "Upr", "Width")), "shareDecrease", "shareIncrease")
  out
}

#' Replicate-run folders that exist for a species (reps_001-050, reps_051-100, ...)
uncRunLabels <- function(cfg, sp) {
  d <- file.path(uncSpDir(cfg, sp), "pred")
  if (!dir.exists(d)) return(character(0))
  sort(grep("^reps_", list.dirs(d, recursive = FALSE, full.names = FALSE), value = TRUE))
}

#' Read one band and year across all replicate runs of a species
#'
#' @param includeZero Keep replicate 0 (the main models; consistency check)? Default FALSE.
#' @return List(idx, P (valid cells x replicates), ids), or NULL if no run has the file.
uncReadBandYear <- function(cfg, sp, year, bandK, includeZero = FALSE) {
  parts <- list()
  for (lab in uncRunLabels(cfg, sp)) {
    f <- file.path(uncSpDir(cfg, sp), "pred", lab, sprintf("%d_band%02d.rds", year, bandK))
    if (file.exists(f)) parts[[lab]] <- uncReadInt16(f)
  }
  if (!length(parts)) return(NULL)
  idx <- parts[[1]]$idx
  if (!all(vapply(parts, function(p) identical(p$idx, idx), logical(1)))) {
    idx <- Reduce(union, lapply(parts, `[[`, "idx")); idx <- sort(idx)
    parts <- lapply(parts, function(p) { m <- matrix(NA_real_, length(idx), ncol(p$P)); m[match(p$idx, idx), ] <- p$P; p$P <- m; p$idx <- idx; p })
  }
  P <- do.call(cbind, lapply(parts, `[[`, "P")); ids <- unlist(lapply(parts, `[[`, "ids"))
  if (anyDuplicated(ids)) stop(sp, ": replicate ids overlap between runs (", paste(unique(ids[duplicated(ids)]), collapse = ","),
                               ") -- two runs must use disjoint ids, e.g. 1:50 and 51:100.")
  if (!includeZero) { keep <- ids != 0; P <- P[, keep, drop = FALSE]; ids <- ids[keep] }
  list(idx = idx, P = P, ids = ids)
}

#' Write a band piece: values for the valid cells of a band template, NA elsewhere
uncWritePiece <- function(template, idx, mat, file) {
  r <- terra::rast(template, nlyrs = ncol(mat))
  vals <- matrix(NA_real_, terra::ncell(template), ncol(mat)); vals[idx, ] <- mat
  terra::values(r) <- vals; names(r) <- colnames(mat)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  part <- sub("[.]tif$", ".part.tif", file)                  # written under a temporary name, renamed when complete
  terra::writeRaster(r, part, datatype = "FLT4S", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  file.rename(part, file)
  invisible(file)
}

#' The change comparisons of the annual report: baseline, 5 years ago, last year (each against the current year)
uncComparisons <- function(cfg) list(vsBaseline = cfg$baselineYear, vs5YearsAgo = cfg$currentYear - 5L, vsLastYear = cfg$currentYear - 1L)

#' Replicate ids (without 0) that exist for a species across all runs, from the stats files of one band
uncRunIds <- function(cfg, sp, bandK = 1L) {
  ids <- integer(0)
  for (lab in uncRunLabels(cfg, sp)) {
    fs <- list.files(file.path(uncSpDir(cfg, sp), "pred", lab), pattern = sprintf("^stats_.*_band%02d[.]rds$", bandK), full.names = TRUE)
    if (length(fs)) ids <- c(ids, readRDS(fs[1])$ids)
  }
  sort(unique(ids[ids != 0]))
}

#' Delete a band's old pieces when the set of replicates changed (e.g. after adding replicates), then record the stamp
uncPieceStamp <- function(dir, bandK, stamp, check = TRUE) {
  f <- file.path(dir, sprintf("stamp_band%02d.txt", bandK))
  if (check) {
    if (!file.exists(f) || !identical(readLines(f, warn = FALSE), stamp)) {
      unlink(Sys.glob(file.path(dir, sprintf("*_band%02d.tif", bandK))))
      message("Replicate set changed (or first run) -- band ", bandK, " summaries are rebuilt.")
    }
  } else writeLines(stamp, f)
}

#' Per-species, per-band summaries -> piece rasters (annual maps, change maps, per-pixel trend)
uncSummarizeSpeciesBand <- function(cfg, sp, band) {
  pieceDir <- uncSpDir(cfg, sp, "pieces")
  tmpl <- band$template
  stamp <- paste(uncRunIds(cfg, sp, band$k), collapse = ",")
  uncPieceStamp(pieceDir, band$k, stamp, check = TRUE)
  on.exit(uncPieceStamp(pieceDir, band$k, stamp, check = FALSE), add = TRUE)
  pf <- function(kind, name) file.path(pieceDir, sprintf("%s_%s_band%02d.tif", kind, name, band$k))
  yrs <- cfg$outYears

  # 1. per-year maps
  for (yr in yrs) {
    if (file.exists(pf("map", yr))) next
    p <- uncReadBandYear(cfg, sp, yr, band$k)
    if (is.null(p) || !length(p$idx) || ncol(p$P) < 2) next
    uncWritePiece(tmpl, p$idx, uncSummarizeMatrix(p$P, cfg$probs), pf("map", yr))
  }

  # 2. change maps (late minus early, per replicate)
  for (comp in names(uncComparisons(cfg))) {
    y0 <- uncComparisons(cfg)[[comp]]; y1 <- cfg$currentYear
    if (!(y0 %in% yrs && y1 %in% yrs) || file.exists(pf("change", comp))) next
    a <- uncReadBandYear(cfg, sp, y0, band$k); b <- uncReadBandYear(cfg, sp, y1, band$k)
    if (is.null(a) || is.null(b)) next
    common <- intersect(a$idx, b$idx)
    cols <- intersect(a$ids, b$ids)
    if (!length(common) || length(cols) < 2) next
    d <- b$P[match(common, b$idx), match(cols, b$ids), drop = FALSE] - a$P[match(common, a$idx), match(cols, a$ids), drop = FALSE]
    uncWritePiece(tmpl, common, uncSummarizeChange(d, cfg$probs), pf("change", comp))
  }

  # 3. per-pixel linear trend over all mapped years (probability per DECADE), per replicate
  if (length(yrs) >= 5 && !file.exists(pf("trend", "all"))) {
    tbar <- mean(yrs); w <- (yrs - tbar) / sum((yrs - tbar)^2) * 10
    cells <- NULL; S <- NULL; ids <- NULL
    for (i in seq_along(yrs)) {
      p <- uncReadBandYear(cfg, sp, yrs[i], band$k)
      if (is.null(p)) { cells <- integer(0); break }
      if (is.null(cells)) { cells <- p$idx; ids <- p$ids; S <- w[i] * p$P }
      else {
        common <- intersect(cells, p$idx); cols <- intersect(ids, p$ids)
        S <- S[match(common, cells), match(cols, ids), drop = FALSE] + w[i] * p$P[match(common, p$idx), match(cols, p$ids), drop = FALSE]
        cells <- common; ids <- cols
      }
    }
    if (length(cells) && length(ids) >= 2) uncWritePiece(tmpl, cells, uncSummarizeChange(S, cfg$probs, prefix = "slope"), pf("trend", "all"))
  }
  invisible(TRUE)
}

#' Community layers per band (needs all species): expected richness (sum of probabilities) and mean change in probability
uncCommunityBand <- function(cfg, band) {
  pieceDir <- file.path(uncRoot(cfg), "community", "pieces"); dir.create(pieceDir, recursive = TRUE, showWarnings = FALSE)
  tmpl <- band$template
  stamp <- paste(vapply(cfg$species, function(s) paste(uncRunIds(cfg, s, band$k), collapse = ","), character(1)), collapse = "|")
  uncPieceStamp(pieceDir, band$k, stamp, check = TRUE)
  on.exit(uncPieceStamp(pieceDir, band$k, stamp, check = FALSE), add = TRUE)
  pf <- function(kind, name) file.path(pieceDir, sprintf("%s_%s_band%02d.tif", kind, name, band$k))
  yrs <- cfg$outYears

  accumulate <- function(items) {              # items: list of list(idx, P, ids); sum & count over species, common replicate ids
    items <- items[!vapply(items, is.null, logical(1))]
    if (!length(items)) return(NULL)
    ids <- Reduce(intersect, lapply(items, `[[`, "ids"))
    if (length(ids) < 2) return(NULL)
    cells <- sort(Reduce(union, lapply(items, `[[`, "idx")))
    S <- matrix(0, length(cells), length(ids)); N <- integer(length(cells))
    for (it in items) {
      m <- match(it$idx, cells); P <- it$P[, match(ids, it$ids), drop = FALSE]
      S[m, ] <- S[m, ] + P; N[m] <- N[m] + 1L
    }
    list(cells = cells, S = S, N = N, ids = ids, nSpecies = length(items))
  }

  y1 <- cfg$currentYear
  if (y1 %in% yrs && !file.exists(pf("SRprob", y1))) {
    acc <- accumulate(lapply(cfg$species, function(sp) uncReadBandYear(cfg, sp, y1, band$k)))
    if (!is.null(acc)) uncWritePiece(tmpl, acc$cells, uncSummarizeMatrix(acc$S, cfg$probs), pf("SRprob", y1))
  }
  for (comp in names(uncComparisons(cfg))) {
    y0 <- uncComparisons(cfg)[[comp]]
    if (!(y0 %in% yrs && y1 %in% yrs) || file.exists(pf("meanDeltaP", comp))) next
    items <- lapply(cfg$species, function(sp) {
      a <- uncReadBandYear(cfg, sp, y0, band$k); b <- uncReadBandYear(cfg, sp, y1, band$k)
      if (is.null(a) || is.null(b)) return(NULL)
      common <- intersect(a$idx, b$idx); cols <- intersect(a$ids, b$ids)
      if (!length(common) || length(cols) < 2) return(NULL)
      list(idx = common, ids = cols,
           P = b$P[match(common, b$idx), match(cols, b$ids), drop = FALSE] - a$P[match(common, a$idx), match(cols, a$ids), drop = FALSE])
    })
    acc <- accumulate(items)
    if (!is.null(acc)) uncWritePiece(tmpl, acc$cells, uncSummarizeChange(acc$S / acc$N, cfg$probs), pf("meanDeltaP", comp))
  }
  invisible(TRUE)
}

#' Stitch the band pieces of one kind into one compressed GeoTIFF (virtual mosaic -> file)
uncStitch <- function(pieceFiles, out) {
  pieceFiles <- pieceFiles[file.exists(pieceFiles)]
  if (!length(pieceFiles)) return(invisible(NULL))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  v <- terra::vrt(pieceFiles, file.path(tempdir(), paste0("stitch_", basename(out), ".vrt")), overwrite = TRUE)
  names(v) <- names(terra::rast(pieceFiles[1]))      # the virtual mosaic loses the layer names
  terra::writeRaster(v, out, datatype = "FLT4S", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  invisible(out)
}

#' Area-mean series per replicate (all runs) from the per-band sums
#' @return data.frame: species, replicate, year, areaMean, nCells
uncAreaMeans <- function(cfg, sp) {
  rows <- list()
  nBands <- cfg$nBands
  for (lab in uncRunLabels(cfg, sp)) {
    dir <- file.path(uncSpDir(cfg, sp), "pred", lab)
    for (yr in cfg$outYears) {
      fs <- file.path(dir, sprintf("stats_%d_band%02d.rds", yr, seq_len(nBands)))
      fs <- fs[file.exists(fs)]
      if (length(fs) < nBands) { if (length(fs)) warning(sp, " ", lab, " ", yr, ": only ", length(fs), " of ", nBands, " bands -- year skipped"); next }
      st <- lapply(fs, readRDS)
      ids <- st[[1]]$ids
      tot <- Reduce(`+`, lapply(st, `[[`, "sum")); n <- sum(vapply(st, `[[`, numeric(1), "n"))
      rows[[length(rows) + 1]] <- data.frame(species = sp, replicate = ids, year = yr, areaMean = tot / n, nCells = n)
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Compare replicate 0 (the main models run through the uncertainty machinery) with the baseline meta-model maps
#'
#' A consistency check on real output: same valid cells and (almost) the same probabilities. Written to
#' `parity_check.txt` and as a warning if it fails. Skipped when replicate 0 or the baseline file is missing.
uncParityCheck <- function(cfg, sp, kMid = ceiling(cfg$nBands / 2)) {
  spClean <- gsub(" ", "_", sp)
  yr <- if (cfg$currentYear %in% cfg$outYears) cfg$currentYear else max(cfg$outYears)
  f <- file.path(cfg$outputRoot, paste0(metamodelLabel(cfg$sharedRes), if (!is.null(cfg$members)) paste0("_", ensTag(cfg$members)) else ""),
                 sprintf("%s_meta_suitability_%d.tif", spClean, yr))
  p <- uncReadBandYear(cfg, sp, yr, kMid, includeZero = TRUE)
  if (!file.exists(f) || is.null(p) || !(0 %in% p$ids)) {
    message(sp, ": parity check skipped (needs replicate 0, band ", kMid, " of year ", yr, " and the baseline file ", basename(f), ")")
    return(invisible(NULL))
  }
  band <- uncBands(cfg, sp)$bands[[kMid]]
  main <- terra::rast(f)[["meta_prob"]]
  mv <- terra::extract(main, terra::xyFromCell(band$template, p$idx))[, 1]
  r0 <- p$P[, which(p$ids == 0)]
  ok <- !is.na(mv) & !is.na(r0)
  lines <- c(sprintf("%s | year %d | band %d of %d | baseline file: %s", sp, yr, kMid, cfg$nBands, basename(f)),
             sprintf("  valid cells: replicate 0 = %d, baseline = %d, both = %d", sum(!is.na(r0)), sum(!is.na(mv)), sum(ok)),
             sprintf("  max |difference| = %.6f | mean |difference| = %.8f", max(abs(r0[ok] - mv[ok])), mean(abs(r0[ok] - mv[ok]))))
  relDiff <- NA_real_
  am <- file.path(uncSpDir(cfg, sp), "area_mean_replicates.csv")
  if (file.exists(am)) {
    a <- utils::read.csv(am); a0 <- a$areaMean[a$replicate == 0 & a$year == yr]
    if (length(a0)) {
      bm <- as.numeric(terra::global(main, "mean", na.rm = TRUE)[1, 1])
      relDiff <- abs(a0 - bm) / bm
      lines <- c(lines, sprintf("  country area mean: replicate 0 = %.6f, baseline = %.6f (relative difference %.2e)", a0, bm, relDiff))
    }
  }
  bad <- sum(xor(is.na(mv), is.na(r0))) > 0.001 * length(r0) || max(abs(r0[ok] - mv[ok])) > 1e-3 || (!is.na(relDiff) && relDiff > 1e-3)
  lines <- c(lines, if (bad) "  RESULT: MISMATCH -- do not trust the uncertainty layers until this is explained" else "  RESULT: OK")
  writeLines(lines, file.path(uncSpDir(cfg, sp), "parity_check.txt"))
  message(paste(lines, collapse = "\n"))
  if (bad) warning(sp, ": replicate 0 does not reproduce the baseline map (see parity_check.txt)", call. = FALSE)
  invisible(!bad)
}

#' Per-species assembly: stitch the band pieces into maps, write the per-replicate area means, run the parity check
#'
#' The area means (`area_mean_replicates.csv`) are what runIndex_Monitor's `computeIndexUncertainty` turns into index
#' intervals (species index, SBI, Chain); this module does not compute any index.
uncAssembleSpecies <- function(cfg, sp) {
  spClean <- gsub(" ", "_", sp)
  pieceDir <- uncSpDir(cfg, sp, "pieces"); outDir <- uncSpDir(cfg, sp, "maps")
  pcs <- function(kind, name) file.path(pieceDir, sprintf("%s_%s_band%02d.tif", kind, name, seq_len(cfg$nBands)))
  for (yr in cfg$outYears) uncStitch(pcs("map", yr), file.path(outDir, sprintf("%s_unc_%d.tif", spClean, yr)))
  for (comp in names(uncComparisons(cfg))) uncStitch(pcs("change", comp), file.path(outDir, sprintf("%s_unc_change_%s.tif", spClean, comp)))
  uncStitch(pcs("trend", "all"), file.path(outDir, sprintf("%s_unc_trend_per_decade.tif", spClean)))

  am <- uncAreaMeans(cfg, sp)
  if (!is.null(am)) {
    utils::write.csv(am, file.path(uncSpDir(cfg, sp), "area_mean_replicates.csv"), row.names = FALSE)
    message(sp, ": area means of ", length(unique(am$replicate)), " replicates x ", length(unique(am$year)), " years -> area_mean_replicates.csv")
  }
  tryCatch(uncParityCheck(cfg, sp), error = function(e) message(sp, ": parity check could not run: ", conditionMessage(e)))
  invisible(TRUE)
}
