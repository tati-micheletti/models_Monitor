# Fast unit tests (seconds, no project data), from the repo root:
#   Rscript modules/models_Monitor/tests/uncertainty/test_units.R
suppressMessages({library(terra); library(gbm); library(glmnet)})
source("modules/models_Monitor/tests/uncertainty/helper.R")
ok <- TRUE
check <- function(label, cond) { cat(if (isTRUE(cond)) "ok   " else "FAIL ", label, "\n"); if (!isTRUE(cond)) ok <<- FALSE }

# --- seeds: reproducible, independent of the caller's RNG state, top-up safe ---------------------------------------
set.seed(1); before <- runif(1); set.seed(1)
a <- uncWithSeed(99, runif(3)); after <- runif(1)
check("uncWithSeed restores the caller's RNG state", identical(before, after))
check("uncWithSeed is reproducible", identical(a, uncWithSeed(99, runif(3))))
x <- runif(200, 0, 1e5); y <- runif(200, 0, 1e5); occ <- rbinom(200, 1, 0.5)
blk <- uncBlockIds(x, y, 2e4)
i1 <- uncRowIndex(blk, occ, 12345); i2 <- uncRowIndex(blk, occ, 12345); i3 <- uncRowIndex(blk, occ, 12346)
check("block draw reproducible for the same seed", identical(as.integer(i1), as.integer(i2)))
check("block draw differs for another seed", !identical(as.integer(i1), as.integer(i3)))
check("a block drawn twice contributes its rows twice (draw contains duplicates)", anyDuplicated(i1) > 0)
check("all rows of a drawn block enter together", all(vapply(split(i1, blk[i1]), function(v) length(unique(table(v))) == 1, logical(1))))
check("redraw when a draw has too few absences", {
  occ2 <- c(rep(0, 12), rep(1, 188)); r <- tryCatch(uncRowIndex(blk, occ2, 1, minPerClass = 10), error = function(e) NULL); is.null(r) || sum(occ2[r] == 0) >= 10 })
check("stableSeed keys give top-up-safe seeds (replicate 51 does not depend on 1:50)",
      stableSeed(c("bootFit", "A b", "habitat", 51)) == stableSeed(c("bootFit", "A b", "habitat", 51)))

# --- int16 storage round trip ----------------------------------------------------------------------------------
P <- matrix(runif(60), 12, 5); f <- tempfile(fileext = ".rds")
uncWriteInt16(f, idx = 101:112, P = P, ids = 1:5); r <- uncReadInt16(f)
check("int16 round trip within 1/30000", max(abs(r$P - P)) <= 0.5 / 30000 + 1e-12 && identical(r$idx, 101:112) && identical(r$ids, 1:5))
check("int16 file is small (2 bytes/value + gzip)", file.size(f) < 60 * 4 + 400)

# --- ridge by hand = predict.cv.glmnet -------------------------------------------------------------------------
set.seed(3); X <- cbind(runif(300), runif(300), runif(300)); yy <- rbinom(300, 1, plogis(-1 + 2 * X[, 1] - X[, 3]))
colnames(X) <- c("climate_mean_prob", "landscape_mean_prob", "habitat_mean_prob")
fit <- cv.glmnet(X, yy, family = "binomial", alpha = 0, nfolds = 10, standardize = TRUE)
co <- matrix(as.vector(coef(fit, s = fit$lambda.1se)), 1)
S <- list(clim = matrix(X[, 1]), land = matrix(X[, 2]), hab = matrix(X[, 3]))
check("uncApplyRidge equals predict(cv.glmnet, type = 'response')",
      max(abs(uncApplyRidge(S, co)[, 1] - as.vector(predict(fit, X, s = fit$lambda.1se, type = "response")))) < 1e-9)

# --- summaries ---------------------------------------------------------------------------------------------------
m <- matrix(runif(500 * 40), 500, 40); s <- uncSummarizeMatrix(m)
check("summary mean/lwr/upr/width", isTRUE(all.equal(unname(s[, "mean"]), rowMeans(m))) && isTRUE(all.equal(unname(s[, "lwr"]), unname(apply(m, 1, quantile, 0.05)))) &&
      isTRUE(all.equal(unname(s[, "width"]), unname(s[, "upr"] - s[, "lwr"]))))
d <- rbind(c(rep(-1, 20), rep(1, 20)), c(rep(1, 20), rep(-1, 20)), rep(c(-1, 1), 20))
sc <- uncSummarizeChange(d)
check("share decreasing / increasing", isTRUE(all.equal(unname(sc[, "shareDecrease"]), c(0.5, 0.5, 0.5))) && all(sc[, "shareDecrease"] + sc[, "shareIncrease"] == 1))
d2 <- matrix(-abs(rnorm(10 * 30)), 10, 30)
check("all-negative change -> shareDecrease 1", all(uncSummarizeChange(d2)[, "shareDecrease"] == 1))

# --- neighbours / windows -------------------------------------------------------------------------------------------
r <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20, crs = "EPSG:3035")
nb <- uncNeighbourCells(r, cellFromRowCol(r, 10, 10), reach = 2L)
check("+-2 cell neighbourhood has 25 cells", length(nb) == 25)
check("neighbourhood at the corner is clipped, not wrapped", length(uncNeighbourCells(r, 1, reach = 2L)) == 9)

# --- trend weights: slope per decade of an exact line --------------------------------------------------------------
yrs <- 2005:2025; w <- (yrs - mean(yrs)) / sum((yrs - mean(yrs))^2) * 10
check("trend weights recover slope per decade", abs(sum(w * (0.2 + 0.003 * yrs)) - 0.03) < 1e-12)

# --- stamp invalidation -----------------------------------------------------------------------------------------------
td <- tempfile(); dir.create(td); file.create(file.path(td, "map_2020_band01.tif"))
uncPieceStamp(td, 1, "1,2,3", check = TRUE); check("first run clears nothing valid (rebuilds)", !file.exists(file.path(td, "map_2020_band01.tif")))
file.create(file.path(td, "map_2020_band01.tif")); uncPieceStamp(td, 1, "1,2,3", check = FALSE); uncPieceStamp(td, 1, "1,2,3", check = TRUE)
check("same replicate set keeps the pieces", file.exists(file.path(td, "map_2020_band01.tif")))
uncPieceStamp(td, 1, "1,2,3,4", check = TRUE); check("new replicates invalidate the pieces", !file.exists(file.path(td, "map_2020_band01.tif")))

# --- stitching several band pieces keeps values, order and layer names ------------------------------------------------
td2 <- tempfile(); dir.create(td2); pcs <- character(0)
full <- rast(nrows = 30, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 30, crs = "EPSG:3035")
for (k in 1:3) {                                   # three bands of 10 rows, from the top
  tm <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 30 - 10 * k, ymax = 40 - 10 * k, crs = "EPSG:3035")
  m <- cbind(mean = rep(k, 90), lwr = rep(k * 10, 90)); idx <- 1:90   # last 10 cells of each band stay NA
  pcs[k] <- uncWritePiece(tm, idx, m, file.path(td2, sprintf("map_2020_band%02d.tif", k)))
}
st <- uncStitch(pcs, file.path(td2, "out.tif")); sr <- rast(st)
check("stitched mosaic has the full extent", all(as.vector(ext(sr)) == c(0, 10, 0, 30)) && nrow(sr) == 30)
check("stitched layer names are kept", identical(names(sr), c("mean", "lwr")))
check("stitched values sit in the right band", all(values(sr[[1]])[c(1, 101, 201)] == c(1, 2, 3)) && sum(is.na(values(sr[[1]]))) == 30)
check("no partial files left behind", length(list.files(td2, pattern = "part")) == 0)

# --- bootstrap fit with fixed hyperparameters ----------------------------------------------------------------------
set.seed(5); d <- data.frame(occurrence = rbinom(400, 1, 0.4), a = runif(400), b = runif(400)); d$occurrence <- rbinom(400, 1, plogis(3 * d$a - 1.5))
mainM <- list(gbm.call = list(best.trees = 50, learning.rate = 0.1, bag.fraction = 0.75, tree.complexity = 2))
mdl <- uncWithSeed(7, uncFitBRT(d, c("a", "b"), mainM, sample(400, 400, TRUE)))
check("replicate BRT uses the fixed tree count", mdl$n.trees == 50 && mdl$shrinkage == 0.1 && mdl$interaction.depth == 2)
pr <- uncPredictMany(list(mdl, mdl), d[, c("a", "b")], 50, 1L)
check("uncPredictMany returns cells x models", identical(dim(pr), c(400L, 2L)) && all(pr >= 0 & pr <= 1))

# --- index intervals (runIndex_Monitor::computeIndexUncertainty) on synthetic area means --------------------------------
source("modules/runIndex_Monitor/R/computeIndexUncertainty.R")
ud <- tempfile(); rd <- tempfile(); set.seed(11)
for (s in c("Aa bb", "Cc dd", "Ee ff")) {
  dir.create(file.path(ud, gsub(" ", "_", s)), recursive = TRUE)
  g <- expand.grid(replicate = 0:12, year = 2005:2009)
  g$species <- s; g$nCells <- 1e6
  g$areaMean <- 0.5 * (1 + 0.02 * (g$year - 2005)) * exp(rnorm(nrow(g), 0, 0.03))
  g$areaMean[g$replicate == 0] <- 0.5 * (1 + 0.02 * (g$year[g$replicate == 0] - 2005))
  utils::write.csv(g[, c("species", "replicate", "year", "areaMean", "nCells")], file.path(ud, gsub(" ", "_", s), "area_mean_replicates.csv"), row.names = FALSE)
}
res <- computeIndexUncertainty(c("Aa bb", "Cc dd", "Ee ff"), ud, rd, baselineYear = 2005, currentYear = 2009, nBands = 2)
st <- res$species; cb <- res$combined
check("species index: baseline year is exactly 100 with a zero-width interval", all(st$indexMean[st$year == 2005] == 100) && all(st$indexLwr[st$year == 2005] == 100))
check("species index interval brackets the replicate mean", all(st$indexLwr <= st$indexMean & st$indexMean <= st$indexUpr))
check("replicate 0 is excluded (12 replicates counted)", all(st$nReplicates == 12))
check("all four combined methods are reported", setequal(unique(cb$method), c("SBI", "Chain", "Analytical", "MSI")))
check("combined interval brackets the estimate (SBI, Chain, Analytical, MSI)", all(cb$lwr <= cb$estimate + 1e-9 & cb$estimate <= cb$upr + 1e-9))
check("index series rises with the simulated +2%/yr trend", all(diff(cb$estimate[cb$method == "SBI"]) > 0))
check("interval files are written", all(file.exists(file.path(rd, c("species_index_uncertainty.csv", "combined_index_uncertainty.csv")), file.path(ud, "combined_index_replicates.csv"))))
check("baseline year not mapped -> skipped without error", is.null(computeIndexUncertainty(c("Aa bb"), ud, rd, baselineYear = 1999, currentYear = 2009, nBands = 2)))

# --- Bray-Curtis turnover: a species swap (case A) vs no change (case B) ----------------------------------------------------
# 9 species at one pixel. 2005: species 1-5 present (p = 0.8). Case A 2025: species 1 stays, 2-5 lost, 6-9 gained. Case B 2025: unchanged.
p0 <- c(rep(0.8, 5), rep(0, 4)); pA <- c(0.8, rep(0, 4), rep(0.8, 4)); pB <- p0
bc <- function(a, b) uncBrayCurtis(matrix(sum(abs(b - a))), matrix(sum(a + b)))[1, 1]
check("swap of 4 of 5 species: same expected richness as no change ...", isTRUE(all.equal(sum(pA), sum(pB))))
check("... and mean change in probability cancels to 0 (cannot see the swap)", isTRUE(all.equal(mean(pA - p0), 0)))
check("Bray-Curtis: case A = 0.8, case B = 0", isTRUE(all.equal(bc(p0, pA), 0.8)) && bc(p0, pB) == 0)
check("Bray-Curtis = 1 - 2 sum(min) / (sum p0 + sum p1) (the textbook form)", isTRUE(all.equal(bc(p0, pA), 1 - 2 * sum(pmin(p0, pA)) / (sum(p0) + sum(pA)))))
check("Bray-Curtis = 1 when no species is shared, and stays defined when everything is 0",
      isTRUE(all.equal(bc(c(0.5, 0), c(0, 0.5)), 1)) && bc(c(0, 0), c(0, 0)) == 0)
check("uncBrayCurtis keeps the matrix shape (cells x replicates)", identical(dim(uncBrayCurtis(matrix(1:6 / 10, 3, 2), matrix(1, 3, 2))), c(3L, 2L)))
# the baseline community layer (runIndex_Monitor::computeChangeMaps) gives the same numbers on two pixels (cell 1 = case A, cell 2 = case B)
source("modules/runIndex_Monitor/R/combineLayersSafely.R"); source("modules/runIndex_Monitor/R/computeChangeMaps.R")
md <- tempfile(); dir.create(md); spp <- paste("Sp", letters[1:9])
for (i in 1:9) for (yr in c(2005, 2025)) {
  pr <- c(if (yr == 2005) p0[i] else pA[i], p0[i]); r <- terra::rast(nrows = 1, ncols = 2, nlyrs = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1)
  terra::values(r) <- cbind(pr, as.numeric(pr > 0.5)); names(r) <- c("meta_prob", "binary")
  terra::writeRaster(r, file.path(md, paste0(gsub(" ", "_", spp[i]), "_meta_suitability_", yr, ".tif")), overwrite = TRUE)
}
cm <- computeChangeMaps(spp, 2005, 2025, md)$community
check("baseline community change has 5 layers incl. turnoverBC", identical(names(cm), c("meanDeltaP", "netGainLoss", "gainCount", "lossCount", "turnoverBC")))
check("baseline: case A meanDeltaP = 0 but gainCount = lossCount = 4 and turnoverBC = 0.8; case B all 0",
      isTRUE(all.equal(as.numeric(terra::values(cm)[1, ]), c(0, 0, 4, 4, 0.8))) && isTRUE(all.equal(as.numeric(terra::values(cm)[2, ]), c(0, 0, 0, 0, 0))))

# --- regional index uncertainty: the matrix version equals the baseline raster function; summaries and shares -----------------------------
source("modules/runIndex_Monitor/R/aggregateSpeciesToGrid.R"); source("modules/runIndex_Monitor/R/computeGriddedCombinedIndex.R")
source("modules/runIndex_Monitor/R/computeRegionalIndexUncertainty.R")
set.seed(5); rd <- tempfile(); dir.create(rd); yrs <- 2020:2023; spp3 <- c("Aa bb", "Cc dd", "Ee ff")
for (sp in spp3) for (yr in yrs) {                                  # 4 x 4 pixels of 100 m -> 2 x 2 cells of 200 m
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 400, ymin = 0, ymax = 400, crs = "EPSG:3035", nlyrs = 2)
  pr <- runif(16, 0.05, 0.9); if (sp == "Ee ff" && yr == 2020) pr[1:4] <- 1e-9    # species Ee ff is absent from one cell in the baseline year
  terra::values(r) <- cbind(pr, 1); names(r) <- c("meta_prob", "binary")
  terra::writeRaster(r, file.path(rd, paste0(gsub(" ", "_", sp), "_meta_suitability_", yr, ".tif")))
}
want <- computeGriddedCombinedIndex(spp3, yrs, 2020, rd, 200)
arr <- lapply(spp3, function(sp) { g <- aggregateSpeciesToGrid(sp, yrs, rd, 200); a <- array(NA_real_, c(4, length(yrs), 2))
  for (i in seq_along(yrs)) a[, i, ] <- terra::values(g[[as.character(yrs[i])]])[, 1]; a })
got <- regionalCombineSpecies(arr, baseIdx = 1)
check("regional matrix index equals computeGriddedCombinedIndex (incl. a species dropped by the baseline floor)",
      isTRUE(all.equal(got[, , 1], unname(as.matrix(terra::values(want))), check.attributes = FALSE)) && identical(got[, , 1], got[, , 2]))
check("baseline year of the combined index is exactly 100 where species contribute", all(abs(got[, 1, 1] - 100) < 1e-9, na.rm = TRUE))
dd <- cbind(c(-1, -2, 3, 0), c(-1, 1, 3, 0), c(-1, 2, 3, 0))
ch <- regionalChange(dd, c(0.05, 0.95), "delta")
check("shares of replicates decreasing / increasing are exact", isTRUE(all.equal(unname(ch[, "shareDecrease"]), c(1, 1 / 3, 0, 0))) &&
      isTRUE(all.equal(unname(ch[, "shareIncrease"]), c(0, 2 / 3, 1, 0))))
check("change layers are named like the species change maps", identical(colnames(ch), c("deltaMean", "deltaSd", "deltaLwr", "deltaUpr", "deltaWidth", "shareDecrease", "shareIncrease")))
sm <- regionalSummary(matrix(1:10, 2, 5), c(0.05, 0.95))
check("summary: mean, sd, interval ordering and width", isTRUE(all.equal(unname(sm[1, "mean"]), 5)) && all(sm[, "lwr"] <= sm[, "upr"]) && isTRUE(all.equal(sm[, "width"], sm[, "upr"] - sm[, "lwr"])))

# --- German outline without GDAL/PROJ: the closed-form LAEA equals terra's projection (where terra's axis order is fine) ---------------------
check("uncLonLatToLAEA: (11.5E, 48.1N) -> easting 4432769 / northing 2777406 (reference values of tools/reprexTerraSfAxis.R)",
      all(abs(uncLonLatToLAEA(11.5, 48.1) - c(4432769, 2777406)) < 1))
if (requireNamespace("geodata", quietly = TRUE) && dir.exists("inputs/predictors/raw/gadm")) {
  ol <- uncOutlineLAEA("inputs/predictors/raw/gadm", "EPSG:3035")
  tp <- terra::project(geodata::gadm(country = "DEU", level = 0, path = "inputs/predictors/raw/gadm"), "EPSG:3035")
  check("closed-form German outline = terra::project outline (extent within 5 m; nothing swapped)",
        max(abs(as.vector(terra::ext(ol)) - as.vector(terra::ext(tp)))) < 5 && terra::ext(ol)$xmin > 3.9e6)
}

# --- nearest neighbour without sf ----------------------------------------------------------------------------------
set.seed(11); q <- cbind(runif(1300, 4e6, 4.6e6), runif(1300, 2.6e6, 3.6e6)); rf <- cbind(runif(900, 4e6, 4.6e6), runif(900, 2.6e6, 3.6e6))
nn <- uncNearestIdx(q, rf, chunk = 400L)
if (requireNamespace("sf", quietly = TRUE)) {
  sfNN <- as.integer(sf::st_nearest_feature(sf::st_as_sf(data.frame(x = q[, 1], y = q[, 2]), coords = c("x", "y"), crs = 3035),
                                            sf::st_as_sf(data.frame(x = rf[, 1], y = rf[, 2]), coords = c("x", "y"), crs = 3035)))
  check("uncNearestIdx equals sf::st_nearest_feature", identical(nn, sfNN))
}
check("uncNearestIdx brute force", identical(nn[c(1, 700, 1300)], vapply(c(1, 700, 1300), function(i) which.min((rf[, 1] - q[i, 1])^2 + (rf[, 2] - q[i, 2])^2), 1L)))

cat(if (ok) "\nALL UNIT TESTS PASSED\n" else "\nSOME UNIT TESTS FAILED\n")
quit(status = if (ok) 0 else 1)
