# Local test, step 0: a MINI BASELINE for one species with the CURRENT training tables, built with the baseline's own
# functions (fit BRTs -> per-scale predictions -> metaModel()). Output folder outputs/utest. The uncertainty workflow's
# replicate 0 (= main models, no resampling) must reproduce its meta-model maps (step 2 checks that).
Sys.setenv(BIRDMONITOR_RUNNAME = "utest", BIRDMONITOR_SPECIES = "Alauda arvensis", BIRDMONITOR_UNC_REPS = "0:3",
           BIRDMONITOR_UNC_YEARS = "2020:2025", BIRDMONITOR_UNC_BANDS = "40")
suppressMessages({library(terra); library(gbm); library(dismo); library(glmnet); library(reproducible)})
source("modules/models_Monitor/tests/uncertainty/helper.R")
cfg <- testCfg(); sp <- cfg$species; spClean <- gsub(" ", "_", sp)
Y <- cfg$outYears; hy <- cfg$habitatYearsOf(sp)
res <- cfg$resOf(sp)
cat("years", Y, "| habitat years", hy, "\n")
uncCovcache(cfg)    # habitat covariate stacks for all needed years (cached ones are skipped)

lr <- c(climate = 0.01, landscape = 0.05, habitat = 0.05)
models <- list()
for (sc in c("climate", "landscape", "habitat")) {
  out <- file.path(uncMainDir(cfg, sp, sc)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(out, paste0(spClean, "_BRT_", .uncModelSuffix[[sc]], ".rds"))
  spPa <- uncTrainingTable(cfg, sp, sc)
  predSel <- readRDS(file.path("inputs", "model_ready", scaleLabel(res[[sc]]), paste0(spClean, "_predictors.rds")))
  if (file.exists(f)) { models[[sc]] <- readRDS(f); next }
  set.seed(stableSeed(c("gbm.step", sp, sc)))
  t0 <- Sys.time()
  m <- dismo::gbm.step(data = spPa[, c("occurrence", predSel)], gbm.x = predSel, gbm.y = "occurrence", family = "bernoulli",
                       tree.complexity = 2, learning.rate = lr[[sc]], bag.fraction = 0.75, n.folds = 5, silent = TRUE, plot.main = FALSE)
  if (is.null(m)) stop("gbm.step failed for ", sc)
  cat(sc, ": trees", m$gbm.call$best.trees, "| predictors", length(predSel), "|", round(difftime(Sys.time(), t0, units = "secs")), "s\n")
  saveRDS(m, f); models[[sc]] <- m
}

# per-scale predictions with the baseline's predictBRTToRaster()
procRoot <- file.path("inputs", "predictors", "processed")
for (yr in Y) {
  # climate
  cl <- file.path(uncMainDir(cfg, sp, "climate"), sprintf("%s_pred_EU_%d.tif", spClean, yr))
  if (!file.exists(cl)) {
    lab <- scaleLabel(res[["climate"]]); bf <- file.path(procRoot, lab, paste0("bioclim_", yr - 5, "-", yr, "_", lab, ".tif"))
    terra::writeRaster(predictBRTToRaster(terra::rast(bf), models$climate$var.names, models$climate, 0.5), cl, overwrite = TRUE)
  }
  # landscape
  ls <- file.path(uncMainDir(cfg, sp, "landscape"), sprintf("%s_pred_landscape_%d.tif", spClean, yr))
  if (!file.exists(ls)) {
    cv <- loadCovariates(yr, file.path(procRoot, scaleLabel(res[["landscape"]])), NULL); names(cv) <- gsub("_[0-9]{4}$", "", names(cv))
    terra::writeRaster(predictBRTToRaster(cv, models$landscape$var.names, models$landscape, 0.5), ls, overwrite = TRUE)
  }
  # habitat (from the covariate cache = the same layers as loadHabitatCovariates())
  hb <- file.path(uncMainDir(cfg, sp, "habitat"), sprintf("%s_pred_habitat_%d.tif", spClean, yr))
  if (!file.exists(hb)) {
    t0 <- Sys.time()
    terra::writeRaster(predictBRTToRaster(terra::rast(uncCovcacheFile(cfg, sp, yr)), models$habitat$var.names, models$habitat, 0.5), hb, overwrite = TRUE)
    cat("habitat prediction", yr, round(difftime(Sys.time(), t0, units = "secs")), "s\n")
  }
}

# baseline meta-model on the Germany window of the reference grid (same cells as the full reference grid)
habLab <- scaleLabel(res[["habitat"]])
ref <- terra::rast(file.path(procRoot, habLab, paste0("solar_radiation_habitat_", habLab, ".tif")))
win <- terra::crop(ref, terra::ext(terra::rast(uncCovcacheFile(cfg, sp, Y[1]))), snap = "out")
metaDir <- file.path("outputs", "utest", "metamodel_02_1_50")
habTable <- uncTrainingTable(cfg, sp, "habitat")
t0 <- Sys.time()
mm <- metaModel(inputsDataGerHabitat = setNames(list(list(data = habTable)), sp), habitatYears = hy, predictionYears = Y,
                modelDirs = list(europe = uncMainDir(cfg, sp, "climate"), landscape = uncMainDir(cfg, sp, "landscape"), habitat = uncMainDir(cfg, sp, "habitat")),
                refRaster = win, outputDir = metaDir, cachePath = file.path(tempdir(), "utest_cache"))
cat("metaModel secs:", round(difftime(Sys.time(), t0, units = "secs")), "\n")
cat("DONE step 0\n")
