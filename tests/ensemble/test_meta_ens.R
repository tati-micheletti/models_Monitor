# Local test of the ENSEMBLE meta-model on the mini baseline (Alauda arvensis, outputs/utest):
# per-scale models glm + nn (cheap) -> ensemble "ens_brt-glm-nn" at the three scales -> metaModel(scaleSource = "ens_brt-glm-nn").
#   Rscript modules/models_Monitor/tests/ensemble/test_meta_ens.R
# Needs step0_mini_baseline.R first (the BRT baseline) and the REBUILT climate windows.
Sys.setenv(BIRDMONITOR_RUNNAME = "utest", BIRDMONITOR_SPECIES = "Alauda arvensis", BIRDMONITOR_UNC_REPS = "0:3", BIRDMONITOR_UNC_YEARS = "2020:2025", BIRDMONITOR_UNC_BANDS = "40")
suppressMessages({library(terra); library(gbm); library(dismo); library(glmnet); library(reproducible); library(ranger); library(mgcv); library(nnet)})
source("modules/models_Monitor/tests/uncertainty/helper.R")
cfg <- testCfg(); sp <- cfg$species; Y <- cfg$outYears; hy <- cfg$habitatYearsOf(sp); res <- cfg$resOf(sp)
ok <- function(label, cond) { cat(if (isTRUE(cond)) "ok   " else "FAIL ", label, "\n"); if (!isTRUE(cond)) quit(status = 1) }
members <- c("brt", "glm", "nn"); tag <- ensTag(members); years <- sort(unique(c(Y, hy)))
cat("ensemble:", tag, "| years", paste(years, collapse = ","), "\n")
for (scale in c("climate", "landscape", "habitat")) {
  for (a in c("brt", "glm", "nn")) algoScaleRun(cfg, sp, scale, a, threads = 2L, years = years)
  ensembleScaleRun(cfg, sp, scale, members, years = years)
}

# the meta-model of the ensemble, next to the BRT one
habLab <- scaleLabel(res[["habitat"]])
ref <- terra::rast(file.path("inputs", "predictors", "processed", habLab, paste0("solar_radiation_habitat_", habLab, ".tif")))
win <- terra::crop(ref, terra::ext(terra::rast(uncCovcacheFile(cfg, sp, Y[1]))), snap = "out")
metaDir <- file.path("outputs", "utest", paste0("metamodel_02_1_50_", tag))
habTable <- uncTrainingTable(cfg, sp, "habitat")
mm <- metaModel(inputsDataGerHabitat = setNames(list(list(data = habTable)), sp), habitatYears = hy, predictionYears = Y,
                modelDirs = list(europe = uncMainDir(cfg, sp, "climate"), landscape = uncMainDir(cfg, sp, "landscape"), habitat = uncMainDir(cfg, sp, "habitat")),
                refRaster = win, outputDir = metaDir, cachePath = file.path(tempdir(), "utest_cache_ens"), honestCfg = cfg, scaleSource = tag)
r <- mm[[sp]]; m <- readRDS(r$modelPath); cb <- as.vector(stats::coef(m$model, s = m$lambda))
cat("ensemble meta coefficients (intercept, climate, landscape, habitat):", round(cb, 3), "\n")
cat("honest AUC", round(r$perf$AUC, 3), "| in-sample AUC", round(r$perfInSample$AUC, 3), "|", r$perf$evalBasis, "\n")
ok("weights never negative", all(cb[-1] >= 0))
ok("trained and evaluated on out-of-fold inputs", identical(r$perf$evalBasis, "out-of-fold inputs") && identical(m$trainBasis, "out-of-fold inputs"))
ok("one map per prediction year", length(list.files(metaDir, pattern = "meta_suitability.*tif$")) == length(Y))
oof <- metaOutOfFoldEnsemble(cfg, sp, tag)
ok("weights = ridge on the ensemble's out-of-fold inputs", isTRUE(all.equal(cb, oof$coefficients["outOfFold", ], tolerance = 1e-6, check.attributes = FALSE)))
cat("TEST OK\n")
