# Local test: metaModel() with (a) non-negative ridge weights and (b) weights + accuracy on OUT-OF-FOLD inputs ("honest").
# Reads the mini baseline of step0 (outputs/utest, never modified) and writes to outputs/utest_honest.
Sys.setenv(BIRDMONITOR_RUNNAME = "utest", BIRDMONITOR_SPECIES = "Alauda arvensis", BIRDMONITOR_UNC_REPS = "0:3",
           BIRDMONITOR_UNC_YEARS = "2020:2025", BIRDMONITOR_UNC_BANDS = "40")
suppressMessages({library(terra); library(gbm); library(dismo); library(glmnet); library(reproducible)})
source("modules/models_Monitor/tests/uncertainty/helper.R")
cfg <- testCfg(); sp <- cfg$species; spClean <- gsub(" ", "_", sp)
Y <- cfg$outYears; hy <- cfg$habitatYearsOf(sp); res <- cfg$resOf(sp)
procRoot <- file.path("inputs", "predictors", "processed")
habLab <- scaleLabel(res[["habitat"]])
ref <- terra::rast(file.path(procRoot, habLab, paste0("solar_radiation_habitat_", habLab, ".tif")))
win <- terra::crop(ref, terra::ext(terra::rast(uncCovcacheFile(cfg, sp, Y[1]))), snap = "out")
outDir <- file.path("outputs", "utest_honest", "metamodel_02_1_50"); unlink(dirname(outDir), recursive = TRUE)
habTable <- uncTrainingTable(cfg, sp, "habitat")
mm <- metaModel(inputsDataGerHabitat = setNames(list(list(data = habTable)), sp), habitatYears = hy, predictionYears = Y,
                modelDirs = list(europe = uncMainDir(cfg, sp, "climate"), landscape = uncMainDir(cfg, sp, "landscape"), habitat = uncMainDir(cfg, sp, "habitat")),
                refRaster = win, outputDir = outDir, cachePath = file.path(tempdir(), "utest_cache_honest"), honestCfg = cfg)
r <- mm[[sp]]
cat("\n--- RESULT ---\n")
m <- readRDS(r$modelPath)
cb <- as.vector(stats::coef(m$model, s = m$lambda))
cat("coefficients (intercept, climate, landscape, habitat):", round(cb, 3), "\n")
stopifnot(all(cb[-1] >= 0))
cat("honest   :", r$perf$evalBasis, "| AUC", round(r$perf$AUC, 3), "TSS", round(r$perf$TSS, 3), "D2", round(r$perf$D2, 3), "\n")
cat("in-sample:", "AUC", round(r$perfInSample$AUC, 3), "TSS", round(r$perfInSample$TSS, 3), "D2", round(r$perfInSample$D2, 3), "\n")
stopifnot(identical(r$perf$evalBasis, "out-of-fold inputs"), identical(m$trainBasis, "out-of-fold inputs"))
chk <- readRDS(file.path(outDir, paste0(spClean, "_meta_check.rds")))
cat("weights = ridge on the check's out-of-fold inputs:", isTRUE(all.equal(cb, chk$coefficients["outOfFold", ], tolerance = 1e-6, check.attributes = FALSE)), "
")
stopifnot(isTRUE(all.equal(cb, chk$coefficients["outOfFold", ], tolerance = 1e-6, check.attributes = FALSE)))
# second call: nothing is rebuilt (same weights -> maps are reused)
mm2 <- metaModel(inputsDataGerHabitat = setNames(list(list(data = habTable)), sp), habitatYears = hy, predictionYears = Y,
                 modelDirs = list(europe = uncMainDir(cfg, sp, "climate"), landscape = uncMainDir(cfg, sp, "landscape"), habitat = uncMainDir(cfg, sp, "habitat")),
                 refRaster = win, outputDir = outDir, cachePath = file.path(tempdir(), "utest_cache_honest"), honestCfg = cfg)
cat("maps written:", length(list.files(outDir, pattern = "meta_suitability.*tif$")), "\n")
cat("TEST OK\n")
