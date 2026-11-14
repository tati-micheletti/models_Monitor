# Test of R/climateFileCompat.R: climate-scale BRT file names are "climate" now, "EU" in the runs up to v3. Run from the repo root:
#   Rscript modules/models_Monitor/tests/test_climateFileCompat.R
source("modules/models_Monitor/R/climateFileCompat.R")
chk <- function(ok, what) { if (!isTRUE(ok)) stop("FAILED: ", what, call. = FALSE); cat("ok:", what, "\n") }

## 1. the old name of the five kinds of climate BRT files
chk(legacyClimateName("a/Alauda_arvensis_BRT_climate.rds") == "a/Alauda_arvensis_BRT_EU.rds", "BRT model")
chk(legacyClimateName("a/Alauda_arvensis_perf_climate.rds") == "a/Alauda_arvensis_perf_EU.rds", "BRT performance")
chk(legacyClimateName("a/Alauda_arvensis_expl_dev_climate.rds") == "a/Alauda_arvensis_expl_dev_EU.rds", "explained deviance")
chk(legacyClimateName("a/Alauda_arvensis_converged_lr_climate.rds") == "a/Alauda_arvensis_converged_lr_EU.rds", "converged learning rate")
chk(legacyClimateName("a/Alauda_arvensis_pred_climate_2025.tif") == "a/Alauda_arvensis_pred_EU_2025.tif", "prediction map")

## 2. files that were always called "climate" (and the other scales) are left alone
for (f in c("Alauda_arvensis_foldModels_climate.rds", "Alauda_arvensis_glm_climate.rds", "Alauda_arvensis_perf_glm_climate.rds", "Alauda_arvensis_perf_ens_brt-glm-gam-rf-nn_climate.rds",
            "Alauda_arvensis_oofhab_rf_climate.rds", "Alauda_arvensis_cvpred_nn_climate.rds", "Alauda_arvensis_pred_glm_climate_2025.tif", "Alauda_arvensis_pred_ens_brt-glm-gam-rf-nn_climate_2025.tif",
            "Alauda_arvensis_BRT_landscape.rds", "Alauda_arvensis_BRT_habitat.rds", "Alauda_arvensis_pred_landscape_2025.tif"))
  chk(identical(legacyClimateName(f), f), paste("unchanged:", f))

## 3. the file that exists is chosen: new name, else old name, else the new name
d <- tempfile("climcompat"); dir.create(d)
fNew <- file.path(d, "Sp_one_BRT_climate.rds"); fOld <- file.path(d, "Sp_one_BRT_EU.rds")
chk(identical(climateFileExisting(fNew), fNew), "neither exists: the new name (so that messages name it)")
file.create(fOld); chk(identical(climateFileExisting(fNew), fOld), "only the old file exists: the old file")
file.create(fNew); chk(identical(climateFileExisting(fNew), fNew), "both exist: the new file")
p <- file.path(d, "Sp_one_pred_climate_2020.tif"); file.create(file.path(d, "Sp_one_pred_EU_2020.tif"))
chk(identical(climateFileExisting(p), file.path(d, "Sp_one_pred_EU_2020.tif")), "old prediction map found from the new name")
chk(identical(climateFileExisting(file.path(d, "Sp_one_pred_climate_2021.tif")), file.path(d, "Sp_one_pred_climate_2021.tif")), "a year without any file: the new name")
unlink(d, recursive = TRUE)
cat("all passed\n")
