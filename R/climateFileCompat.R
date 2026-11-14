#' Climate-scale file names: "climate" now, "EU" in the runs up to v3 (2026-11-12)
#'
#' Until v3 the files of the climate-scale BRT carried the tag `EU` (`<sp>_BRT_EU.rds`, `<sp>_perf_EU.rds`, `<sp>_expl_dev_EU.rds`, `<sp>_converged_lr_EU.rds`,
#' `<sp>_pred_EU_<year>.tif`) while every other climate file (fold models, the ensemble members, `oofhab`, `cvpred`) and every file of the other two scales
#' used the scale name. New runs write `climate` everywhere. These two functions let the code read the files of finished runs (v1 to v3), which are never renamed.
#' `legacyClimateName()` gives the old name of a new-style path (a path that is not one of those five kinds is returned unchanged);
#' `climateFileExisting()` gives the path that exists: the new name, else the old one, else the new name (so that "missing" messages name the new file).
legacyClimateName <- function(f) {
  g <- sub("_(BRT|perf|expl_dev|converged_lr)_climate[.]rds$", "_\\1_EU.rds", f)
  sub("_pred_climate_([0-9]+)[.]tif$", "_pred_EU_\\1.tif", g)
}

climateFileExisting <- function(f) {
  if (file.exists(f)) return(f)
  g <- legacyClimateName(f)
  if (!identical(g, f) && file.exists(g)) g else f
}
