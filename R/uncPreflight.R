#' Pre-flight check: everything the uncertainty workflow reads must exist BEFORE jobs are submitted
#'
#' Cheap (file existence only, no models read), so it can run on a login node. Prints every problem found.
#'
#' @return Invisibly TRUE if all is fine; stops with a summary otherwise.
uncPreflight <- function(cfg) {
  problems <- character(0)
  add <- function(...) problems <<- c(problems, paste0(...))

  for (p in c("matrixStats", "gbm", "glmnet", "terra", "sf", "blockCV", "parallel"))
    if (!requireNamespace(p, quietly = TRUE)) add("R package not installed: ", p)

  for (sp in cfg$species) {
    res <- cfg$resOf(sp)
    for (sc in c("climate", "landscape", "habitat")) {
      f1 <- file.path(cfg$inputRoot, "model_ready", scaleLabel(res[[sc]]), paste0(gsub(" ", "_", sp), "_inputs.rds"))
      if (!file.exists(f1)) add(sp, ": training table missing: ", f1)
      f2 <- file.path(uncMainDir(cfg, sp, sc), paste0(gsub(" ", "_", sp), "_BRT_", .uncModelSuffix[[sc]], ".rds"))
      if (!file.exists(f2)) add(sp, ": main BRT missing (baseline model array not finished?): ", f2)
    }
    proc <- file.path(cfg$inputRoot, "predictors", "processed")
    for (yr in uncAllYears(cfg, sp)) {
      hab <- file.path(proc, scaleLabel(res[["habitat"]]), paste0("landuse_", yr, "_habitat_", scaleLabel(res[["habitat"]]), ".tif"))
      lan <- file.path(proc, scaleLabel(res[["landscape"]]), paste0("landuse_", yr, "_landscape_", scaleLabel(res[["landscape"]]), ".tif"))
      cli <- file.path(proc, scaleLabel(res[["climate"]]), paste0("bioclim_", yr - (cfg$climateWindowLength - 1), "-", yr, "_", scaleLabel(res[["climate"]]), ".tif"))
      need <- c(lan, cli)
      if (!file.exists(uncCovcacheFile(cfg, sp, yr))) need <- c(hab, need)   # the habitat stack is only rebuilt from these if not cached
      for (f in need) if (!file.exists(f)) add(sp, " ", yr, ": covariate file missing: ", f)
    }
    ref <- file.path(proc, scaleLabel(res[["habitat"]]), paste0("solar_radiation_habitat_", scaleLabel(res[["habitat"]]), ".tif"))
    if (!file.exists(ref)) add(sp, ": reference grid missing: ", ref)
  }
  habRes <- unique(vapply(cfg$species, function(s) cfg$resOf(s)[["habitat"]], numeric(1)))
  if (length(habRes) > 1) add("Species have different habitat resolutions (", paste(habRes, collapse = ", "), "): the band layout assumes one.")

  if (length(problems)) {
    problems <- unique(problems)
    message("PRE-FLIGHT FOUND ", length(problems), " PROBLEM(S):\n  - ", paste(utils::head(problems, 40), collapse = "\n  - "),
            if (length(problems) > 40) paste0("\n  ... and ", length(problems) - 40, " more") else "")
    stop("Pre-flight failed -- fix the above, then submit again.", call. = FALSE)
  }
  message("Pre-flight OK: ", length(cfg$species), " species, replicates ", cfg$repLabel, ", years ", min(cfg$outYears), "-", max(cfg$outYears),
          " (", length(cfg$outYears), " mapped), ", cfg$nBands, " bands.")
  invisible(TRUE)
}
