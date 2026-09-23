#' Build the full annual multi-species reporting package from meta-model outputs
#'
#' Ties together every output discussed with the co-authors, all derived
#' from already-computed `metaModel()` rasters -- no new model fitting:
#' (1) per-species indices and the Zbinden-et-al.-style combined index
#' (never presented without the per-species series alongside it,
#' per their own caution), (2) a species richness map for the current
#' year, and (3) prevalence-change / gain-loss / stable-increase-decrease
#' maps for three comparisons: vs. baseline, vs. 5 years ago, vs. last
#' year -- each with both a per-species and a community layer.
#'
#' @param species Character vector of Latin species names.
#' @param baselineYear Integer. Index baseline (2005 -- the first year with
#'   genuinely real landscape-scale training data in this pipeline).
#' @param currentYear Integer. The year this annual report is "for".
#' @param allYears Integer vector. Full year range to build the index time
#'   series over (typically `predictionYears`, e.g. 2005:currentYear).
#' @param metaDir Character. `metaModel()`'s output directory
#'   (`outputs/<runName>/models_Monitor/meta`).
#' @param outputDir Character. Directory to save this report's outputs in
#'   (e.g. `outputs/<runName>/models_Monitor/indicators`).
#' @param changeThresh Numeric. Passed to `computeChangeMaps()`.
#' @param nBoot Integer. Passed to `computeCombinedIndex()`.
#' @return Invisibly, a list with `speciesIndex` (matrix), `combinedIndex`
#'   (data.frame), `srMap` (SpatRaster), and `changeMaps` (named list by
#'   comparison: `vsBaseline`, `vs5YearsAgo`, `vsLastYear`, each as
#'   returned by `computeChangeMaps()`). Everything is also written to
#'   `outputDir`.
computeAnnualReport <- function(species, baselineYear, currentYear, allYears,
                                 metaDir, outputDir, changeThresh = 0.05, nBoot = 999) {

  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  ## 1. Per-species + combined index -----------------------------------------
  message("Extracting per-species occurrence-probability series (", min(allYears),
          "-", max(allYears), ")...")

  rawList <- lapply(species, extractMetaProbSeries, years = allYears, metaDir = metaDir)
  names(rawList) <- species

  indexList <- lapply(rawList, computeSpeciesIndex, baselineYear = baselineYear)
  indexMat <- do.call(rbind, indexList)
  rownames(indexMat) <- species

  utils::write.csv(indexMat, file.path(outputDir, "species_index.csv"))
  message("Saved -> species_index.csv")

  message("Computing combined (Zbinden-style) multi-species index...")
  combinedIndex <- computeCombinedIndex(indexMat, nBoot = nBoot)
  utils::write.csv(combinedIndex, file.path(outputDir, "combined_index.csv"), row.names = FALSE)
  message("Saved -> combined_index.csv")

  ## 2. Species richness map for the current year ------------------------------
  message("Computing species richness map for ", currentYear, "...")
  srMap <- computeSRMap(species, currentYear, metaDir)
  if (!is.null(srMap)) {
    terra::writeRaster(srMap, file.path(outputDir, paste0("SR_map_", currentYear, ".tif")),
                        overwrite = TRUE)
    message("Saved -> SR_map_", currentYear, ".tif")
  }

  ## 3. Change maps: vs. baseline, vs. 5 years ago, vs. last year ---------------
  comparisons <- list(
    vsBaseline  = baselineYear,
    vs5YearsAgo = currentYear - 5,
    vsLastYear  = currentYear - 1
  )

  changeMaps <- list()
  for (compName in names(comparisons)) {
    yearRef <- comparisons[[compName]]
    message("Computing change maps: ", yearRef, " -> ", currentYear, " (", compName, ")...")

    cm <- computeChangeMaps(species, yearRef, currentYear, metaDir, changeThresh = changeThresh)
    changeMaps[[compName]] <- cm

    if (!is.null(cm$community)) {
      terra::writeRaster(cm$community,
                          file.path(outputDir, paste0("change_", compName, "_community.tif")),
                          overwrite = TRUE)
    }
    for (sp in names(cm$perSpecies)) {
      if (is.null(cm$perSpecies[[sp]])) next
      spClean <- gsub(" ", "_", sp)
      terra::writeRaster(cm$perSpecies[[sp]],
                          file.path(outputDir, paste0("change_", compName, "_", spClean, ".tif")),
                          overwrite = TRUE)
    }
    message("Saved -> change_", compName, "_* (community + per-species)")
  }

  invisible(list(speciesIndex = indexMat, combinedIndex = combinedIndex,
                  srMap = srMap, changeMaps = changeMaps))
}
