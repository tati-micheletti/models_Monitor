#' Build the full annual multi-species reporting package from meta-model outputs
#'
#' Ties together every output discussed with the co-authors, all derived
#' from already-computed `metaModel()` rasters -- no new model fitting:
#' (1) per-species indices and THREE combined-index methods (never
#' presented without the per-species series alongside them, per Zbinden et
#' al.'s own caution) -- see the "Combined-index methods" note below,
#' (2) a species richness map for the current year, and (3) prevalence-
#' change / gain-loss / stable-increase-decrease maps for three
#' comparisons: vs. baseline, vs. 5 years ago, vs. last year -- each with
#' both a per-species and a community layer.
#'
#' **Combined-index methods** (all three use the same per-species inputs,
#' so they're directly comparable -- report all three, not just one, since
#' they can legitimately disagree at the margins and that disagreement is
#' itself informative):
#' - **SBI** (`computeCombinedIndex()`): Zbinden et al. (2005)/Swiss Bird
#'   Index -- geometric mean of baseline-normalized levels, bootstrap CI.
#' - **Analytical** (`computeCombinedIndexAnalytical()`): same geometric
#'   mean, but the CI is the DDA/PECBMS analytical delta-method formula
#'   (Gregory et al. 2005, Appendix A) instead of a bootstrap.
#' - **Chain** (`computeChainIndex()`): Living Planet Index chain method --
#'   geometric mean of year-to-year log-ratios, then cumulated. Baseline-
#'   year-independent by construction (see that function's docstring).
#'
#' @param species Character vector of Latin species names.
#' @param baselineYear Integer. Index baseline (2005 -- the first year with
#'   genuinely real landscape-scale training data in this pipeline). Used
#'   by the SBI and Analytical methods; the Chain method ignores it.
#' @param currentYear Integer. The year this annual report is "for".
#' @param allYears Integer vector. Full year range to build the index time
#'   series over (typically `predictionYears`, e.g. 2005:currentYear).
#' @param metaDir Character. `metaModel()`'s output directory.
#' @param outputDir Character. Directory to save this report's outputs in.
#' @param changeThresh Numeric. Passed to `computeChangeMaps()`.
#' @param nBoot Integer. Passed to `computeCombinedIndex()`/`computeChainIndex()`.
#' @return Invisibly, a list with `speciesIndex` (matrix), `combinedIndexSBI`,
#'   `combinedIndexAnalytical`, `combinedIndexChain` (data.frames), `srMap`
#'   (SpatRaster), and `changeMaps` (named list by comparison: `vsBaseline`,
#'   `vs5YearsAgo`, `vsLastYear`, each as returned by `computeChangeMaps()`).
#'   Everything is also written to `outputDir`.
computeAnnualReport <- function(species, baselineYear, currentYear, allYears,
                                 metaDir, outputDir, changeThresh = 0.05, nBoot = 999) {

  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  ## 1. Per-species + combined index -----------------------------------------
  message("Extracting per-species occurrence-probability series (", min(allYears),
          "-", max(allYears), ")...")

  rawList <- lapply(species, extractMetaProbSeries, years = allYears, metaDir = metaDir)
  names(rawList) <- species
  rawSDList <- lapply(species, extractMetaProbSD, years = allYears, metaDir = metaDir)
  names(rawSDList) <- species

  baseKey <- as.character(baselineYear)

  indexList <- lapply(rawList, computeSpeciesIndex, baselineYear = baselineYear)
  indexMat <- do.call(rbind, indexList)
  rownames(indexMat) <- species

  seIndexList <- lapply(species, function(sp) {
    computeIndexSE(rawSDList[[sp]], rawBaseline = rawList[[sp]][[baseKey]])
  })
  names(seIndexList) <- species
  seIndexMat <- do.call(rbind, seIndexList)
  rownames(seIndexMat) <- species

  utils::write.csv(indexMat, file.path(outputDir, "species_index.csv"))
  message("Saved -> species_index.csv")

  message("Computing combined multi-species index -- three methods...")

  message("  (1/3) SBI (Zbinden et al. 2005) -- geometric mean, bootstrap CI")
  combinedIndexSBI <- computeCombinedIndex(indexMat, nBoot = nBoot)
  combinedIndexSBI$SBI_smooth <- smoothIndex(combinedIndexSBI$year, combinedIndexSBI$SBI)
  utils::write.csv(combinedIndexSBI, file.path(outputDir, "combined_index_sbi.csv"), row.names = FALSE)
  message("  Saved -> combined_index_sbi.csv")

  message("  (2/3) Analytical (DDA/PECBMS-style) -- same geometric mean, delta-method CI")
  combinedIndexAnalytical <- computeCombinedIndexAnalytical(indexMat, seIndexMat)
  combinedIndexAnalytical$SBI_smooth <- smoothIndex(combinedIndexAnalytical$year, combinedIndexAnalytical$SBI)
  utils::write.csv(combinedIndexAnalytical, file.path(outputDir, "combined_index_analytical.csv"), row.names = FALSE)
  message("  Saved -> combined_index_analytical.csv")

  message("  (3/3) Chain (Living Planet Index method) -- geometric mean of year-to-year log-ratios")
  combinedIndexChain <- computeChainIndex(indexMat, nBoot = nBoot)
  combinedIndexChain$chainIndex_smooth <- smoothIndex(combinedIndexChain$year, combinedIndexChain$chainIndex)
  utils::write.csv(combinedIndexChain, file.path(outputDir, "combined_index_chain.csv"), row.names = FALSE)
  message("  Saved -> combined_index_chain.csv")

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

  invisible(list(speciesIndex = indexMat, combinedIndexSBI = combinedIndexSBI,
                  combinedIndexAnalytical = combinedIndexAnalytical,
                  combinedIndexChain = combinedIndexChain,
                  srMap = srMap, changeMaps = changeMaps))
}
