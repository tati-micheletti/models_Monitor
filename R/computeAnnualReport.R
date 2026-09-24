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
#' **Combined-index methods** (all use the same per-species inputs, so
#' they're directly comparable -- report all, not just one, since they can
#' legitimately disagree at the margins and that disagreement is itself
#' informative):
#' - **SBI** (`computeCombinedIndex()`): Zbinden et al. (2005)/Swiss Bird
#'   Index -- geometric mean of baseline-normalized levels, bootstrap CI.
#' - **Analytical** (`computeCombinedIndexAnalytical()`): same geometric
#'   mean, but the CI is the OLD DDA/PECBMS analytical delta-method formula
#'   (Gregory et al. 2005, Appendix A).
#' - **MSI** (`computeCombinedIndexMSI()`): same geometric mean again, but
#'   the CI is the CURRENT DDA/PECBMS standard -- the MSI-tool's Monte
#'   Carlo method (Soldaat et al. 2017). Reported alongside Analytical, not
#'   instead of it, since both are legitimate DDA/PECBMS-lineage methods
#'   (old vs. current standard) and stakeholders may want to see both.
#' - **Chain** (`computeChainIndex()`): Living Planet Index chain method --
#'   geometric mean of year-to-year log-ratios, then cumulated. Baseline-
#'   year-independent by construction (see that function's docstring).
#'   Computed twice: once over `allYears` (matches the other methods, but
#'   rests on the same hindcast assumption they do for pre-`habitatYears`
#'   years) and once restricted to `restrictedYears` (if supplied) -- a
#'   robustness check using only non-hindcast years, at the cost of not
#'   reaching back to `baselineYear`.
#'
#' @param species Character vector of Latin species names.
#' @param baselineYear Integer. Index baseline (2005 -- the first year with
#'   genuinely real landscape-scale training data in this pipeline). Used
#'   by the SBI/Analytical/MSI methods; the Chain method ignores it.
#' @param currentYear Integer. The year this annual report is "for".
#' @param allYears Integer vector. Full year range to build the index time
#'   series over (typically `predictionYears`, e.g. 2005:currentYear).
#' @param restrictedYears Integer vector or NULL. If supplied (e.g.
#'   `habitatYears`, 2022:2025 -- the meta-model's real, non-hindcast
#'   training years), also computes the Chain index restricted to just
#'   these years, as a hindcast-robustness check. NULL skips this.
#' @param metaDir Character. `metaModel()`'s output directory.
#' @param outputDir Character. Directory to save this report's outputs in.
#' @param changeThresh Numeric. Passed to `computeChangeMaps()`.
#' @param nBoot Integer. Passed to `computeCombinedIndex()`/`computeChainIndex()`.
#' @param nSim Integer. Passed to `computeCombinedIndexMSI()`.
#' @return Invisibly, a list with `speciesIndex` (matrix), `combinedIndexSBI`,
#'   `combinedIndexAnalytical`, `combinedIndexMSI`, `combinedIndexChain`,
#'   `combinedIndexChainRestricted` (NULL if `restrictedYears` is NULL)
#'   (data.frames), `srMap` (SpatRaster), and `changeMaps` (named list by
#'   comparison: `vsBaseline`, `vs5YearsAgo`, `vsLastYear`, each as returned
#'   by `computeChangeMaps()`). Everything is also written to `outputDir`.
computeAnnualReport <- function(species, baselineYear, currentYear, allYears,
                                 metaDir, outputDir, restrictedYears = NULL,
                                 changeThresh = 0.05, nBoot = 999, nSim = 1000) {

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

  message("Computing combined multi-species index -- four methods...")

  message("  (1/4) SBI (Zbinden et al. 2005) -- geometric mean, bootstrap CI")
  combinedIndexSBI <- computeCombinedIndex(indexMat, nBoot = nBoot)
  combinedIndexSBI$SBI_smooth <- smoothIndex(combinedIndexSBI$year, combinedIndexSBI$SBI)
  utils::write.csv(combinedIndexSBI, file.path(outputDir, "combined_index_sbi.csv"), row.names = FALSE)
  message("  Saved -> combined_index_sbi.csv")

  message("  (2/4) Analytical (DDA/PECBMS OLD method, Gregory et al. 2005) -- delta-method CI")
  combinedIndexAnalytical <- computeCombinedIndexAnalytical(indexMat, seIndexMat)
  combinedIndexAnalytical$SBI_smooth <- smoothIndex(combinedIndexAnalytical$year, combinedIndexAnalytical$SBI)
  utils::write.csv(combinedIndexAnalytical, file.path(outputDir, "combined_index_analytical.csv"), row.names = FALSE)
  message("  Saved -> combined_index_analytical.csv")

  message("  (3/4) MSI (DDA/PECBMS CURRENT standard, Soldaat et al. 2017) -- Monte Carlo CI")
  combinedIndexMSI <- computeCombinedIndexMSI(indexMat, seIndexMat, nSim = nSim)
  combinedIndexMSI$SBI_smooth <- smoothIndex(combinedIndexMSI$year, combinedIndexMSI$SBI)
  utils::write.csv(combinedIndexMSI, file.path(outputDir, "combined_index_msi.csv"), row.names = FALSE)
  message("  Saved -> combined_index_msi.csv")

  message("  (4/4) Chain (Living Planet Index method) -- geometric mean of year-to-year log-ratios")
  combinedIndexChain <- computeChainIndex(indexMat, nBoot = nBoot)
  combinedIndexChain$chainIndex_smooth <- smoothIndex(combinedIndexChain$year, combinedIndexChain$chainIndex)
  utils::write.csv(combinedIndexChain, file.path(outputDir, "combined_index_chain.csv"), row.names = FALSE)
  message("  Saved -> combined_index_chain.csv")

  combinedIndexChainRestricted <- NULL
  if (!is.null(restrictedYears)) {
    restrictedKeys <- as.character(restrictedYears)
    availableKeys <- intersect(restrictedKeys, colnames(indexMat))
    if (length(availableKeys) >= 2) {
      message("  (bonus) Chain restricted to ", min(restrictedYears), "-", max(restrictedYears),
              " (non-hindcast years only) -- hindcast-robustness check")
      combinedIndexChainRestricted <- computeChainIndex(indexMat[, availableKeys, drop = FALSE], nBoot = nBoot)
      combinedIndexChainRestricted$chainIndex_smooth <- smoothIndex(
        combinedIndexChainRestricted$year, combinedIndexChainRestricted$chainIndex)
      utils::write.csv(combinedIndexChainRestricted,
                        file.path(outputDir, "combined_index_chain_restricted.csv"), row.names = FALSE)
      message("  Saved -> combined_index_chain_restricted.csv")
    } else {
      warning("restrictedYears has fewer than 2 years overlapping allYears -- skipping restricted chain index.")
    }
  }

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
                  combinedIndexMSI = combinedIndexMSI,
                  combinedIndexChain = combinedIndexChain,
                  combinedIndexChainRestricted = combinedIndexChainRestricted,
                  srMap = srMap, changeMaps = changeMaps))
}
