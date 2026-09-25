#' Resolve the starting learning rate for a species' BRT optimization
#'
#' Priority: (1) an explicit per-species override in `perSpeciesLR` (intended
#' to eventually be fed by a per-species config, e.g. a future CSV -- see
#' improvements.md), (2) the learning rate a PRIOR run actually converged on
#' for this species (persisted by `persistConvergedLR()` after a fresh fit,
#' read back here so a forced refit -- new year's data appended, cache
#' invalidated -- doesn't blindly restart from `defaultLR` and re-walk the
#' same halving/doubling steps `optimizeBRT()` already found last time), (3)
#' `defaultLR` when neither of the above applies.
#'
#' @param sp Character. Species Latin name.
#' @param spClean Character. Filesystem-safe species name (`gsub(" ", "_", sp)`).
#' @param defaultLR Numeric. Fallback starting LR (the scale's existing
#'   `initialLR` default, e.g. 0.08 for habitat/landscape, 0.01 for Europe).
#' @param lrStateDir Character. Directory `persistConvergedLR()` writes to,
#'   typically the scale's own `outputDir`.
#' @param perSpeciesLR Named numeric vector/list, or NULL (default). Species
#'   not named here fall through to the persisted/default value.
#' @param lrStateSuffix Character. Distinguishes state files by scale (e.g.
#'   "habitat"/"landscape"/"EU") so the three scales' persisted LRs for the
#'   same species don't collide in the same directory.
#' @return Numeric. The learning rate to start `optimizeBRT()` from.
resolveStartingLR <- function(sp, spClean, defaultLR, lrStateDir, perSpeciesLR = NULL,
                               lrStateSuffix = "") {
  if (!is.null(perSpeciesLR) && sp %in% names(perSpeciesLR)) {
    return(as.numeric(perSpeciesLR[[sp]]))
  }

  lrStateFile <- file.path(lrStateDir, paste0(spClean, "_converged_lr", lrStateSuffix, ".rds"))
  if (file.exists(lrStateFile)) {
    persisted <- tryCatch(readRDS(lrStateFile), error = function(e) NULL)
    if (is.numeric(persisted) && length(persisted) == 1) {
      message("  Using previously converged LR for ", sp, ": ", persisted)
      return(persisted)
    }
  }

  defaultLR
}

#' Persist the learning rate a BRT fit converged on
#'
#' Called after a fresh (non-cached) `optimizeBRT()` fit, so a future forced
#' refit for this species can start from this value instead of `defaultLR`
#' via `resolveStartingLR()`.
#'
#' @param brtM A fitted `dismo::gbm.step()` model (`optimizeBRT()`'s return).
#' @param spClean Character. Filesystem-safe species name.
#' @param lrStateDir Character. Directory to write the state file in.
#' @param lrStateSuffix Character. See `resolveStartingLR()`.
persistConvergedLR <- function(brtM, spClean, lrStateDir, lrStateSuffix = "") {
  dir.create(lrStateDir, recursive = TRUE, showWarnings = FALSE)
  lrStateFile <- file.path(lrStateDir, paste0(spClean, "_converged_lr", lrStateSuffix, ".rds"))
  saveRDS(brtM$gbm.call$learning.rate, lrStateFile)
}
