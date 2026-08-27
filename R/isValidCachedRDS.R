#' Check a cached RDS object exists and is non-empty
#'
#' Unlike dataPrep_Monitor's `isValidRDSFile()` (which checks `nrow(x) > 0`,
#' appropriate for cached occurrence tables), this checks `length(x) > 0` --
#' appropriate here because models_Monitor caches non-tabular objects too
#' (a fitted `gbm.step()` model is a list with no meaningful `nrow`; calling
#' `nrow()` on one returns NULL, and `NULL > 0` errors). Deliberately named
#' differently from dataPrep_Monitor's version so the two are never mistaken
#' for duplicates of the same function -- see tools/check_duplicated_functions.R.
#'
#' @param path Character. Path to an RDS file.
#' @return Logical.
isValidCachedRDS <- function(path) {
  if (!file.exists(path)) return(FALSE)
  tryCatch({
    x <- readRDS(path)
    length(x) > 0
  }, error = function(e) FALSE)
}
