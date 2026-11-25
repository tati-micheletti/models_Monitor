#' Get the CORINE snapshot year covering a given occurrence year
#'
#' The latest snapshot that is not after the occurrence year (the first snapshot for years before it);
#' a year after the last snapshot, e.g. 2026, uses the last one (see `corineSnapshotYears()`).
#'
#' @param occYear Integer. Occurrence year.
#' @return Integer CORINE snapshot year.
corineYear <- function(occYear) {
  snaps <- corineSnapshotYears()
  snaps[max(findInterval(occYear, snaps), 1L)]
}
