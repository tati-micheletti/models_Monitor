#' Temporal assignment of CORINE land cover snapshots to occurrence years
#'
#' Method: nearest preceding snapshot -- use the most recent CORINE map
#' available at the time of occurrence; never use a "future" land cover
#' map for historical data. Years before the first snapshot use the first one.
#' 2006 map: occurrence years 2005-2011
#' 2012 map: occurrence years 2012-2017
#' 2018 map: occurrence years 2018 and every later year (no end year written into the code, 2026-11-25)
#' Only the snapshot years are listed here; a new CORINE release is added to `corineSnapshotYears()`.
corineSnapshotYears <- function() c(2006L, 2012L, 2018L)

#' Snapshot -> occurrence years, derived from `corineSnapshotYears()` (names = snapshot years, as before)
#' @param years Integer vector of occurrence years to assign (default 2005:2100).
corineYearMap <- function(years = 2005:2100) {
  snaps <- corineSnapshotYears()
  assigned <- vapply(years, function(y) snaps[max(findInterval(y, snaps), 1L)], integer(1))
  split(years, factor(assigned, levels = snaps))
}
