#' Temporal assignment of CORINE land cover snapshots to occurrence years
#'
#' Method: nearest preceding snapshot -- use the most recent CORINE map
#' available at the time of occurrence; never use a "future" land cover
#' map for historical data.
#' 2006 map: occurrence years 2005-2011
#' 2012 map: occurrence years 2012-2017
#' 2018 map: occurrence years 2018-2025
corineYearMap <- function() {
  list("2006" = 2005:2011,
       "2012" = 2012:2017,
       "2018" = 2018:2025)
}
