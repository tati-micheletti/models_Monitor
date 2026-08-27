#' Get the CORINE snapshot year covering a given occurrence year
#'
#' @param occYear Integer. Occurrence year.
#' @return Integer CORINE snapshot year. Defaults to the most recent
#'   snapshot (2018) if `occYear` is beyond the mapped ranges.
corineYear <- function(occYear) {
  yearMap <- corineYearMap()
  for (corineYr in names(yearMap)) {
    if (occYear %in% yearMap[[corineYr]]) {
      return(as.integer(corineYr))
    }
  }
  return(2018L)
}
