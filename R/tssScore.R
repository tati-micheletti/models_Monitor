#' True Skill Statistic from a confusion matrix
#'
#' @param cmx A `PresenceAbsence::cmx()` confusion matrix.
#' @return Numeric, sensitivity + specificity - 1.
tssScore <- function(cmx) {
  PresenceAbsence::sensitivity(cmx, st.dev = FALSE) +
    PresenceAbsence::specificity(cmx, st.dev = FALSE) - 1
}
