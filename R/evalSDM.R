#' Evaluate SDM predictions: AUC, TSS, Kappa, sensitivity/specificity, D2
#'
#' @param observation Numeric vector of observed 0/1 responses.
#' @param predictions Numeric vector of predicted probabilities.
#' @param thresh Numeric, optional fixed threshold. If NULL, chosen via
#'   `PresenceAbsence::optimal.thresholds()`.
#' @param threshMethod Character. Threshold selection method, default
#'   "MaxSens+Spec".
#' @param weights Numeric vector of observation weights.
#' @return One-row data.frame with AUC, TSS, Kappa, Sens, Spec, PCC, D2, thresh.
evalSDM <- function(observation, predictions, thresh = NULL,
                     threshMethod = "MaxSens+Spec", weights = rep(1, length(observation))) {
  threshDat <- data.frame(ID = seq_along(observation), obs = observation, pred = predictions)

  if (is.null(thresh)) {
    threshMat <- PresenceAbsence::optimal.thresholds(DATA = threshDat, req.sens = 0.85,
                                                       req.spec = 0.85, FPC = 1, FNC = 1)
    thresh <- threshMat[threshMat$Method == threshMethod, 2]
  }

  cmxOpt <- PresenceAbsence::cmx(DATA = threshDat, threshold = thresh)
  data.frame(AUC = PresenceAbsence::auc(threshDat, st.dev = FALSE),
             TSS = tssScore(cmxOpt),
             Kappa = PresenceAbsence::Kappa(cmxOpt, st.dev = FALSE),
             Sens = PresenceAbsence::sensitivity(cmxOpt, st.dev = FALSE),
             Spec = PresenceAbsence::specificity(cmxOpt, st.dev = FALSE),
             PCC = PresenceAbsence::pcc(cmxOpt, st.dev = FALSE),
             D2 = explDeviance(observation, predictions, weights = weights),
             thresh = thresh)
}
