#' Explained deviance of a model's predictions against observed values
#'
#' Used both as a block-CV variable-importance score (see inputs_Monitor's
#' `computeUnivarCV()`) and as a standalone performance metric (see
#' models_Monitor's model-fitting functions). Follows Wiedenroth et al.
#' NOTE: does not clamp to non-negative itself -- callers that need that
#' (e.g. `computeUnivarCV()`) do it themselves.
#'
#' NOTE: this exact function is deliberately duplicated verbatim in
#' models_Monitor (same filename) so neither module has a load-time
#' dependency on the other for it. Keep both copies byte-identical --
#' tools/check_duplicated_functions.R checks this.
#'
#' @param obs Numeric vector of observed 0/1 responses.
#' @param pred Numeric vector of predicted probabilities.
#' @param family Character. GLM family, default "binomial".
#' @param weights Numeric vector of observation weights.
#' @return Numeric, explained deviance (D2).
explDeviance <- function(obs, pred, family = "binomial", weights = rep(1, length(obs))) {
  if (family == "binomial") {
    pred <- ifelse(pred < 0.00001, 0.00001, ifelse(pred > 0.9999, 0.9999, pred))
  }
  nullPred <- rep(mean(obs), length(obs))
  1 - (dismo::calc.deviance(obs, pred, family = family, weights = weights) /
         dismo::calc.deviance(obs, nullPred, family = family, weights = weights))
}
