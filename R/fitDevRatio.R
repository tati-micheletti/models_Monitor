#' Fit a ridge regression and return its explained deviance ratio
#'
#' @param X Numeric matrix of predictors.
#' @param y Numeric vector of 0/1 responses.
#' @param lambda Numeric. Ridge penalty.
#' @return Numeric, `glmnet`'s `dev.ratio`.
fitDevRatio <- function(X, y, lambda) {
  fit <- glmnet::glmnet(x = X, y = y, family = "binomial", alpha = 0,
                         lambda = lambda, standardize = TRUE)
  fit$dev.ratio
}
