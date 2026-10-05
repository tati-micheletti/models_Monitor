#' Deterministic seed derived from a key (e.g. species + scale + purpose)
#'
#' Random steps (CV block-to-fold assignment, BRT bagging, CV folds of the ridge) are seeded from a
#' KEY such as `c("gbm.step", species, n, ...)`, not from the order in which species happen to run, so
#' the same data always give the same result, on a laptop and on a cluster task alike. Not meant to be
#' cryptographic: a simple, stable hash.
#'
#' @param key Character (or coercible) vector; pasted together.
#' @return Integer seed in [1, 2^31 - 2].
stableSeed <- function(key) {
  chars <- utf8ToInt(paste(key, collapse = "|"))
  as.integer(sum(chars * ((seq_along(chars) %% 97) + 1)) %% 2147483629 + 1)
}
