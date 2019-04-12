
ieqTEvol <- function(w, M2, wref, maxTEvol) {

  ### input
  # w         : weight vector
  # M2        : covariance matrix
  # wref      : weight vector of reference portfolio
  # maxTEvol  : maximum tracking error volatility
  #
  ### output
  # objective : tracking error volatility minus maxTEvol
  # gradient  : the gradient of the objective with respect to w

  # objective value
  wdiff <- w - wref
  M2wdiff <- M2 %*% wdiff
  TEvol <- sqrt(sum(wdiff * M2wdiff))
  obj <- TEvol - maxTEvol

  # gradient
  gr <- M2wdiff / TEvol

  return (list("objective" = obj, "gradient" = gr))
}
