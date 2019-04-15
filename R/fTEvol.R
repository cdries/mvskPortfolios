
fTEvol <- function(w, M2, wref) {

  ### input
  # w         : weight vector
  # M2        : covariance matrix
  # wref      : weight vector of reference portfolio
  #
  ### output
  # objective : tracking error volatility minus maxTEvol
  # gradient  : the gradient of the objective with respect to w

  # objective value
  wdiff <- w - wref
  M2wdiff <- M2 %*% wdiff
  TEvol <- sqrt(sum(wdiff * M2wdiff))
  obj <- TEvol

  # gradient
  if (obj < 1e-10) gr <- rep(0, length(w)) else gr <- M2wdiff / TEvol

  list(objective = obj, gradient = gr)
}
