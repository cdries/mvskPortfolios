
fERC <- function(w, M2) {

  # TODO

  ### input
  # w         : weight vector
  # M2        : covariance matrix
  #
  ### output
  # objective : minus the diversification ratio
  # gradient  : the gradient of the objective with respect to w

  # objective value
  s <- sqrt(diag(M2))
  wSigw <- sum(w * M2 %*% w)
  divratio <- sum(w * s) / sqrt(wSigw)
  obj <- -divratio

  # gradient
  gr <- -(s / sqrt(wSigw) - sum(w * s) / (wSigw)^1.5 * M2 %*% w)

  return (list("objective" = obj, "gradient" = gr))
}
