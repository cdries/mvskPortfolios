
fERC <- function(w, M2) {

  ### input
  # w         : weight vector
  # M2        : covariance matrix
  #
  ### output
  # objective : equal-risk (standard deviation) contribution objective function
  # gradient  : the gradient of the objective with respect to w

  # objective value
  p <- length(w)
  Sigw <- M2 %*% w
  pctRC <- w * Sigw / sum(w * Sigw)
  obj <- sum((pctRC - 1 / p)^2)

  # gradient
  sp <- sqrt(sum(w * Sigw))

  dff <- pctRC - 1 / p
  gr <- 2 * (sp^2 * (M2 %*% (w * dff) + dff * Sigw) -
               2 * sum(w * dff * Sigw) * Sigw) / sp^4

  return (list("objective" = obj, "gradient" = gr))
}
