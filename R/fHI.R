
fHI <- function(w) {

  ### input
  # w         : weight vector
  #
  ### output
  # objective : normalized Herfindahl index
  # gradient  : the gradient of the objective with respect to w

  # objective value
  p <- length(w)
  obj <- (sum(w^2) - 1 / p) / (1 - 1 / p)

  # gradient
  gr <- 2 * w / (1 - 1 / p)

  return (list("objective" = obj, "gradient" = gr))
}
