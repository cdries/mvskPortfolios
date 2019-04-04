
fHI <- function(w) {

  ### input
  # w         : weight vector
  #
  ### output
  # objective : Herfindahl index
  # gradient  : the gradient of the objective with respect to w

  # objective value
  obj <- sum(w^2)

  # gradient
  gr <- 2 * w

  return (list("objective" = obj, "gradient" = gr))
}
