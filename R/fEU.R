
fEU <- function(w, gamma, M2, M3, M4, m1 = NULL) {

  ### input
  # w         : weight vector
  # gamma     : coefficient of risk aversion
  # M2        : covariance matrix
  # M3        : coskewness matrix
  # M4        : cokurtosis matrix
  # m1        : optional vector with expected returns
  #
  ### output
  # objective : minus the expected utility (CRRA preferences)
  # gradient  : the gradient of the objective with respect to w

  if (is.null(m1)) m1 <- rep(0, length(w))

  # objective value
  mom1 <- sum(w * m1)
  mom2 <- portm2(w, M2)
  mom3 <- portm3(w, M3)
  mom4 <- portm4(w, M4)
  obj <- -mom1 + gamma * mom2 / 2 - gamma * (gamma + 1) * mom3 / 6 +
    gamma * (gamma + 1) * (gamma + 2) * mom4 / 24

  # gradient
  momsgrad1 <- m1
  momsgrad2 <- derportm2(w, M2)
  momsgrad3 <- derportm3(w, M3)
  momsgrad4 <- derportm4(w, M4)
  gr <- -m1 + gamma * momsgrad2 / 2 - gamma * (gamma + 1) * momsgrad3 / 6 +
    gamma * (gamma + 1) * (gamma + 2) * momsgrad4 / 24

  list(objective = obj, gradient = gr)
}
