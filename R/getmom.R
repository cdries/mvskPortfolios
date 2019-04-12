
getmom <- function(indmom, w, m1, M2, M3, M4) {

  ### input
  # indmom    : vector of length 4 with booleans selecting the
  #             order of the portfolio moments to compute
  # w         : weight vector
  # m1        : vector with expected returns
  # M2        : covariance matrix
  # M3        : coskewness matrix
  # M4        : cokurtosis matrix
  #
  ### output
  # vector with portfolio moments of order selected by indmom

  moms <- NULL
  if (indmom[1]) moms <- c(moms, sum(w * m1))
  if (indmom[2]) moms <- c(moms, PerformanceAnalytics:::portm2(w, M2))
  if (indmom[3]) moms <- c(moms, PerformanceAnalytics:::portm3(w, M3))
  if (indmom[4]) moms <- c(moms, PerformanceAnalytics:::portm4(w, M4))

  return (moms)
}
