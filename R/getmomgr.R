
getmomgr <- function(indmom, w, m1, M2, M3, M4) {

  ### input
  # indmom    : vector of length 4 with booleans selecting the
  #             order of the gradient of the portfolio moments to compute
  # w         : weight vector
  # m1        : vector with expected returns
  # M2        : covariance matrix
  # M3        : coskewness matrix
  # M4        : cokurtosis matrix
  #
  ### output
  # jacobian or gradient of the portfolio moments with respect to w

  momsgr <- matrix(NA, nrow = 4, ncol = length(w))
  if (indmom[1]) momsgr[1,] <- m1
  if (indmom[2]) momsgr[2,] <- derportm2(w, M2)
  if (indmom[3]) momsgr[3,] <- derportm3(w, M3)
  if (indmom[4]) momsgr[4,] <- derportm4(w, M4)
  momsgr <- momsgr[indmom,]

  momsgr
}
