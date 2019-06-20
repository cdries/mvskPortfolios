###############################################################################
# Functions copied from PerformanceAnalytics in order to pass CRAN checks
###############################################################################

# Wrapper function for casting the coskewness matrix into the vector of unique elements
M3.mat2vec <- function(M3) {
  # M3        : numeric matrix of dimension p x p^2

  if (is.null(dim(M3))) stop("M3 must be a matrix")

  .Call('M3mat2vec', as.numeric(M3), NROW(M3), PACKAGE="mvskPortfolios")
}

# Wrapper function for casting the cokurtosis matrix into the vector of unique elements
M4.mat2vec <- function(M4) {
  # M4        : numeric matrix of dimension p x p^3

  if (is.null(dim(M4))) stop("M4 must be a matrix")

  .Call('M4mat2vec', as.numeric(M4), NROW(M4), PACKAGE="mvskPortfolios")
}

# portfolio covariance and gradient
portm2 <- function(w, sigma) {
  return(as.numeric(t(w) %*% sigma %*%w))
}

derportm2 <- function(w, sigma) {
  return(2 * sigma %*% w)
}

# portfolio third order central moment and gradient
#' @useDynLib mvskPortfolios
portm3 <- function(w, M3) {
  w <- as.numeric(w)
  if (NCOL(M3) != 1) M3 <- M3.mat2vec(M3)
  .Call('M3port', w, M3, length(w), PACKAGE="mvskPortfolios")
}

derportm3 <- function(w, M3) {
  w <- as.numeric(w)
  if (NCOL(M3) != 1) M3 <- M3.mat2vec(M3)
  as.matrix(.Call('M3port_grad', w, M3, length(w), PACKAGE="mvskPortfolios"), ncol = 1)
}

# portfolio fourth order central moment and gradient
portm4 <- function(w, M4) {
  w <- as.numeric(w)
  if (NCOL(M4) != 1) M4 <- M4.mat2vec(M4)
  .Call('M4port', w, M4, length(w), PACKAGE="mvskPortfolios")
}

derportm4 <- function(w, M4) {
  w <- as.numeric(w)
  if (NCOL(M4) != 1) M4 <- M4.mat2vec(M4)
  as.matrix(.Call('M4port_grad', w, M4, length(w), PACKAGE="mvskPortfolios"), ncol = 1)
}
