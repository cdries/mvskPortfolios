
get_href <- function(href, m1, M2, M3, M4, param) {

  ### input
  # href      : name of the reference function ("DR", "ERC", "HI"/"EW", "TEvol)
  # m1        : vector with expected returns
  # M2        : covariance matrix
  # M3        : coskewness matrix
  # M4        : cokurtosis matrix
  # param     : list with additional parameters
  #
  ### output
  # objective function to minimize that retruns "objective" and "gradient"

  if (href == "DR") {
    fn <- function(w) fDR(w, M2)
  } else if (href == "ERC") {
    fn <- function(w) fERC(w, M2)
  } else if (href %in% c("HI", "EW")) {
    fn <- function(w) fHI(w)
  } else if (href == "EU") {
    fn <- function(w) fEU(w, param$gamma, M2, M3, M4, m1)
  } else if (href == "TEvol") {
    fn <- function(w) fTEvol(w, M2, param$wref)
  } else {
    stop("Choose a valid reference function")
  }

  fn
}
