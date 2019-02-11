#' Risk-based higher-order moment efficient portfolios
#'
#' function determining risk based portfolios that are mean-variance-skewness-kurtosis efficient.
#'
#' moments with NULL are ignored; default is maxDiv as starting point; unconstrained, just set kappa to 1
#'
#' possible initial portfolios are "maxDiv" for the maximum diversified portfolio of REF.
#'
#' the vector g should have positive values, the default value is g = (0, sigma_^2_w0, abs(phi_w0), psi_w0)
#'
#' lin_eq and lin_eqC default to the full investment constraint
#'
#' @name mvskPortfolio
#' @encoding UTF-8
#' @concept risk-based portfolio
#' @param m1 mean vector; not used if NULL
#' @param M2 covariance matrix
#' @param M3 coskewness matrix; not used if NULL
#' @param M4 cokurtosis matrix; not used if NULL
#' @param w0 weight of the benchmark portfolio / initial portfolio name, see details
#' @param g vector with preferences for the moments, see details
#' @param kappa vector of values indicating maximum deviation of the benchmark portfolio
#' @param lb lower bound for the weights, default rep(0, p)
#' @param ub upper bound for the weights, default rep(1, p)
#' @param lin_eq equality constraints: eq w = eqC (should be matrix!), see details
#' @param lin_eqC equality constraints: eq w = eqC, see details
#' @param nlin_eq function with non-linear equality constraint (returns objective value and jacobian)
#' @param lin_ieq inequality constraints: ieq w leq ieqC (should be matrix!)
#' @param lin_ieqC inequality constraints: ieq w leq ieqC
#' @param nlin_ieq function with non-linear inequality constraints (returns objective value and jacobian)
#' @param riskcriterion optimal value of kappa minimizes the risk criterion function
#' @param options optimization options
#' @param relative steps of kappa are procentual if TRUE, absolute if FALSE
#' @author Dries Cornilly
#' @references
#' Choueifaty, Y., & Coignard, Y. (2008). Toward maximum diversification.
#' Journal of Portfolio Management, 35(1), 40.
#'
#' Briec, W., Kerstens, K., & Jokung, O. (2007). Mean-variance-skewness portfolio
#' performance gauging: a general shortage function and dual approach.
#' Management science, 53(1), 135-149.
#'
#' @examples
#' # load data
#' data(indices)
#' p <- 5
#' x <- as.matrix(indices[, 1:p])
#'
#' # estimate moments
#' M2 <- cov(x)
#' M3 <- PerformanceAnalytics:::M3.MM(x, as.mat = FALSE)
#' M4 <- PerformanceAnalytics:::M4.MM(x, as.mat = FALSE)
#'
#' # benchmark portfolio
#' w0 <- "maxDiv"
#'
#' # optimal VSK portfolio
#' res <- mvskPortfolio(M2 = M2, M3 = M3, M4 = M4, w0 = w0)
#'
#'
#' @export mvskPortfolio
mvskPortfolio <- function(m1 = NULL, M2 = NULL, M3 = NULL, M4 = NULL,
                          w0 = NULL, g = NULL, kappa = NULL, lb = NULL, ub = NULL,
                          lin_eq = NULL, lin_eqC = NULL, nlin_eq = NULL,
                          lin_ieq = NULL, lin_ieqC = NULL, nlin_ieq = NULL,
                          riskcriterion = NULL, options = list(), relative = TRUE) {

  p <- nrow(M2)

  # initial portfolio
  if (is.null(w0)) w0 <- "maxDiv"
  if (is.null(lin_eq) || is.null(lin_eqC)) {
    lin_eq <- matrix(1, nrow = 1, ncol = p)
    lin_eqC <- 1
  }
  if (!is.numeric(w0)) {
    initport <- solvePortfolio(p, w0, m1, M2, M3, M4, lb, ub, lin_eq, lin_eqC,
                               nlin_eq, lin_ieq, lin_ieqC, nlin_ieq, options)
  } else {
    initport <- NULL
  }

  # efficient update - for each value of kappa
  indmom <- !c(is.null(m1), is.null(M2), is.null(M3), is.null(M4))
  if (is.null(kappa)) kappa <- 1
  if (is.null(g)) g <- abs(getmom(indmom, initport$wopt, m1, M2, M3, M4))

  wopt <- matrix(NA, nrow = length(kappa), ncol = p)
  delta <- rep(NA, length(kappa))
  moms <- matrix(NA, nrow = length(kappa), ncol = sum(indmom))
  for (ii in 1:length(kappa)) {
    effport <- solveMVSKPortfolio(p, initport, kappa[ii], g, m1, M2, M3, M4, indmom, lb, ub, lin_eq,
                                  lin_eqC, nlin_eq, lin_ieq, lin_ieqC, nlin_ieq, options, relative)
    wopt[ii,] <- effport$wopt
    delta[ii] <- effport$delta
    moms[ii,] <- effport$moms
  }

  # select optimal value of kappa depending on some other criterium
  if (is.null(riskcriterion)) riskcriterion <- function(w) fERC(w, M2)$objective
  critvals <- apply(wopt, 1, riskcriterion)
  indopt <- which.min(critvals)

  summ_list <- list("w" = wopt, "kappa" = kappa, "delta" = delta, "moms" = moms,
                    "critvals" = critvals, "indopt" = indopt)
  w <- wopt[indopt,]

  return (list("w" = w, "summ" = summ_list))
}
