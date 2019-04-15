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
#' @param lb lower bound for the weights, default rep(0, p)
#' @param ub upper bound for the weights, default rep(1, p)
#' @param lin_eq equality constraints: eq w = eqC (should be matrix!), see details
#' @param lin_eqC equality constraints: eq w = eqC, see details
#' @param nlin_eq function with non-linear equality constraint (returns objective value and jacobian)
#' @param lin_ieq inequality constraints: ieq w leq ieqC (should be matrix!)
#' @param lin_ieqC inequality constraints: ieq w leq ieqC
#' @param nlin_ieq function with non-linear inequality constraints (returns objective value and jacobian)
#' @param href reference function to keep into account while tilting the portfolios (name or function)
#' @param kappa vector of values indicating maximum deviation of the reference objective
#' @param relative boolean indicating if kappa is absolute or relative
#' @param param list with extra arguments for the href function
#' @param options optimization options
#' @param mompref direction of preference for the moments, defaults to higher mean and skewnes,
#' lower variance and kurtosis
#' @author Dries Cornilly
#' @references
#' Boudt, K., Cornilly, D., Van Holle, F., & Willems, J. (2019). Algorithmic portfolio
#' tilting to harvest higher moment gains. working paper
#'
#' Briec, W., Kerstens, K., & Jokung, O. (2007). Mean-variance-skewness portfolio
#' performance gauging: a general shortage function and dual approach.
#' Management science, 53(1), 135-149.
#'
#' Choueifaty, Y., & Coignard, Y. (2008). Toward maximum diversification.
#' Journal of Portfolio Management, 35(1), 40.
#'
#' @examples
#' # load data
#' library(PerformanceAnalytics)
#' data(edhec)
#' x <- edhec[, 1:5]
#'
#' # estimate moments
#' m1 <- colMeans(x)
#' M2 <- cov(x)
#' M3 <- PerformanceAnalytics:::M3.MM(x, as.mat = FALSE)
#' M4 <- PerformanceAnalytics:::M4.MM(x, as.mat = FALSE)
#'
#' # optimal MVSK portfolio
#' resMVSK <- mvskPortfolio(m1 = m1, M2 = M2, M3 = M3, M4 = M4, w0 = "DR",
#'                          g = "mvsk", ub = rep(0.3, 5), href = "DR",
#'                          kappa = c(0, 0.01, 0.025, 0.05))
#'
#' # show weights
#' barplot(resMVSK$w, beside = TRUE)
#'
#' @export mvskPortfolio
mvskPortfolio <- function(m1 = NULL, M2 = NULL, M3 = NULL, M4 = NULL, w0 = NULL, g = NULL,
                          lb = NULL, ub = NULL, lin_eq = NULL, lin_eqC = NULL, nlin_eq = NULL,
                          lin_ieq = NULL, lin_ieqC = NULL, nlin_ieq = NULL, href = NULL, kappa = NULL,
                          relative = FALSE, param = NULL, options = list(), mompref = NULL) {

  p <- nrow(M2)

  # default constraints
  if (is.null(lin_eq) || is.null(lin_eqC)) {
    lin_eq <- matrix(1, nrow = 1, ncol = p)
    lin_eqC <- 1
  }
  if (is.null(lb)) lb <- rep(0, p)
  if (is.null(ub)) ub <- rep(1, p)

  # initial portfolio
  if (is.null(w0)) w0 <- rep(1 / p, p)
  if (!is.numeric(w0)) {
    w0 <- solvePortfolio(p, w0, m1, M2, M3, M4, lb, ub, lin_eq, lin_eqC,
                         nlin_eq, lin_ieq, lin_ieqC, nlin_ieq, options, param)
  }

  # initialize direction g of moment improvement
  indmom <- !c(is.null(m1), is.null(M2), is.null(M3), is.null(M4))
  if (is.null(g)) g <- abs(getmom(indmom, w0, m1, M2, M3, M4))
  if (is.character(g) && g[1] == "mvsk") {
    g <- abs(getmom(indmom, w0, m1, M2, M3, M4))
    if (indmom[1]) g[1] <- 0
  }

  # MVSK efficient portfolio
  if (is.null(href)) {
    # unrestricted MVSK efficient portfolio
    effport <- solveMVSKPortfolio(p, w0, g, m1, M2, M3, M4, indmom, lb, ub, lin_eq, lin_eqC, nlin_eq,
                                  lin_ieq, lin_ieqC, nlin_ieq, options, NULL, 0, FALSE, NULL, mompref)

  } else {
    # efficient portfolio with restriction on href
    wopt <- matrix(NA, nrow = length(kappa), ncol = p)
    delta <- rep(NA, length(kappa))
    moms <- matrix(NA, nrow = length(kappa), ncol = sum(indmom))
    constr <- NULL
    if (is.null(kappa)) kappa <- 0
    for (ii in 1:length(kappa)) {
      effport <- solveMVSKPortfolio(p, w0, g, m1, M2, M3, M4, indmom, lb, ub, lin_eq, lin_eqC,
                                    nlin_eq, lin_ieq, lin_ieqC, nlin_ieq, options,
                                    href, kappa[ii], relative, param, mompref)
      wopt[ii,] <- effport$w
      delta[ii] <- effport$delta
      moms[ii,] <- effport$moms
      constr <- rbind(constr, effport$constr)
    }
    effport <- list(w = wopt, delta = delta, moms = moms, constr = constr)
  }

  effport
}
