
getmomgr <- function(indmom, w, m1, M2, M3, M4) {

  momsgr <- matrix(NA, nrow = 4, ncol = length(w))
  if (indmom[1]) momsgr[1,] <- m1
  if (indmom[2]) momsgr[2,] <- PerformanceAnalytics:::derportm2(w, M2)
  if (indmom[3]) momsgr[3,] <- PerformanceAnalytics:::derportm3(w, M3)
  if (indmom[4]) momsgr[4,] <- PerformanceAnalytics:::derportm4(w, M4)
  momsgr <- momsgr[indmom,]

  return (momsgr)
}
