#' Lambert W function
#' 
#' Implementation of the Lambert W function
#' 
#' @param x Value
#' @param tolerance Acceptable precision
#' @param maxit Maximum number of iterations to find the answer


f.lambertW <- function (x, tolerance = 1e-10, maxit = 50) 
{
  if (any(Im(x) != 0)) 
    stop("argument 'x' must be real, not complex!")
  ans <- x
  ans[!is.na(x) & x < -exp(-1)] <- NA
  ans[!is.na(x) & x >= -exp(-1)] <- log1p(x[!is.na(x) & x >= 
                                              -exp(-1)])
  ans[!is.na(x) & x >= 0] <- sqrt(x[!is.na(x) & x >= 0])/2
  cutpt <- 3
  if (any(myTF <- !is.na(x) & x > cutpt)) {
    L1 <- log(x[!is.na(x) & x > cutpt])
    L2 <- log(L1)
    wzinit <- L1 - L2 + (L2 + (L2 * (-2 + L2)/(2) + (L2 * 
                                                       (6 + L2 * (-9 + L2 * 2))/(6) + L2 * (-12 + L2 * 
                                                                                              (36 + L2 * (-22 + L2 * 3)))/(12 * L1))/L1)/L1)/L1
    ans[myTF] <- wzinit
  }
  for (ii in 1:maxit) {
    exp1 <- exp(ans)
    exp2 <- ans * exp1
    delta <- (exp2 - x)/(exp2 + exp1 - ((ans + 2) * (exp2 - 
                                                       x)/(2 * (ans + 1))))
    ans <- ans - delta
    if (all(is.na(delta)) || max(abs(delta), na.rm = TRUE) < 
        tolerance) 
      break
    if (ii == maxit) 
      warning("did not converge")
  }
  ans[x == Inf] <- Inf
  ans
}