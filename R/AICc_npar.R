#' Calculate AICc using a user-specified number of parameters
#'
#' This function is used to estimate AICc in cases where it is necessary to pass estimated fixed values of parameters into a model in order to estimate other parameters. Arises when using step functions with a discontinuous likelihood profile.
#'
#' @param object mle2 model object
#' @param npar Number of parameters in the model. If NA, will use the length of fullcoef(object)
#' @param use.AIC Use AIC instead of AICc?
#' @param npar.out Returns number of parameters and AICc.
#'
#' @author S.K. Rohan \email{skrohan@@uw.edu}

AICc_npar <- function(object, npar = NA, npar.out = FALSE) {
  nn <- length(object@data[[1]])

  if(is.na(npar)) {
    npar <- length(object@fullcoef)
  }
  out <- as.numeric((2*npar - 2*logLik(object)) + (2*npar^2 + 2*npar)/(nn - npar -1))
  if(npar.out){
    out <- c(AICc = out, df = npar)
  }

  return(out)

}
