#' Broken stick reaction distance model 
#'
#' Continuous version modified from Hansen et al. (2013, 2015)
#' 
#' @param rr Observed reaction distance
#' @param Eb Light intensity
#' @param cc Effective attenuation coefficient, beam attenuation coefficient, or NTU
#' @param vv Light function intercept
#' @param ww Light function slope
#' @param yy Turbidity function slope
#' @param Eb.break Saturation intensity threshold for light
#' @param cc.break Minimum turbidity threshold
#' @param rr.max Maximum reaction distance
#' @param sigma.rr Standard deviation of reaction distance
#' @param fit Use the model in fitting mode?
#' @return If fit == TRUE, returns the negative log likelihood. Otherwise, returns a list containing model parameters, negative log-likelihood, and fitted reaction distances.

hansen_broken_stick <- function(rr, # Reaction distance
                                Eb, # Light intensity
                                cc, # Beam attenuation
                                vv, # Light function intercept
                                ww, # Light function slope
                                # xx = NA, 
                                yy, # Turbidity function slope
                                Eb.break, # Saturation intensity threshold
                                cc.break, # Minimum turbidity threshold
                                rr.max, # Maximum reaction distance
                                sigma.rr = NULL, # Standard deviation
                                fit = TRUE) {
  rr.fit = vector(length = length(rr))
  
  # Light function
  rr.fit[Eb <= Eb.break] <- ww * Eb[Eb <= Eb.break] + vv
  rr.fit[Eb > Eb.break] <- rr.max
  
  # Turbidity function 
  rr.fit[cc > cc.break] <- rr.fit[cc > cc.break] * exp(yy * (cc.break- cc[cc > cc.break]))
  
  nll <- sum(-1*dnorm(log(rr.fit), log(rr), exp(sigma.rr), log = TRUE))
  if(!fit) {
    return(list(out = rr.fit, 
                nll = nll, 
                rr = rr, 
                Eb = Eb, 
                cc = cc, 
                ww = ww, 
                yy = yy, 
                Eb.break = Eb.break, 
                cc.break = cc.break, 
                rr.max = rr.max, 
                sigma.rr = sigma.rr))
  } else {
    if(is.na(nll)) {
      return(1e32)
    } else {
      return(nll)
    }
  }
}