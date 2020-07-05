#' Modified Naka-Rushton function
#' 
#' Contrast function (Naka-Rushton) with the flexibility to fit additional shapes, as described in Peirce (2007). Either of the function arguments nn or ss can be considered the equivalent of the alpha parameter.
#' 
#' @param Rmax Maximum retinal response (default = 1)
#' @param Eb Incident light
#' @param Ke Half-saturation parameter
#' @param nn Shape parameter (Peirce 2007)
#' @param ss Shape parameter 
#' @param BB Baseline signal without any incident light.
#' @return Returns a vector of values of the Naka-Rushton function.
#' @references Peirce, J.W., 2007. The potential importance of saturating and supersaturating contrast response functions in visual cortex. J. Vis. 7, 1–10. https://doi.org/10.1167/7.6.13

naka_rushton <- function(Rmax = 1, Eb, Ke, nn = 1, ss = 1, BB = 0) {
  return(Rmax * Eb^nn /(Ke^(nn*ss) + Eb^(nn*ss)) + BB)
}