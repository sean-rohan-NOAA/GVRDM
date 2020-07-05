#' Dynamic scaling function
#' 
#' Function to plot a dynamic scaling function. Parameters correspond with model parameters.
#' 
#' @param beta Dynamic scaling function beta
#' @param delta Dynamic scaling function delta
#' @param kk Dynamic scaling function h
#' @param cc Effective attenuation coefficient or beam attenuation coefficient
#' @param return A vector of dynamic scaling function values.

dynamic_scaling <- function(beta, delta, kk, cc) {
  return(beta + dgamma(x = cc, shape = kk, scale = delta))
}
