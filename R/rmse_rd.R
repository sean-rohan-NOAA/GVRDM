#' Root mean square error lor root mean square log error
#' 
#' Function to calculate the root mean square error (RMSE) or root mean square log error (RMSLE).
#' 
#' @param obs Vector of observed values
#' @param fit Vector of fitted values
#' @param log.trans Should observed and fitted vectors be log-transformed to calculate RMSLE?
#' @return Vector of RMSE or RMSLE.


rmse_rd <- function(obs, fit, log.trans = FALSE) {
  if(log.trans) {
    obs <- log(obs)
    fit <- log(fit)
  }  
  return(sum(sqrt((obs-fit)^2))/length(obs))
    
}