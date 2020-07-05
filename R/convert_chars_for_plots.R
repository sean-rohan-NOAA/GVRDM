#' Convert model parameters to expressions for plotting
#' 
#' Function used in Rohan et al. GVRDM paper
#' @param x Data frame with a variable column that includes parameters to be converted to expressions
#' @return Data frame with character string that can be converted to an expression.

convert_chars_for_plots <- function(x) {
  # Transform parameters for plotting
  x$value[x$variable %in% c("sigma", "NVsigma", "sigma.rr")] <- exp(x$value[x$variable %in% c("sigma", "NVsigma", "sigma.rr")])
  
  x$Lwr_25[x$variable %in% c("sigma", "NVsigma", "sigma.rr")] <- exp(x$Lwr_25[x$variable  %in% c("sigma", "NVsigma", "sigma.rr")])
 
  x$Lwr_250[x$variable %in% c("sigma", "NVsigma", "sigma.rr")] <- exp(x$Lwr_250[x$variable %in% c("sigma", "NVsigma", "sigma.rr")])
 
  x$Upr_750[x$variable  %in% c("sigma", "NVsigma", "sigma.rr")] <- exp(x$Upr_750[x$variable  %in% c("sigma", "NVsigma", "sigma.rr")])
   
  x$Upr_975[x$variable  %in% c("sigma", "NVsigma", "sigma.rr")] <- exp(x$Upr_975[x$variable  %in% c("sigma", "NVsigma", "sigma.rr")])
    
  # Change parameter names for plotting
  x$variable <- as.character(x$variable)
  
  if(any(x$variable == "tt1")) {
    x$variable[x$variable == "tt"] <- "T[1]"
    x$variable[x$variable == "tt1"] <- "T[2]"
    x$variable[x$variable == "tt2"] <- "T[3]"
    x$variable[x$variable == "tt3"] <- "T[4]"
    x$variable[x$variable == "tt4"] <- "T[5]"
  }
  
  x$variable[x$variable == "Ke"] <- "K[e]"
  x$variable[x$variable == "kk"] <- "h"
  x$variable[x$variable == "NVrd"] <- "D[NV]"
  x$variable[x$variable == "NVsigma"] <- "sigma[NV]"
  x$variable[x$variable == "sigma"] <- "sigma[V]"
  x$variable[x$variable == "tt"] <- " T "
  x$variable[x$variable == "NVthreshold"] <- "q[NV]"
  x$variable[x$variable == "rr.max"] <- "r[max]"
  x$variable[x$variable == "sigma.rr"] <- "sigma"
  x$variable[x$variable == "kd.mult"] <- "m"
  x$variable[x$variable == "vv"] <- "v"
  x$variable[x$variable == "ww"] <- "w"
  x$variable[x$variable == "yy"] <- "y"
  


  x$variable <- factor(x$variable, levels = c(" T ", "T[1]", "T[2]", "T[3]", "T[4]", "T[5]",
                                              "K[e]", "sigma[V]", "D[NV]", "q[NV]", "sigma[NV]", 
                                              "h", "delta", "beta", 
                                              "r[max]", "sigma", "m", "v", "w", "y"))
  return(x)
}
