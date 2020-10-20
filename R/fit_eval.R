#' Function to fit models
#' 
#' Uses model and predictors to generate fit. Also returns root means squared log error and R-squared.
#' 
#' @param data Data set as a list, passed in the same way as for an mle2 model.
#' @param mod mle2 model object.
#' @param Ke Optional. Default (NA) results in searching for the parameter in the model object.
#' @param alpha Optional. Default (NA) results in searching for the parameter in the model object.
#' @param sigma Optional. Default (NA) results in searching for the parameter in the model object.
#' @param tt Optional. Default (NA) results in searching for the parameter in the model object.
#' @param NVthreshold Optional. Default (NA) results in searching for the parameter in the model object.
#' @param silent Should diagnostic plots be generated? Passed to diagnostic_plots function.
#' @return A list with the input, model fits, R-squared, and root mean square log error


fit_eval <- function(dat, mod, Ke = NA, alpha = NA, sigma = NA, NVthreshold = NA, silent = F) {
  parid <- table(c("tt", "Ke", "sigma", "gamma0", "beta", "kk", "delta", "NVrd", "NVthreshold", "NVsigma", "cc.mult", "tt1", "tt2", "tt3", "tt4", "alpha", "angle", "kd", "ccoffset", "kd.mult"))
  parid[1:length(parid)] <- NA

  parid[match(names(mod@fullcoef), names(parid))] <- mod@fullcoef
  
  
  if("prey.types" %in% names(mod@data)) {
    prey.types <- dat$prey.types
    prey <- dat$prey
  } else {
    prey.types <- NA
    prey <- NA
  }
  
  if(is.na(alpha)) {
  alpha <- parid["alpha"]
  if(is.na(alpha)) {
    alpha <- 1
  }
  }
  
  if(is.na(Ke)) {
    Ke <- parid["Ke"]
  }
  
  if(is.na(sigma)) {
    sigma <- parid["sigma"]
  }

  # Fill out data list if kd and angle aren't in dat 
  if(!("kd" %in% names(dat))) {
    dat$kd <- NA
  }
  
  if(!("angle" %in% names(dat))) {
    dat$angle <- NA
  }
  
  if(is.null(dat$Ke.trans)) {
    warning("No Ke.trans found in dat. Setting Ke.trans to FALSE")
    dat$Ke.trans <- FALSE
  }
  
  mod_fit <- NLL_dynamic_contrast(rr = dat$rr, 
                       cc = dat$cc,
                       tt = parid["tt"],
                       tt1 = parid["tt1"],
                       tt2 = parid["tt2"],
                       tt3 = parid["tt3"],
                       tt4 = parid["tt4"],
                       alpha = alpha, 
                       Ke = Ke,
                       Eb = dat$Eb,
                       angle = dat$angle,
                       kd = dat$kd,
                       kd.mult = parid["kd.mult"],
                       Ke.trans = dat$Ke.trans,
                       beta = parid["beta"],
                       kk = parid["kk"],
                       delta = parid["delta"],
                       sigma = sigma,
                       NVrd = parid["NVrd"],
                       NVthreshold = parid["NVthreshold"],
                       NVsigma = parid["NVsigma"],
                       prey.types = prey.types,
                       ccoffset = parid["ccoffset"],
                       prey = prey,
                       fit.model = F,
                       rr.range = dat$rr.range,
                       silent = silent)
  
  mod_fit$R2 <- cor(log(mod_fit$rr), log(mod_fit$out))^2
  mod_fit$log_rmse <- rmse_rd(obs = mod_fit$rr, fit = mod_fit$out, log.trans = TRUE)
  
  return(mod_fit)
}
