#' Function for fitting the GVRDM and Aksnes and Utne model
#' 
#' Generalized visual reaction distance model and Aksnes and Utne model.
#' 
#' @param rr Reaction distance in meters
#' @param cc Effective attenuation coefficient or beam attenuation coefficient
#' @param Ke Half-saturation constant
#' @param Ke.trans Should Ke be transformed in the model?
#' @param Eb # Light level
#' @param Ap Prey area. If Ap and C0 are not provided, tt will be estimated.
#' @param C0 Prey inherent contrast. If Ap and C0 are not provided, tt will be estimated.
#' @param Eprime Composite saturation parameter.
#' @param tt T parameter for prey type 0. Default = NA
#' @param tt1 T parameter for prey type 1. Default = NA
#' @param tt2 T parameter for prey type 2. Default = NA
#' @param tt3 T parameter for prey type 3. Default = NA
#' @param tt4 T parameter for prey type 4. Default = NA
#' @param angle Nadir viewing angle in degrees
#' @param kd Diffuse attenuation coefficient of downwelling irradiance
#' @param beta Dynamic scaling function intercept parameter
#' @param hh Dynamic scaling function rate parameter
#' @param delta Dynamic scaling function shape parameter
#' @param alpha Naka-Rushton exponent. Default = 1
#' @param sigma Standard deviation of visual reaction distance
#' @param NVrd Non-visual reaction distance. Default = NA (only fit visual component of the model)
#' @param NVthreshold Light threshold for non-visual reaction. Default = NA(only fit visual component of the model)
#' @param NVsigma Standard deviation of non-visual reaction distance. Default = NA(only fit visual component of the model)
#' @param kd.mult Multiplier to convert beam attenuation 
#' @param prey.types A vector specifying unique prey name categories that are passed to 'prey.' Used if parameters are simultaneously estimated for multiple prey types.
#' @param prey Vector of prey names for each observation.
#' @param prey.size Vector of prey size
#' @param ccoffset Shift kappa or beam attenuation? Default = NA results in no shift.
#' @param rr.log Are visual errors lognormal distribution? Default = TRUE
#' @param fit.model Logical.Should the model be fitted? If FALSE, produces model diagnostics and predictions based on provided parameters.
#' @param fit.obs Logical. Should observations be used to calculate a likelihood?
#' @param silent Logical. Should diagnostic plots be produced?
#' @param return If fit.model is TRUE, returns the negative log likelihood. If fit is FALSE, returns model diagnostics and predictions using parameters that are passed to the model.
#' @details The order of composite tt parameters should match the order of prey.types in the model. By default, the model is set up to estimate five prey parameters (tt, tt1, tt2, tt3, tt4). Additional parameters can be added by modifying the function with additional numbered parameters that have names beginning with tt. Regardless of how many tt parameters are added, only parameters that are assigned values and have a corresponding prey category will be estimated by the model.

fit_gvrdm <- function(rr,
                                 cc,
                                 Ap = NA,
                                 C0 = NA,
                                 Eprime = NA,
                                 prey.types = NA,
                                 prey = NA,
                                 prey.size = NA,
                                 tt = NA,
                                 tt1 = NA,
                                 tt2 = NA,
                                 tt3 = NA,
                                 tt4 = NA,
                                 Ke,
                                 Ke.trans = FALSE,
                                 Eb,
                                 angle = NA,
                                 kd = NA,
                                 kd.mult = NA,
                                 beta = NA,
                                 kk = NA,
                                 delta = NA,
                                 alpha = 1,
                                 sigma = NA,
                                 rr.log = TRUE,
                                 NVrd = NA,
                                 NVthreshold = 0,
                                 NVsigma = NA,
                                 fit.model = TRUE,
                                 fit.obs = TRUE,
                                 ccoffset = NA,
                                 silent = F, ...) {
  # Initialize output vectors
  rhs <- vector(length = length(cc))
  lhs <- vector(length = length(cc))
  out <- vector(length = length(cc))
  VIS <- 1:length(cc)
  
  # Transform Ke?
  if(Ke.trans == TRUE) {
    Ke <- 10^Ke
  }
  
  # Offset beam attenuation?
  if(is.na(ccoffset)) {
    ccoffset <- 0
  }
  
  # Angular dependence?
  if(is.na(angle[1])) {
    kd <- rep(0, length(rhs))
    angle <- rep(90, length(rhs))
  }
  
  # Diffuse downwelling attenuation coefficient multiplier?
  if(!is.na(kd.mult)) {
    kd <- kd*kd.mult
  }
  
  # Split visual and Non-visual reaction distance
  if(!is.na(NVrd)) {
    VIS <- which(Eb > NVthreshold)
  }
  
  # Non-visual reaction distance
  out[-VIS] <- NVrd
  
  # Assign separate tt for each prey type
  # Multiple prey types
  if(length(prey.types) > 1) {
    tt_vars <- ls()[grepl("tt", ls())]
    tt_vars <- sapply(tt_vars, function(x) eval(parse(text=x)))
    tt <- rep(tt, length(rr))
    tt_vec <- match(prey, prey.types)
    tt_vec <- tt_vars[tt_vec]
    tt <- tt_vec
  } else {
    tt <- rep(tt, length(rr))
  }
  
  # Size multiplier
  if(!is.na(prey.size[1])) {
    tt <- tt * prey.size^2
  }
  
  # Lambert W function
  lambert.fn <- function(cc, xx) {
    return((2/cc)*f.lambertW(cc*sqrt(xx)/2))
  }
  
  if(!is.na(C0)) {
    C0 <- trans.atan(val = C0, a = 0, b = 1)
  }
  
  # Left-hand side
  if(fit.obs) {
    lhs[VIS] <- 2 * log(rr[VIS]) + (cc[VIS]-kd[VIS]*round(cos(angle[VIS]*pi/180), 3))*rr[VIS]
  }
  
  # Start rhs
  rhs[VIS] <- Eb[VIS]^alpha/(Eb[VIS]^alpha + Ke^alpha)
  
  
  if(is.na(Ap) & is.na(C0)) {
    if(is.na(Eprime)) {
      # Aksnes and Utne model
      rhs[VIS] <- tt[VIS]*rhs[VIS]
    } else {
      # Cases with prey type variation
    }
  }
  
  if(is.na(beta) & is.na(kk)) {
    # Aksnes and Utne model - Do nothing
  } else {
    if(is.na(beta)) {
      cont_shape <- dgamma(x = abs(cc[VIS]-ccoffset), shape = kk, scale = delta)
    } else {
      cont_shape <- (beta + dgamma(x = abs(cc[VIS]-ccoffset), shape = kk, scale = delta))
    }
    
    rhs[VIS] <- rhs[VIS]*cont_shape
    
    if(fit.model == T) {
      # Avoid cont_shape being out of bounds
      if(any(is.na(cont_shape))) {
        return(1e7)
      } else if(any(cont_shape <= 0)) {
        return(1e7)
      } 
    }
  }
  
  # Lambert W function to find root
  out[VIS] <- lambert.fn(cc = (cc[VIS]-kd[VIS]*round(cos(angle[VIS]*pi/180), 3)), xx = rhs[VIS])
  
  # NLL for visual
  if(!is.na(sigma)) {
    if(rr.log) { # Log-transformed?
      sigma <- exp(sigma) # Must >0
      NLL <- -1*sum(log(dnorm(log(rr[VIS]), mean = log(out[VIS]), sd = sigma)))
    } else {
      NLL <- -1*sum(log(dnorm(rr[VIS], mean = out[VIS], sd = sigma)))
    }
  }
  # NLL for nonvisual
  if(!is.na(NVrd)) {
    if(is.na(NVsigma)) {
      NVsigma <- sigma
    } else {
      NVsigma <- exp(NVsigma)
    }
    NLL <- NLL + -1*sum(log(dnorm(log(rr[-VIS]), mean = log(NVrd), sd = NVsigma)))
    
  }
  
  
  if(is.infinite(NLL) | is.na(NLL)) {
    NLL <- 1e32
  }
  
  if(fit.model) {
    return(NLL)
  } else {
    if(fit.obs) {
      if(!is.na(NVrd)) {
        out[Eb <= NVthreshold] <- NVrd
      }
      if(!silent) {
        out <- list(out = out, 
                    vis.diagnostic.plots = diagnostic_plots(rr = rr[VIS], cc = cc[VIS], Eb = Eb[VIS], out = out[VIS], cont_shape = cont_shape[VIS], sigma = sigma, NVrr = rr[-VIS], NVrd = NVrd, NVsigma = NVsigma))
      } else {
        out <- list(out = out)
      }
    } else {
      if(!is.na(NVrd)) {
        out[Eb <= NVthreshold] <- NVrd
      }
      prediction_plots(cc = cc, Eb = Eb, out = out)
    }
    
    if(!is.na(NVrd)) {
      out$Eb <- Eb
      out$cc <- cc
      out$rr <- rr
      out$angle <- angle
      out$kd <- kd
    } else {
      out$Eb <- Eb
      out$cc <- cc
      out$rr <- rr      
      out$angle <- angle
      out$kd <- kd
    }
    
    return(out)
  }
  
}