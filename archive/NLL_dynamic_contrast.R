# Function for fitting the GVRDM and Aksnes and Utne model
# Sean K. Rohan <sean.rohan@noaa.gov>
# Last update: January 22, 2020

NLL_dynamic_contrast <- function(rr, # Distance
                                 cc, # Beam attenuation
                                 Ap = NA, # Prey area
                                 C0 = NA, # Inherent contrast
                                 Eprime = NA, # E' parameter
                                 prey.types = NA, # A vector with named prey types
                                 prey = NA,
                                 prey.size = NA,
                                 tt = NA, # T parameter for prey type 0
                                 tt1 = NA, # T parameter for prey type 1
                                 tt2 = NA, # T parameter for prey type 2
                                 tt3 = NA, # T parameter for prey type 3
                                 tt4 = NA, # T parameter for prey type 4
                                 Ke, # Half saturation constant
                                 Eb, # Light
                                 angle = NA, # Sagittal attack angle
                                 kd = NA, # Downwelling attenuation
                                 kd.mult = NA,
                                 beta = NA, # Dynamic contrast beta
                                 kk = NA, # Dynamic contrast k
                                 delta = NA, # Dynamic contrast delta
                                 alpha = 1, # Naka-Rushton exponent
                                 sigma = NA, # Standard deviation
                                 rr.log = TRUE, # Should the MLE for all values of rr be done in log space
                                 NVrd = NA, # Non-visual reaction distance threshold
                                 NVthreshold = 0, # Non-visual light threshold
                                 NVsigma = NA, # Non-visual sigma
                                 cont_shape = NA, # Deprecated
                                 fit.model = TRUE, # Logical: should the model be fitted or plotted
                                 fit.obs = TRUE,
                                 ccoffset = NA,
                                 silent = F, # 
                                 ...) {  # Estimated parameter
  # Initialize output vectors
  rhs <- vector(length = length(cc))
  lhs <- vector(length = length(cc))
  out <- vector(length = length(cc))
  VIS <- 1:length(cc)
  
  # Offset beam attenuation?
  if(is.na(ccoffset)) {
    ccoffset <- 0
  }
  
  # Angular dependence?
  if(is.na(angle[1])) {
    kd <- rep(0, length(rhs))
    angle <- rep(0, length(rhs))
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
  if(is.na(tt1)) {
    tt <- rep(tt, length(rr))
  } else {
    tvec <- rep(0, length(prey))
    tvec[prey == prey.types[1]] <- tt
    tvec[prey == prey.types[2]] <- tt1
    if(!is.na(tt2)) {
      tvec[prey == prey.types[3]] <- tt2
      if(!is.na(tt3)) {
        tvec[prey == prey.types[4]] <- tt3
      }
      if(!is.na(tt4)) {
        tvec[prey == prey.types[5]] <- tt4
      }
    }
    tt <- tvec
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
    lhs[VIS] <- 2 * log(rr[VIS]) + (cc[VIS]-kd[VIS]*round(cos(angle[VIS]), 3))*rr[VIS]
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
      #cont_shape <- rhs[VIS] * dgamma(x = abs(cc[VIS]-ccoffset), shape = kk, scale = delta)
      cont_shape <- dgamma(x = abs(cc[VIS]-ccoffset), shape = kk, scale = delta)
    } else {
      #cont_shape <- rhs[VIS] * (beta + dgamma(x = abs(cc[VIS]-ccoffset), shape = kk, scale = delta))
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
  out[VIS] <- lambert.fn(cc = (cc[VIS]-kd[VIS]*round(cos(angle[VIS]), 3)), xx = rhs[VIS])
  
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
      }
      if(!silent) {
        out <- list(out = out, 
                    vis.diagnostic.plots = diagnostic_plots(rr = rr[VIS], cc = cc[VIS], Eb = Eb[VIS], out = out[VIS], cont_shape = cont_shape[VIS], sigma = sigma, NVrr = rr[-VIS], NVrd = NVrd, NVsigma = NVsigma))
      } else {
        out <- list(out = out)
      }
    } else {
      if(!is.na(NVrd)) {
        out[out < NVrd & Eb <= NVthreshold] <- NVrd
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
