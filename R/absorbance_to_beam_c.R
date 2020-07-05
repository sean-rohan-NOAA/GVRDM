#' Convert absorbance measurements to beam attenuation
#' 
#' Function converts absorbance measurements obtained using spectophotometer to beam attenuation based on the acceptance angle of the spectrophotometer, turbidity standard, measurement units for turbidity, cuvette length, and transmittance wavelength.
#' 
#' @param ANGLE Acceptance angle of the spectrophotometer
#' @param TURB Turbidity measurement in units of JTU, NTU
#' @param L Length of spectrophotometer cuvette
#' @param TRANSMITTANCE Percentage of transmitted light reaching sensor
#' @param UNITS Turbidity measurement units
#' @param WAVELENGTH Transmitted wavelength for spectrophotometer

absorbance_to_beam_c <- function(ANGLE, 
                                 TURB, 
                                 L, 
                                 TRANSMITTANCE, 
                                 UNITS = "JTU", 
                                 WAVELENGTH = 800){
  if(UNITS == "JTU" & WAVELENGTH == 800) {
    B_theta <- 50.93*TURB + 39 # Zaneveld et al. 1979
    return(2*pi*B_theta*(1-cos(ANGLE*pi/180))-log(TRANSMITTANCE)/L) # Aksnes and Utne (1997)    
  }
}