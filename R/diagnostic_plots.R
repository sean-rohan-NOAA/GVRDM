#' Diagnostic plots for the AUM and GVRDM
#' 
#' Function called by fit_gvrdm to generate model diagnostic plots.
#' 
#' @param rr Observed reaction distances for the visual component
#' @param out Fitted reaction distance for the visual component
#' @param cc Effective attenuation coefficient or beam attenuation coefficient
#' @param Eb light
#' @param sigma standard deviation for the visual component
#' @param cont_shape Optional. Manually-specified dynamic scaling function values.
#' @param NVrd Observed non-visual reaction distances for the non-visual component
#' @param NVrr Fitted non-visual reaction distances for the non-visual component
#' @param NVsigma Standard deviation for non-visual reaction distance
#' @return Returns a list containing diagnostic plots.

diagnostic_plots <- function(rr = NA, out = NA, cc = NA, Eb = NA, sigma = NA, cont_shape = NA, NVrd = NA, NVrr = NA, NVsigma = NA) {
  qqplot.data <- function(vec) {
    y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    
    d <- data.frame(resids = vec)
    
    return(ggplot(d, aes(sample = resids)) + stat_qq() + 
             geom_abline(slope = slope, intercept = int) + 
             scale_x_continuous(name = "Theoretical")
           + scale_y_continuous(name = "Sample") + theme_bw())
    
  }
  RD.stresid <- log(rr)-log(out)
  RD.stresid <- rr-out
  
  p0 <- qqplot.data(RD.stresid)
  
  p1 <- ggplot() + 
    geom_point(aes(x = out, y = RD.stresid))+
    scale_x_continuous(name = expression(hat(r))) + 
    scale_y_continuous(name = expression(ln(r)-ln(hat(r)))) + 
    theme_bw()
  
  p2 <- ggplot() + 
    geom_point(aes(x = Eb, y = RD.stresid)) +
    geom_smooth(aes(x = Eb, y = RD.stresid)) +
    scale_y_continuous(name = expression(ln(r)-ln(hat(r)))) + 
    scale_x_log10(name = expression(E[b])) + theme_bw()
  
  p3 <- ggplot() + 
    geom_point(aes(x = cc, y = RD.stresid)) +
    geom_smooth(aes(x = cc, y = RD.stresid)) +
    scale_y_continuous(name = expression(ln(r)-ln(hat(r)))) + 
    scale_x_continuous(name = "c") + theme_bw()
  
  p4 <- ggplot() + 
    geom_histogram(aes(x = RD.stresid), bins = 10) +
    scale_x_continuous(name = expression(ln(r)-ln(hat(r)))) + 
    scale_y_continuous(name = "Frequency") + 
    theme_bw()
  
  p5 <- ggplot() + 
    geom_point(aes(x = rr, y = out, color = cc), size = rel(2)) + 
    scale_x_continuous(name = "r") + 
    scale_y_continuous(name = expression(hat(r))) + 
    scale_color_continuous(name = "c") + 
    theme_bw() + 
    theme(legend.position = "bottom")
  
  if(is.numeric(NVrr)) {
    rr.vals <- c(NVrr, rr)
  } else{
    rr.vals <- rr
  }
  
  fit_by_light <- data.frame(rd = c(rr, out), 
                             cc = cc, 
                             Eb = Eb, 
                             Type = c(rep("Observed", length(rr)), rep("Fitted", length(out)))) 
  
  p6 <- ggplot() + 
    geom_point(data = fit_by_light, aes(x = cc, y = rd, color = Eb, shape = Type)) + 
    scale_x_continuous(name = "c") +
    scale_y_continuous(name = "r") + 
    scale_color_continuous(name = expression(E[b])) + theme_bw() + theme(legend.position = "bottom")
  
  if(!is.na(NVrd)) {
    NVRD.resid <- NVrr - NVrd
    nv0 <- qqplot.data(NVRD.resid)
    nv1 <- ggplot() + 
      geom_histogram(aes(x = NVRD.resid), bins = 10) + 
      scale_x_continuous(name = "NV Residual") + 
      scale_y_continuous(name = "Frequency") + 
      theme_bw()
    print(grid.arrange(nv0, nv1, ncol = 2))
  }
  
  if(!any(is.na(cont_shape))) {
    p7 <- ggplot() + 
      geom_point(aes(x = cc, y = cont_shape)) + 
      scale_x_continuous(name = "c") + 
      scale_y_continuous(name = expression(omega(c))) + 
      theme_bw()
    print(grid.arrange(p0, p4, p2, p3, nrow = 2, ncol = 2))
    print(p1)
    print(p5)
    print(p6)
    
    if(!is.na(NVrd)) {
      return(list(p0 = p0, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6, nv0 = nv0, nv1 = nv1))
    } else {
      return(list(p0 = p0, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6))
    }
    
  } else {
    print(grid.arrange(p0, p4, p2, p3, nrow = 2, ncol = 2))
    print(grid.arrange(p1, p5, p6, ncol = 2, nrow = 2))
    
    if(!is.na(NVrd)) {
      return(list(p0 = p0, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6, nv0 = nv0, nv1 = nv1))
    } else {
      return(list(p0 = p0, p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6))
    }
    return(list(p0 = p0,p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6))
  }
  
}