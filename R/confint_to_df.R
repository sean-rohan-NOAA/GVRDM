# Make output of confint into a data frame
# Sean K. Rohan <skrohan@uw.edu>
# Last update: January 21, 2020

confint_to_df <- function(conf.dat, Model) {
  out <- as.data.frame(conf.dat)
  names(out) <- c("Lwr_250", "Upr_750", "Lwr_25", "Upr_975")
  out$variable <- rownames(out)
  out$Model <- Model
  return(out)
}