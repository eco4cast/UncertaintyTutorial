#Set number of samples
ne = 100
lag = 3

#Parameters
library(mvtnorm)
paramSamples <- rmvnorm(ne, fit$coef, fit$var.coef)

#IC (want first output from forecast)
ICuncert <- forecast(fit, h=0, level=0.68, xreg = as.matrix(site_target_past[, c('air_temperature', 'relative_humidity')]))
IC.mean = tail(ICuncert$mean,n=lag)
IC.sd   = (tail(ICuncert$upper,n=lag) - tail(ICuncert$lower,n=lag))/2 
IC      = rmvnorm(ne,IC.mean,diag(as.vector(IC.sd)^2))
