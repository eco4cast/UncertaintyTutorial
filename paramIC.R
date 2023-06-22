#Set number of samples
n = 100

#Parameters
library(mvtnorm)
paramSamples <- rmvnorm(n, fit$coef, fit$var.coef)

#IC (want first output from forecast)
ICuncert <- forecast(fit, h=0, level=0.68, xreg = as.matrix(site_target_past[, c('air_temperature', 'relative_humidity')]))
