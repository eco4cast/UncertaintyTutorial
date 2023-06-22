#Sobol calculation of indices for ARIMA model

library('sensitivity')
library('mvtnorm')

#Number of design samples
n <- 100
lag = 3

#Design X1 and X2
par1 <- rmvnorm(n, fit$coef, fit$var.coef)
par2 <- rmvnorm(n, fit$coef, fit$var.coef)

d1map <- as.matrix(sample(31, n, replace = TRUE))
colnames(d1map) <- c('noaa_ensemble_member')
d2map <- as.matrix(sample(31, n, replace = TRUE))
colnames(d2map) <- c('noaa_ensemble_member')

ICuncert <- forecast(fit, h=0, level=0.68, xreg = as.matrix(site_target_past[, c('air_temperature', 'relative_humidity')]))
IC.mean = head(ICuncert$mean,n=lag)
IC.sd   = (head(ICuncert$upper,n=lag) - head(ICuncert$lower,n=lag))/2 
ic1 <- rmvnorm(n,IC.mean,diag(as.vector(IC.sd)^2))
colnames(ic1) <- c('lag2', 'lag1', 'dayof')
ic2 <- rmvnorm(n, IC.mean, diag(as.vector(IC.sd)^2))
colnames(ic2) <- c('lag2', 'lag1', 'dayof')

x1 <- data.frame(cbind(par1, d1map, ic1))
x2 <- data.frame(cbind(par2, d2map, ic2))

#Sobol function call
x <- sensitivity::soboljansen(model = , X1 = x1, X2 = x2, conf = 0.95)
