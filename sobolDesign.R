#Sobol calculation of indices for ARIMA model

library('sensitivity')
library('mvtnorm')

#Number of design samples
n <- 100000
ne <- 100000
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
ic1 <- rmvnorm(n,as.vector(IC.mean),diag(as.vector(IC.sd)^2,nrow = lag))
colnames(ic1) <- tail(c('lag2', 'lag1', 'dayof'),n = lag)
ic2 <- rmvnorm(n, as.vector(IC.mean), diag(as.vector(IC.sd)^2, nrow = lag))
colnames(ic2) <- tail(c('lag2', 'lag1', 'dayof'), n =lag)

sigma.ens1 = data.frame(sigma.ens = 1:ne)
sigma.ens2 = data.frame(sigma.ens = sample(1:ne,ne,replace=TRUE))

x1 <- data.frame(cbind(par1, d1map, ic1, sigma.ens1))
x2 <- data.frame(cbind(par2, d2map, ic2, sigma.ens2))

epsilon = rmvnorm(ne,rep(0,horiz+lag),diag(rep(fit$sigma2,horiz+lag)))

#Sobol function call
sobolOuts <- soboljansen(model = arima.fx, X1 = x1, X2 = x2, drivers = noaa_future_daily, horiz = 35, epsilon = epsilon, conf = 0.95)

#Clean up negatives (set to 0)
sobolOuts$S[sobolOuts$S < 0] <- 0
sobolOuts$T[sobolOuts$T < 0] <- 0

#Plots
plot(sobolOuts$S[1,], type = 'l', ylim = c(0,1.05), lwd = 3, main = 'Sobol First Order Indices',
     ylab = 'Index', xlab = 'Day')
for(i in 2:nrow(sobolOuts$S)){lines(sobolOuts$S[i,], ylim = c(0,1.1), col = i, lwd = 3)}
legend('topright', legend = rownames(sobolOuts$S), col = 1:7, lty=1, lwd=3, cex = 0.8)

plot(sobolOuts$T[1,], type = 'l', ylim = c(0,1.05), lwd = 3, main = 'Sobol Total Order Indices',
     ylab = 'Index', xlab = 'Day')
for(i in 2:nrow(sobolOuts$T)){lines(sobolOuts$S[i,], ylim = c(0,1.1), col = i, lwd = 3)}
legend('topright', legend = rownames(sobolOuts$T), col = 1:7, lty=1, lwd=3, cex = 0.8)
