#One-at-a-time analysis for ARIMA model

library(tidyverse)
library(distributional)

#reset the matrix
A <- data.frame(cbind(par1, d1map, ic1,sigma.ens1))

#forecast with the ensemble
y1 <- arima.fx(x1,drivers,epsilon,horiz=35)


#Creates matrices that remove uncertainty by replacing with means
A_means <- x1 %>% 
              mutate(air_temperature = mean(x1$air_temperature)) %>% 
              mutate(relative_humidity = mean(x1$relative_humidity)) %>% 
              mutate(intercept = mean(x1$intercept)) %>% 
              mutate(lag1 = mean(x1$lag1)) %>% 
              mutate(lag2 = mean(x1$lag2))  %>% 
              mutate(dayof = mean(x1$dayof)) 

driver_means <- drivers %>% 
              mutate(air_temperature = mean(drivers$air_temperature)) %>% 
              mutate(surface_downwelling_shortwave_flux_in_air = mean(drivers$surface_downwelling_shortwave_flux_in_air)) %>% 
              mutate(precipitation_flux = mean(drivers$precipitation_flux)) %>% 
              mutate(relative_humidity = mean(drivers$relative_humidity)) 

epsilon.none = 0*epsilon


#forecast with no uncertainty

y.none = arima.fx(A_means,driver_means,epsilon.none,horiz=35)
var.none = apply(y.none,2,var,na.rm=TRUE)
plot(var.none)


#replace variable mean with data to test for parameter (air temp)
A_AirTemp <- A_means %>% mutate(air_temperature = x1$air_temperature)

#forecast with uncertainty in one parameter (air temp)

y.none = arima.fx(A_means,driver_means,epsilon.none,horiz=35)
var.none = apply(y.none,2,var,na.rm=TRUE)
plot(var.none) 


