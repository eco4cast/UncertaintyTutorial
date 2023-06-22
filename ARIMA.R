# tg_arima model
# written by ASL, 21 Jan 2023



#### Step 0: load packages

library(tidyverse)
library(neon4cast)
library(lubridate)
library(rMR)
library(glue)
#source("ignore_sigpipe.R")
library(tsibble)
library(fable)
library(arrow)
#source("download_target.R")
library(forecast)
library(here)

download_target <- function(theme = c("aquatics", "beetles",
                                      "phenology", "terrestrial_30min",
                                      "terrestrial_daily","ticks")){  
  theme <- match.arg(theme)
  
  target_file <- switch(theme,
                        aquatics = "aquatics-targets.csv.gz",
                        beetles = "beetles-targets.csv.gz",
                        phenology = "phenology-targets.csv.gz",
                        terrestrial_daily = "terrestrial_daily-targets.csv.gz",
                        terrestrial_30min = "terrestrial_30min-targets.csv.gz",
                        ticks = "ticks-targets.csv.gz"
  )
  download_url <- paste0("https://data.ecoforecast.org/neon4cast-targets/",
                         theme, "/", target_file)
  
  readr::read_csv(download_url, show_col_types = FALSE,
                  lazy = FALSE, progress = FALSE)#%>% 
  #as_tibble(index=time, key=siteID)
}
target = download_target(theme="phenology")
target = target |> filter(site_id == "HARV", variable == "gcc_90")

#### Step 1: Define team name, team members, and theme
model_id = "tg_arima"
model_themes = c("phenology") #This model is only relevant for three themes. I am registered for all three
model_types = c("phenology") #Replace terrestrial daily and 30min with terrestrial
#Options: aquatics, beetles, phenology, terrestrial_30min, terrestrial_daily, ticks


#### Step 2: Get NOAA driver data
forecast_date <- as.Date("2023-05-04")
noaa_date <- forecast_date - lubridate::days(1)  #Need to use yesterday's NOAA forecast because today's is not available yet

#We're going to get data for all sites relevant to this model, so as to not have to re-load data for the same sites
site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") %>%
  filter(if_any(matches(model_types),~.==1))
all_sites = site_data$field_site_id

# specify meteorological variables needed to make predictions
variables <- c('air_temperature',
               "surface_downwelling_shortwave_flux_in_air",
               "precipitation_flux",
               "relative_humidity")

# Load stage 2 data
endpoint = "data.ecoforecast.org"
use_bucket <- paste0("neon4cast-drivers/noaa/gefs-v12/stage2/parquet/0/", noaa_date)
use_s3 <- arrow::s3_bucket(use_bucket, endpoint_override = endpoint, anonymous = TRUE)
noaa_future <- arrow::open_dataset(use_s3) |>
  dplyr::filter(site_id %in% 'HARV',
                datetime >= forecast_date,
                #                reference_datetime == lubridate::as_datetime(forecast_date),
                variable %in% variables) |>
  dplyr::collect()

# Format met forecasts
noaa_future_daily <- noaa_future |> 
  mutate(datetime = lubridate::as_date(datetime)) |> 
  # mean daily forecasts at each site per ensemble
  group_by(datetime, parameter, variable) |> 
  summarize(prediction = mean(prediction)) |>
  pivot_wider(names_from = variable, values_from = prediction) |>
  # convert to Celsius
  mutate(air_temperature = air_temperature - 273.15) |> 
  select(datetime, all_of(variables), parameter)

## grab past met data
met = neon4cast::noaa_stage3() |>
  filter(site_id == "HARV",variable %in% variables,datetime < forecast_date) |>
  collect()
met_daily = met |>
  mutate(datetime = lubridate::as_date(datetime)) |> 
  # mean daily forecasts at each site per ensemble
  group_by(datetime, variable) |> 
  summarize(prediction = mean(prediction)) |>
  pivot_wider(names_from = variable, values_from = prediction) |>
  # convert to Celsius
  mutate(air_temperature = air_temperature - 273.15) |> 
  select(datetime, all_of(variables))

#### Step 3.0: Define the forecasts model for a site
site = "HARV"
target_variable = "gcc_90"
horiz = 35
#forecast_site <- function(site, target_variable, horiz,step) {

message(paste0("Running site: ", site))

# Get site information for elevation
#site_info <- site_data |> dplyr::filter(field_site_id == site)

# Format site data for arima model
site_target_raw <- target |>
  dplyr::select(datetime, site_id, variable, observation) |>
  dplyr::filter(variable == target_variable, 
                site_id == site) |> 
  tidyr::pivot_wider(names_from = "variable", values_from = "observation")

#  if(!target_variable%in%names(site_target_raw)||sum(!is.na(site_target_raw[target_variable]))==0){
#    message(paste0("No target observations at site ",site,". Skipping forecasts at this site."))
#    return()

# } else {

# if(theme %in% c("ticks","beetles")){
#   site_target = site_target_raw %>%
#     filter(wday(datetime,label = T)=="Mon")|>
#     complete(datetime = full_seq(datetime,step),site_id)
#   #Find the most recent Monday
#   mon = Sys.Date()-abs(1-as.numeric(strftime(Sys.Date(), "%u")))
#   h = as.numeric(floor((mon-max(site_target$datetime))/step)+horiz)
# } else {
site_target = site_target_raw |>
  complete(datetime = full_seq(datetime,1),site_id)
h = as.numeric(forecast_date-max(site_target$datetime)+horiz)
#}

site_target_past = site_target |> filter(datetime < forecast_date) |>
  right_join(met_daily,"datetime") ## merge in covariate data

# Fit arima model
if(sum(site_target[target_variable]<0,na.rm=T)>0){#If there are any negative values, don't consider transformation
  fit = forecast::auto.arima(site_target_past[target_variable])
} else {
  fit = forecast::auto.arima(site_target_past[target_variable], 
                             xreg = as.matrix(site_target_past[,c("air_temperature","relative_humidity")]),
                             #                       lambda = "auto",
                             max.d = 0, max.D = 0,max.q=0)
}

lambda = fit$lambda
boxcox = function(x,lambda){(x^lambda-1)/lambda}
inv.boxcox = function(x,lambda){
  ((x * lambda) + 1)^(1 / lambda) - 1
}

## not currently generalized to multiple types of arima model
arima.fx <- function(inputs,drivers,epsilon,horiz=35,lag=1){
  IC = inputs[,"dayof"]
  met.ens = inputs[,"noaa_ensemble_member"]
  sig.ens = inputs[,"sigma.ens"]
  param   = inputs[,1:3]
  betas = param[,which(colnames(param) %in% variables)]
  if(is.null(dim(IC))) IC = as.matrix(IC,ncol=1)
  IC.bc = IC - param[,"intercept"] #boxcox(IC,lambda)
  X = matrix(NA,nrow=nrow(IC),ncol=horiz+lag+1)
  X[,seq_len(lag)] = IC.bc
  dates = sort(unique(drivers$datetime))
  for(t in lag + (0:horiz)){
    met = drivers |>
      filter(datetime == dates[t]) |>
      ungroup() |>
      select(air_temperature,relative_humidity)
    met = met[met.ens,]
    ## gap filling
    met$air_temperature[is.na(met$air_temperature)]     = mean(met$air_temperature,na.rm = TRUE)
    met$relative_humidity[is.na(met$relative_humidity)] = mean(met$relative_humidity,na.rm = TRUE)
    #met[,nrow(met)+1] = apply(met,2,mean,na.rm=TRUE)
    
    X[,t+1] = unlist(X[,t] + met[,1] * betas[,1] + met[,2]*betas[,2] + epsilon[sig.ens,t] )
    

  }
  y = X + param[,"intercept"]#inv.boxcox(X[,-lag],lambda)
  return(y)
}

## build ensembles
source("paramIC.R")
lag = 1
param = paramSamples
IC    = IC[,lag]
met.ens = sample(1:31,ne,replace=TRUE)
drivers = noaa_future_daily
epsilon = rmvnorm(ne,rep(0,horiz+lag),diag(rep(fit$sigma2,horiz+lag)))

source("sobolDesign.R")
y1 = arima.fx(x1,drivers,epsilon,horiz=35)
y2 = arima.fx(x2,drivers,epsilon,horiz=35)

plot(y1[1,],type='l')
for(i in 1:ne){lines(y1[i,])}
for(i in 1:ne){lines(y2[i,],col=2)}

#### Validation uncertainty #### 

# model time series with CI and data
ybar = apply(y1,2,quantile,c(0.025,0.5,0.975),na.rm=TRUE)

apply(is.na(y1),2,sum)

# raw error vs lead time

# quantile error vs lead time

# crps vs lead time

#### Sobol Analysis ####

