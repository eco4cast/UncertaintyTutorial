# tg_arima model
# written by ASL, 21 Jan 2023



#### Step 0: load packages

library(tidyverse)
library(neon4cast)
library(lubridate)
library(rMR)
library(glue)
source("ignore_sigpipe.R")
library(tsibble)
library(fable)
library(arrow)
source("download_target.R")
library(forecast)


#### Step 1: Define team name, team members, and theme

team_name <- "EFI Theory"

team_list <- list(list(individualName = list(givenName = "Abby", 
                                             surName = "Lewis"),
                       organizationName = "Virginia Tech",
                       electronicMailAddress = "aslewis@vt.edu")
)

model_id = "tg_arima"
model_themes = c("terrestrial_daily","aquatics","phenology","beetles","ticks") #This model is only relevant for three themes. I am registered for all three
model_types = c("terrestrial","aquatics","phenology","beetles","ticks") #Replace terrestrial daily and 30min with terrestrial
#Options: aquatics, beetles, phenology, terrestrial_30min, terrestrial_daily, ticks


#### Step 2: Get NOAA driver data
forecast_date <- Sys.Date()
noaa_date <- Sys.Date() - lubridate::days(1)  #Need to use yesterday's NOAA forecast because today's is not available yet

#We're going to get data for all sites relevant to this model, so as to not have to re-load data for the same sites
site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") %>%
  filter(if_any(matches(model_types),~.==1))
all_sites = site_data$field_site_id

# specify meteorological variables needed to make predictions
variables <- c('air_temperature',
"surface_downwelling_longwave_flux_in_air",
"surface_downwelling_shortwave_flux_in_air",
"precipitation_flux",
"air_pressure",
"relative_humidity",
"air_temperature",
"northward_wind",
"eastward_wind")

# Load stage 2 data
endpoint = "data.ecoforecast.org"
use_bucket <- paste0("neon4cast-drivers/noaa/gefs-v12/stage2/parquet/0/", noaa_date)
use_s3 <- arrow::s3_bucket(use_bucket, endpoint_override = endpoint, anonymous = TRUE)
noaa_future <- arrow::open_dataset(use_s3) |>
  dplyr::collect() |>
  dplyr::filter(site_id %in% all_sites,
                datetime >= forecast_date,
                variable == variables) 

# Format met forecasts
noaa_future_daily <- noaa_future |> 
  mutate(datetime = lubridate::as_date(datetime)) |> 
  # mean daily forecasts at each site per ensemble
  group_by(datetime, site_id, parameter, variable) |> 
  summarize(prediction = mean(prediction)) |>
  pivot_wider(names_from = variable, values_from = prediction) |>
  # convert to Celsius
  mutate(air_temperature = air_temperature - 273.15) |> 
  select(datetime, site_id, all_of(variables), parameter)

#### Step 3.0: Define the forecasts model for a site
forecast_site <- function(site,noaa_future_daily,target_variable) {
  message(paste0("Running site: ", site))
  
  
  # Get site information for elevation
  site_info <- site_data |> dplyr::filter(field_site_id == site)
  
  mod_file <- list.files(here("Forecast_submissions/Generate_forecasts/tg_randfor/trained_models/"), pattern = paste(theme, site, target_variable, sep = "-"))
  
  if(!file.exists(here(paste0("Forecast_submissions/Generate_forecasts/tg_randfor/trained_models/",mod_file)))){
    message(paste0("No trained model for site ",site,". Skipping forecasts at this site."))
    return()
    
} else {

    #  Get 30-day predicted temperature ensemble at the site
    noaa_future <- noaa_future_daily%>%
      filter(site_id == site)

#### Step 3.0: Define the forecast model for a site
forecast_site <- function(site, target_variable, horiz,step) {
  
  message(paste0("Running site: ", site))
  
  # Get site information for elevation
  #site_info <- site_data |> dplyr::filter(field_site_id == site)
  
  # Format site data for arima model
  site_target_raw <- target |>
    dplyr::select(datetime, site_id, variable, observation) |>
    dplyr::filter(variable == target_variable, 
                  site_id == site) |> 
    tidyr::pivot_wider(names_from = "variable", values_from = "observation")
  
  if(!target_variable%in%names(site_target_raw)||sum(!is.na(site_target_raw[target_variable]))==0){
    message(paste0("No target observations at site ",site,". Skipping forecasts at this site."))
    return()
    
  } else {
    
    if(theme %in% c("ticks","beetles")){
      site_target = site_target_raw %>%
        filter(wday(datetime,label = T)=="Mon")|>
        complete(datetime = full_seq(datetime,step),site_id)
      #Find the most recent Monday
      mon = Sys.Date()-abs(1-as.numeric(strftime(Sys.Date(), "%u")))
      h = as.numeric(floor((mon-max(site_target$datetime))/step)+horiz)
    } else {
      site_target = site_target_raw |>
        complete(datetime = full_seq(datetime,1),site_id)
      h = as.numeric(Sys.Date()-max(site_target$datetime)+horiz)
    }
    
    # Fit arima model
    if(sum(site_target[target_variable]<0,na.rm=T)>0){#If there are any negative values, don't consider transformation
      fit = auto.arima(site_target[target_variable])
    } else {
      fit = auto.arima(site_target[target_variable], lambda = "auto")
    }
    
    # use the model to forecast target variable
    forecast_raw <- as.data.frame(forecast(fit,h=h,level=0.68))%>% #One SD
      mutate(sigma = `Hi 68`-`Point Forecast`)
    
    forecast = data.frame(datetime = (1:h)*step+max(site_target$datetime),
                          reference_datetime = Sys.Date(),
                          site_id = site,
                          family = "normal",
                          variable = target_variable,
                          mu = as.numeric(forecast_raw$`Point Forecast`),
                          sigma = as.numeric(forecast_raw$sigma),
                          model_id = model_id)%>%
      pivot_longer(cols = c(mu,sigma), names_to = "parameter",values_to = "prediction")%>%
      select(model_id, datetime, reference_datetime,
             site_id, family, parameter, variable, prediction)
    return(forecast)
  }
}

#Quick function to repeat for all variables
run_all_vars = function(var,sites,forecast_site,horiz,step){
  
  message(paste0("Running variable: ", var))
  forecast <- map_dfr(sites,forecast_site,var,horiz,step)
  
}

### AND HERE WE GO! We're ready to start forecasting ### 
for (theme in model_themes) {
  if(!theme%in%c("beetles","ticks") | wday(Sys.Date(), label=TRUE)=="Sun"){ #beetles and ticks only want forecasts every Sunday
    #Step 1: Download latest target data and site description data
    target = download_target(theme)
    type = ifelse(theme%in% c("terrestrial_30min", "terrestrial_daily"),"terrestrial",theme)
    
    if("siteID" %in% colnames(target)){ #Sometimes the site is called siteID instead of site_id. Fixing here
      target = target%>%
        rename(site_id = siteID)
    }
    if("time" %in% colnames(target)){ #Sometimes the datetime column is instead labeled "time"
      target = target%>%
        rename(datetime = time)
    }
    
    site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") %>%
      filter(get(type)==1)
    sites = site_data$field_site_id
    
    #Set target variables
    if(type == "aquatics")           {vars = c("temperature","oxygen","chla")
                                                horiz = 30
                                                step = 1
                                                }
    if(type == "ticks")              {vars = c("amblyomma_americanum")
                                                horiz = 52 #52 weeks
                                                step = 7
                                                }
    if(type == "phenology")          {vars = c("gcc_90","rcc_90")
                                                horiz = 30 
                                                step = 1}
    if(type == "beetles")            {vars = c("abundance","richness")
                                                horiz = 52 #52 weeks
                                                step = 7}
    if(theme == "terrestrial_daily")  {vars = c("nee","le")
                                                horiz = 30
                                                step = 1}
    if(theme == "terrestrial_30min")  {vars = c("nee","le")
                                                horiz = 30
                                                step = 1/24/2}
  
    ## Test with a single site first!
    #forecast <- map_dfr(vars,run_all_vars,sites[1],forecast_site,horiz,step)
    
    #Visualize the ensemble predictions -- what do you think?
    #forecast %>%
    #  pivot_wider(names_from = parameter,values_from = prediction)%>%
    #  ggplot(aes(x = datetime)) +
    #  geom_ribbon(aes(ymax = mu + sigma, ymin = mu-sigma))+
    #  geom_line(aes(y = mu),alpha=0.3) +
    #  facet_wrap(~variable, scales = "free")
    
    # Run all sites -- may be slow!
    forecast <- map_dfr(vars,run_all_vars,sites,forecast_site,horiz,step)
    
    #Forecast output file name in standards requires for Challenge.
    # csv.gz means that it will be compressed
    file_date <- Sys.Date() #forecast$reference_datetime[1]
    model_id = "tg_arima"
    forecast_file <- paste0(theme,"-",file_date,"-",model_id,".csv.gz")
    
    forecast <- forecast%>%
      filter(datetime>=file_date)
    
    #Write csv to disk
    write_csv(forecast, forecast_file)
    
    # Step 5: Submit forecast!
    neon4cast::submit(forecast_file = forecast_file, metadata = NULL, ask = FALSE)
  }
}
