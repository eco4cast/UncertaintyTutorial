#Script creates trained model for each target variable using Random Forest in regression mode

#### NOTE: Re-running this script will NOT overwrite previously saved models - if they are not deleted, the downstream forecasts
# using the saved model objects will not run correctly

#### Step 1: Load libraries
library(here)
library(ranger) #needed for random forest implementation
library(doParallel)
library(tidyverse)
library(tidymodels)
library(vip)
library(butcher)
library(bundle)
library(neon4cast)
library(lubridate)
library(rMR)
library(glue)
library(decor)
#source("ignore_sigpipe.R")
library(tsibble)
library(fable)
library(arrow)
library(sensitivity)


here::i_am("randfor.R")
source(here("R/download_target.R"))

# Set model types
model_themes = c("phenology") #This model is only relevant for three themes
model_types = c("phenology")


#### Step 2: Get NOAA driver data

forecast_date <- lubridate::date("2023-05-04")
noaa_date <- lubridate::date("2023-05-04") - lubridate::days(1)  #Need to use yesterday's NOAA forecast because today's is not available yet
site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") %>%
  filter(field_site_id %in% c("HARV"))|> # can be useful for testing
  filter(if_any(matches(model_types),~.==1))
all_sites = site_data$field_site_id


load_stage2 <- function(site, endpoint, variables){
  message('run ', site)
  use_bucket <- paste0("neon4cast-drivers/noaa/gefs-v12/stage3/parquet/", site)
  use_s3 <- arrow::s3_bucket(use_bucket, endpoint_override = endpoint, anonymous = TRUE)
  parquet_file <- arrow::open_dataset(use_s3) |>
    dplyr::collect() |>
    dplyr::filter(datetime >= lubridate::ymd('2017-01-01'),
                  variable %in% variables)|> #It would be more efficient to filter before collecting, but this is not running on my M1 mac
    na.omit() |> 
    mutate(datetime = lubridate::as_date(datetime)) |> 
    group_by(datetime, site_id, variable) |> 
    summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
    pivot_wider(names_from = variable, values_from = prediction) |> 
    # convert air temp to C
    mutate(air_temperature = air_temperature - 273.15)
}


# Specify desired met variables - all meteo
variables <- c('air_temperature',
               "surface_downwelling_longwave_flux_in_air",
               "surface_downwelling_shortwave_flux_in_air",
               "precipitation_flux",
               "air_pressure",
               "relative_humidity",
               "air_temperature",
               "northward_wind",
               "eastward_wind")

#Code from Freya Olsson to download and format meteorological data (had to be modified to deal with arrow issue on M1 mac). Major thanks to Freya here!!

# Load stage 3 data
endpoint = "data.ecoforecast.org"
use_bucket <- paste0("neon4cast-drivers/noaa/gefs-v12/stage2/parquet/0/", noaa_date)


load_stage3 <- function(site, endpoint, variables){
  message('run ', site)
  use_bucket <- paste0("neon4cast-drivers/noaa/gefs-v12/stage3/parquet/", site)
  use_s3 <- arrow::s3_bucket(use_bucket, endpoint_override = endpoint, anonymous = TRUE)
  parquet_file <- arrow::open_dataset(use_s3) |>
    dplyr::collect() |>
    dplyr::filter(datetime >= noaa_date,
                  variable %in% variables)|> #It would be more efficient to filter before collecting, but this is not running on my M1 mac
    na.omit() |> 
    mutate(datetime = lubridate::as_date(datetime)) |> 
    group_by(datetime, site_id, variable, parameter) |> 
    summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
    pivot_wider(names_from = variable, values_from = prediction) |> 
    # convert air temp to C
    mutate(air_temperature = air_temperature - 273.15)
}


if(file.exists(here("Forecast_submissions/Generate_forecasts/noaa_downloads/past_allmeteo.csv"))) {
  noaa_past <- read_csv(here("Forecast_submissions/Generate_forecasts/noaa_downloads/past_allmeteo.csv"))
} else {
  noaa_past <- map_dfr(all_sites, load_stage2, endpoint, variables)
  noaa_future <- map_dfr(all_sites, load_stage3, endpoint, variables)
}
noaa_past$site_id <- 1
noaa_future$site_id <- 1

############################################ SET UP TRAINING LOOPS ###################################

##### Training function ##########

train_site <- function(sites, noaa_past, target_variable) {
  message(paste0("Running ",target_variable," at all sites"))
  
  # browser()
  # Merge in past NOAA data into the targets file, matching by date.
  site_target <- target |>
    dplyr::select(datetime, site_id, variable, observation) |>
    dplyr::filter(variable %in% c(target_variable), 
                  site_id %in% sites) |>
    tidyr::pivot_wider(names_from = "variable", values_from = "observation") |>
    dplyr::left_join(noaa_past%>%
                       filter(site_id %in% sites), by = c("datetime", "site_id"))|>
    drop_na() |> #removes non-complete cases - BEWARE
    select(-datetime, -site_id)
  
  if(!target_variable%in%names(site_target)){
    message(paste0("No target observations at site ",site,". Skipping forecasts at this site."))
    return()
    
  } else if(sum(!is.na(site_target$air_temperature)&!is.na(site_target[target_variable]))==0){
    message(paste0("No historical air temp data that corresponds with target observations at site ",site,". Skipping forecasts at this site."))
    return()
    
  } else {
    # Tune and fit lasso model - making use of tidymodels
    
    
    #Recipe for training models
    rec_base <- recipe(site_target)|>
      #step_rm(c("datetime", "site_id", "parameter"))|>
      update_role(everything(), new_role = "predictor")|>
      update_role({{target_variable}}, new_role = "outcome")|>
      #step_dummy(site_id)|> #Random forest handles categorical predictor without need to convert to dummy 
      step_normalize(all_numeric(), -all_outcomes())
    
    ## Set up tuning and fitting engine
    
    tune_randfor <- rand_forest(
      mtry = tune(),
      trees = 500,
      min_n = tune()) |>
      set_mode("regression") %>%
      set_engine("ranger", importance = "impurity", keep.inbag = TRUE) 
    
    #k-fold cross-validation
    randfor_resamp <- vfold_cv(site_target, v = 10, repeats = 5)# define k-fold cross validation procedure 
    ## Assemble workflow and tune
    wf <- workflow() %>%
      add_recipe(rec_base)
    
    #Tune models
    #If running in parallel  
    library(doParallel)
    cl <- makePSOCKcluster(16) #SET 
    registerDoParallel(cl) 
    randfor_grid <- 
      tune_grid(
        wf %>% add_model(tune_randfor),
        resamples = randfor_resamp,
        grid = 20 
      )
    
    ## Select best model via RMSE
    best_mod<-randfor_grid|>
      select_best("rmse")
    
    #select model with best tuning parameter by RMSE, cross-validation approach
    final_mod <- finalize_workflow(
      wf %>% add_model(tune_randfor),
      best_mod
    )
    
    final_fit <- fit(final_mod, site_target)
    
    vip <- final_fit|>extract_fit_parsnip()|>vip()|>pluck("data")|>
      pivot_wider(names_from = "Variable", values_from = "Importance", names_prefix = "importance_")
    
    
    final_preds <- predict(final_fit, site_target)|>
      bind_cols(site_target)
    
    final_rmse<-rmse(final_preds, estimate = .pred, truth = {{target_variable}})
    #try to extract fit and write to tibble variable importance as columns bind_cols
    
    #save model fit in minimal form
    res_bundle <-
      final_fit %>%            
      bundle()
    
    predictor_formula <- formula(rec_base|>prep())|>as.character()|>pluck(3)
    
    saveRDS(res_bundle, here(paste0("trained_models/rf/", paste(theme, target_variable,"trained",Sys.Date(), sep = "-"), ".Rds")))
    tibble(theme = theme, site = "all sites", 
           predictor_formula = predictor_formula, n_obs = nrow(site_target), 
           n_vfolds = "10", target_variable = target_variable, 
           rmse = final_rmse$.estimate, mtry = best_mod$mtry, min_n = best_mod$min_n)|>
      bind_cols(vip)
    
  }
}


######### Loop to train all sites ########

for (theme in model_themes) {
  target = download_target(theme)
  target$site_id <- 1
  type = ifelse(theme%in% c("terrestrial_30min", "terrestrial_daily"),"terrestrial",theme)
  
  if("siteID" %in% colnames(target)){ #Sometimes the site is called siteID instead of site_id. Fixing here
    target = target%>%
      rename(site_id = siteID)
  }
  if("time" %in% colnames(target)){ #Sometimes the time column is instead labeled "datetime"
    target = target%>%
      rename(datetime = time)
  }
  
  site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv")|>
    filter(field_site_id %in% c("HARV"))|> # can be useful for testing
    filter(get(type)==1) 
  
  #sites = site_data$field_site_id
  sites = 1
  
  #Set target variables for each theme
  if(theme == "aquatics")           {vars = c("temperature","oxygen","chla")}
  if(theme == "phenology")          {vars = c("gcc_90")}
  if(theme == "terrestrial_daily")  {vars = c("nee","le")}
  if(theme == "beetles")            {vars = c("abundance","richness")}
  if(theme == "ticks")              {vars = c("amblyomma_americanum")}
  
  
  mod_summaries <- map(vars, ~train_site(sites = sites, target_variable = ., noaa_past = noaa_past))|> #possibly only accepts static values so can't map '.x' into site or target_variable
    compact()|>
    list_rbind()|>
    drop_na()
  
  assign(x = paste0(theme, "_mod_summaries"), value = mod_summaries)
  
}


mod_sums_all <- syms(apropos("_mod_summaries"))|>
  map_dfr(~eval(.)|>bind_rows())|>
  write_csv(here(paste0("trained_models/rf/model_training_summaries-", Sys.Date(),".csv")))

gcc_model <- readRDS("trained_models/rf/phenology-gcc_90-trained-2023-06-22.Rds") |> unbundle()
#predict(gcc_model$object$fit$fit$object, noaa_future, type = "raw") -> outputs
#predict(gcc_model$object$fit$fit$object, noaa_future_a[,2:11], type = "conf_int") -> outputs_unc
#predict(gcc_model$object$fit$fit$object, noaa_future_b, type = "numeric") -> outputs_numeric

rec <- prep(extract_preprocessor(gcc_model))
N_ens <- 50000
noaa_futures <- noaa_future |>
  filter(datetime >= noaa_date & datetime <= as.Date("2021-06-30")) |> 
  group_by(datetime) |> 
  group_split()

sobol_time <- tibble(
  forecast_date = seq.Date(noaa_date, as.Date("2021-06-30"), by = "day"),
  ensemble = (noaa_futures))

library(future)
library(furrr)
library(future.callr)
plan(callr, workers = 8)


source("R/predict_rf_ensemble.R")
model_test <- function(X, mod = gcc_model$fit$fit$fit) {
  y <- predict_tidy_ensemble(mod, X)
  #  y <- y$.pred
  y
}

sobol_time_fit <- sobol_time |> 
  mutate(sobol = future_map(noaa_futures, function(x) {
    noaa_future2 <- x |> 
      select(-datetime, -site_id, -parameter) |> 
      sample_n(N_ens, replace = TRUE)
    
    noaa_future2 <- bake(rec, noaa_future2) |> 
      mutate(parameter_seed = sample.int(n()))
    
    noaa_future_a <- noaa_future2 |>
      filter(parameter_seed %in% 1:round(N_ens/2))
    noaa_future_b <- noaa_future2 |>
      filter(parameter_seed %in% round(N_ens/2 + 1):N_ens) 
    
    for(i in seq_len(ncol(noaa_future_b))) {
      noaa_future_b[[i]] <- sample(noaa_future_b[[i]], length(noaa_future_b[[i]]), replace = FALSE)
    }
    
    s2 <- soboljansen(model = model_test, X1 = noaa_future_a, X2 = noaa_future_b, nboot = 0, conf = 0.95)
    
    return(s2)
    
  }))


sobol_time_fit_data <- sobol_time_fit |> 
  mutate(T = map(sobol, function(x) {
    x$T |> rownames_to_column(var = "variable")
  })) |> 
  select(forecast_date, T) |> 
  unnest(T) |> 
  mutate(variable = fct_recode(variable, tree_uncertainty = "parameter_seed")) |> 
  mutate(variable = fct_relevel(variable, "tree_uncertainty", after = Inf))

ggplot(sobol_time_fit_data, aes(x = forecast_date, fill = variable, y = original)) +
  geom_area()



