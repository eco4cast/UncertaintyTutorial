#Neon4cast has functions for a lot of what we are hoping to do! 
#Here I am reviving a deprecated function from that package to download data for a given target 
#(why was this function deprecated?? unsure)
#https://github.com/eco4cast/neon4cast/blob/main/R/download_target.R

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