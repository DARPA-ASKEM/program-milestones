# August 23, 2023
# File to read in data for T2 cleanly
# Goal is to export data to a script to compute metrics and make plots

# Updated February 7, 2024
# Remove specific working directory references to allow TA4/others to use them

# This file requires input CSVs from CIEMSS (T2_ensemble_of1) which are also uploaded
################################################################################

# Set WD to folder with ensemble output
setwd()...

# Load libraries; (install if you do not have them)
library(evalcast)
library(covidcast)
library(tibble)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readr)
library(lubridate)
library(dplyr)

################################################################################
# Function to convert_performer_input assuming it is of the form found here
# https://raw.githubusercontent.com/reichlab/covid19-forecast-hub/master/data-processed/COVIDhub-ensemble/2023-06-26-COVIDhub-ensemble.csv

# Function to convert performer input to the format the package expects
# This has two arguments
# 1) Input (CSV form of predictions. Downloaded from Reichlab Github for this example).
# 2) Metric (whether they are forecasting hospitalizations, cases, or deaths)

convert_performer_input <- function(input, metric = c("hosp","case","death")){
  
  # Find which rows are for death, hosps, or cases
  case.ind <- grep(pattern = "inc case", x = input$target)
  # Changed this to incident instead of cumulative
  
  # THIS DOES NOT MATCH WITH THE SIGNALS SINCE THE DATA ARE MISLABELED
  death.ind <- grep(pattern = "cum death", x = input$target)
  hosp.ind <- grep(pattern = "inc hosp", x = input$target)
  
  # Enumerate data sources and signals
  # These are added to the data frame based on which one it is
  data.sources <- c("jhu-csse", "jhu-csse", "hhs")
  
  # Changed this to incidence
  data.signals <- c("confirmed_incidence_num", "deaths_cumulative_num","confirmed_admissions_covid_1d")
  
  # Select the columns we want, compute "ahead", and add incidence_period, data_source, and signal
  
  # If we are interested in cases
  if (metric == "case") {
    
    output <- input %>%
      
      # Filter to only include the rows with case in them
      filter(row_number() %in% case.ind) %>%
      
      # Select the variables we want
      select(forecaster, quantile, value, 
             forecast_date, target_end_date, location) %>%
      
      # Compute ahead and also add incidence_period, data_source,  and signal
      mutate(ahead = as.numeric(target_end_date - forecast_date)) %>%
      mutate(incidence_period = "day") %>% #or epiweek
      mutate(data_source = data.sources[1]) %>%
      mutate(signal = data.signals[1])
    
    # Do the same for death
  } else if (metric == "death"){
    output <- input %>%
      filter(row_number() %in% death.ind) %>%
      select(forecaster, quantile, value, 
             forecast_date, target_end_date, location) %>%
      mutate(ahead = as.numeric(target_end_date - forecast_date)) %>%
      mutate(incidence_period = "day") %>% #or epiweek instead of day
      mutate(data_source = data.sources[2]) %>%
      mutate(signal = data.signals[2])
    
    # And do the same for hospitalizations
  } else if (metric == "hosp") {
    output <- input %>%
      filter(row_number() %in% hosp.ind) %>%
      select(forecaster, quantile, value, 
             forecast_date, target_end_date, location) %>%
      mutate(ahead = as.numeric(target_end_date - forecast_date)) %>%
      mutate(incidence_period = "day") %>% # or epiweek
      mutate(data_source = data.sources[3]) %>%
      mutate(signal = data.signals[3])
  }
  
  # Get a df with state:fips map
  # This is because the Github data is in FIPS but we need them in states
  geo.df <- covidcast::state_census %>% select(STATE, ABBR)
  
  # Merge geo.df onto output so that these have the two character state names
  # Ch
  geo = merge(x = geo.df, y = output, by.x = "ABBR", by.y = "location")
  
  # Convert it to lowercase
  geo$ABBR = tolower(geo$ABBR)
  
  # Rename and select columns
  # Also arrange them in the order that matches the data from the packages
  geo2 <- geo %>% 
    select(-"STATE") %>%
    rename(geo_value = "ABBR") %>%
    group_by(geo_value) %>%
    arrange(desc(geo_value), forecast_date, target_end_date, quantile) %>%
    select(ahead, geo_value, quantile, value, forecaster, forecast_date, data_source, signal, target_end_date, incidence_period)
  
}

################################################################################

# Set API key (get an API key from covicast if you do not have one)
# options(covidcast.auth = "<your-API-key-here>")

# Make signal deaths_cumulative_num
signal.var = "deaths_cumulative_num"

# Define signals
signals <- tibble(data_source = "jhu-csse", 
                  signal = signal.var, 
                  start_day = "2021-05-12",
                  geo_type = "state",
                  geo_value = "ny")

# Define forecast dates
forecast_dates_dec <- seq.Date(from = as.Date("2021-05-13"),
                               to = as.Date("2021-08-17"),
                               by = 1)

################################################################################
# Pull down data from other forecasters

# Grab all the ones we compare against in one call
predictions <- get_covidhub_predictions(
  covidhub_forecaster_name = 
    c("COVIDhub-4_week_ensemble", 
      "CU-select",
      "UChicagoCHATTOPADHYAY-UnIT",
      "LNQ-ens1",
      "UMass-MechBayes",
      "Karlen-pypm",
      "MOBS-GLEAM_COVID"),
  as_date(forecast_dates_dec), 
  ahead = NULL,
  geo_values = "ny",
  signal = "deaths_cumulative_num") %>%
  filter(nchar(geo_value)==2)
      
################################################################################
# Read in T2 data

# Find the files which are single and which are ensemble
T2.single <- list.files(pattern = "T2_ensemble_of1")
T2.ensemble <- list.files(pattern = "T2_ensemble_of2")

# Read in the single forecasts
single.forecasts <- NULL
for(i in T2.single){
  a <- read_csv(i)
  single.forecasts <- rbind(single.forecasts, a)
}

# Read in the ensemble forecasts
ensemble.forecasts <- NULL
for(i in T2.ensemble){
  a <- read_csv(i)
  ensemble.forecasts <- rbind(ensemble.forecasts, a)
}

#### Amend them to the right format

# Drop the first column of both files
single.forecasts$...1 <- ensemble.forecasts$...1 <- NULL

# Fix the location in both files
single.forecasts$location <- ensemble.forecasts$location <- "NY"

# Add forecaster name
single.forecasts$forecaster <- "performer-single"
ensemble.forecasts$forecaster <- "performer-ensemble"

# Pull out cumulative deaths
single.deaths <- convert_performer_input(single.forecasts, metric = "death")
ensemble.deaths <- convert_performer_input(ensemble.forecasts, metric = "death")

# Remove forecasts from the performers that aren't done every week
# This means picking only those whose (days ahead remainder 7) = 0
single.deaths.week <- 
  single.deaths %>%
  filter(ahead %% 7 == 0)

ensemble.deaths.week <- 
  ensemble.deaths %>%
  filter(ahead %% 7 == 0)

# Bind them all together
all.predictions <- bind_rows(predictions, single.deaths.week, ensemble.deaths.week)

# T2: Save all predictions as an RDS object
saveRDS(all.predictions, file = "T2-cum-deaths-predictions.RDS")
