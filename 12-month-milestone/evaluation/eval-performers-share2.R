# July 10, 2023
# Eval notebook script
# Augmented from https://cmu-delphi.github.io/covidcast/evalcastR/articles/evalcast.html

################################################################################
# Packages to install about covidcast and evalcast
# The "pak" package below may be necessary to install evalcast in the ways specified below

install.packages("pak")
pak::pkg_install("cmu-delphi/covidcast/R-packages/evalcast")
install.packages("covidcast")

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
# This function currently selects incident cases and hospitalizations and cumulative deaths

convert_performer_input <- function(input, metric = c("hosp","case","death")){
  
  # Find which rows are for death, hosps, or cases
  case.ind <- grep(pattern = "inc case", x = input$target)
  death.ind <- grep(pattern = "cum death", x = input$target)
  hosp.ind <- grep(pattern = "inc hosp", x = input$target)
 
  # Enumerate data sources and signals
  # These are added to the data frame based on which one it is
  data.sources <- c("jhu-csse", "jhu-csse", "hhs")
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
  geo = merge(x = geo.df, y = output, by.x = "STATE", by.y = "location")
  
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

# Access API; define dates used for this forecast; define signals

# Use API key to pull data from the Github
# Register for an API key here: https://cmu-delphi.github.io/delphi-epidata/
# (Click the "request a registered API key" link)
options(covidcast.auth = "put-real-API-key-here")

# Make signal confirmed_admissions_covid_1d
signal.var = "confirmed_admissions_covid_1d"

# Note that the possible options for signal here are
# signal = c("confirmed_incidence_num", "deaths_incidence_num", 
# "deaths_cumulative_num", "confirmed_admissions_covid_1d")

# Define signals to be used later on.
signals <- tibble(data_source = "hhs", 
                  signal = signal.var, 
                  start_day = "2020-06-15",
                  geo_type = "state")

# Get covidhub forecaster names (for reference)
# This is a list of 127 submissions (some from the same universities)
covidhub.names <- get_covidhub_forecaster_names()

# Specify forecast dates
# Get forecast dates
forecast_dates <- get_covidhub_forecast_dates("CMU-TimeSeries")

# Set this to only be 12/07 for the sake of this demo and plotting
forecast_dates_dec <- forecast_dates[forecast_dates >= "2020-12-01" & 
                                       forecast_dates <= "2020-12-08"]

################################################################################
# Get predictions

# Get ensemble predictions for hospitalizations
# This is for this forecasted date in December 2020
# This pulls down multiple files programmatically. 
predictions_cards_ens <- get_covidhub_predictions(
  "COVIDhub-ensemble", as_date(forecast_dates_dec), ahead = 3, 
  signal = "confirmed_admissions_covid_1d") %>%
  filter(nchar(geo_value)==2) # remove counties

# Check that these populate and have numeric values that aren't 0
head(predictions_cards_ens)
summary(predictions_cards_ens$value)

# Get the "baseliner" for hospitalizations
# This is where the signals object gets used
predictions_cards <- get_predictions(baseline_forecaster,
                                     name_of_forecaster = "baseline",
                                     signals = signals,
                                     forecast_dates = forecast_dates_dec,
                                     incidence_period = "epiweek",
                                     forecaster_args = list(
                                       ahead = 3
                                     )
)

# Check output
head(predictions_cards)

# Bind hospital predictions
predictions <- bind_rows(predictions_cards,
                         predictions_cards_ens)

# Note that these only have dates from December in them (can ignore June part)

###############################################################################
# Try an old dataset with hosps, cases, and deaths

# Using GT deep learning submission downloaded from Github
ex2 <- read_csv("gt-12-07-20.csv")

# Call the forecaster "performer" (normally this is pulled from Github)
ex2$forecaster <- "performer"

# Get inputs for hospitalizations, cases, and deaths
# We only use hospitalizations below, but wanted to check that this functionality works
# This seems to work and separate out the datasets into the different possible metrics
out2.hosp <- convert_performer_input(ex2, metric = "hosp")
out2.case <- convert_performer_input(ex2, metric = "case")
out2.death <- convert_performer_input(ex2, metric = "death")

# Limit the forecast dates to those in the ensemble to which we are comparing
predictions.target.dates <- unique(predictions$target_end_date)
predictions.forecast.dates <- unique(predictions$forecast_date)

# Filter out hospitalizations just to the dates we are comparing to (from the package)
out2.hosp.f <- out2.hosp %>%
  filter(forecast_date %in% predictions.forecast.dates) %>%
  filter(target_end_date %in% predictions.target.dates)

# Bind all predictions for December
dec.preds <- bind_rows(predictions, out2.hosp.f)

###############################################################################
# Define error measures
err_measures <- list(wis = weighted_interval_score,
                     ae = absolute_error,
                     spread = sharpness,
                     coverage_80 = interval_coverage(coverage = 0.8))

# Create the scorecards
# Use December predictions via dec.preds
# These compute WIS, absolute error, coverage, and spread for all the forecasters
scorecards <- evaluate_covid_predictions(
  dec.preds,
  err_measures = err_measures,
  backfill_buffer = 10,
  geo_type = "state"
)

# Look at mean WIS, AE, spread
scorecards %>%
  group_by(forecaster) %>%
  summarise(mean.wis = mean(wis, na.rm=T),
            mean.ae = mean(ae,na.rm=T),
            mean.spread = mean(spread, na.rm=T))

# How many rows do these have? Seem comparable
table(scorecards$forecaster)

# Note that this is because we only pulled down data for one forecast date
# We can pull down additional data and get rows for the "performer" for other dates
table(scorecards$forecaster, scorecards$forecast_date)

# The performer has end dates for all the end date in December
# The covid-hub ensemble is missing Dec 26 but otherwise is the same
# The baseline is sporadic
table(scorecards$forecaster, scorecards$target_end_date)

###############################################################################
# Add some example plots. Additional plots may be added later.
scorecards %>% filter(geo_value != "us") %>%
  ggplot(aes(y=forecaster, x=wis, fill=forecaster)) +
  geom_boxplot() +
  facet_wrap(~forecast_date) +
  scale_x_log10() +
  theme_bw() +
  scale_fill_viridis_d() +
  theme(legend.position = "bottom")

# Plot predictions
# Note that "3" and "19" indicate the time between forecast date (12/07) and the target date
# Target dates are 12/10 and 12/26
plot_width(predictions_cards = dec.preds) + scale_y_log10() 

