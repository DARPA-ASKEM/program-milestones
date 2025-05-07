# April 10, 2025
# Code to read in ensemble forecast challenge predictions
################################################################################
# Load libraries
library(here)
library(readr)
library(stringr)
library(covidHubUtils)
library(dplyr)
library(zoltr)
library(hubEvals)
library(httr)
library(evalcast)
library(covidcast)
library(ggplot2)

# Read in data
pa1 <- read_csv(here("CIEMSS ensembles/Pennsylvania_single_flu_ensemble_results.csv"))
pa.ens <- read_csv(here("CIEMSS ensembles/Pennsylvania_flu_ensemble_results.csv"))

fl1 <- read_csv(here("CIEMSS ensembles/Florida_single_flu_ensemble_results.csv"))
fl.ens <- read_csv(here("CIEMSS ensembles/Florida_flu_ensemble_results.csv"))

################################################################################
# Fix target column in PA data

pa1$target <- ifelse(str_detect(pa1$target, "1.0 week"),"1 wk inc hosp",
              ifelse(str_detect(pa1$target, "2.0 week"), "2 wk inc hosp",
              ifelse(str_detect(pa1$target, "3.0 week"), "3 wk inc hosp",
                     "4 wk inc hosp")))
head(pa1$target)
table(pa1$target)

pa.ens$target <- ifelse(str_detect(pa.ens$target, "1.0 week"),"1 wk inc hosp",
                 ifelse(str_detect(pa.ens$target, "2.0 week"), "2 wk inc hosp",
                 ifelse(str_detect(pa.ens$target, "3.0 week"), "3 wk inc hosp",
                                   "4 wk inc hosp")))
table(pa.ens$target)

# Fix additional columns in PA data
pa.ens$unit <- pa1$unit <- "42"
pa.ens$model <- "PAEns"
pa1$model <- "PA1"
pa1$timezero <- pa1$forecast_date
pa.ens$timezero <- pa.ens$forecast_date

# Reorder columns to put model first
pa1 <- pa1 %>% relocate(model)
pa.ens <- pa.ens %>% relocate(model)
################################################################################
# Fix target column in FL data
fl1$target <- ifelse(str_detect(fl1$target, "1.0 week"),"1 wk inc hosp",
              ifelse(str_detect(fl1$target, "2.0 week"), "2 wk inc hosp",
              ifelse(str_detect(fl1$target, "3.0 week"), "3 wk inc hosp",
                                   "4 wk inc hosp")))
head(fl1$target)
table(fl1$target)

fl.ens$target <- ifelse(str_detect(fl.ens$target, "1.0 week"),"1 wk inc hosp",
                 ifelse(str_detect(fl.ens$target, "2.0 week"), "2 wk inc hosp",
                 ifelse(str_detect(fl.ens$target, "3.0 week"), "3 wk inc hosp",
                                      "4 wk inc hosp")))
table(fl.ens$target)

# Fix additional columns in florida data
fl1$unit <- fl.ens$unit <- "12"

fl1$model <- "FL1"
fl.ens$model <- "FLEns"

fl1$timezero <- fl1$forecast_date
fl.ens$timezero <- fl.ens$forecast_date

# Reorder columns to put model first
fl1 <- fl1 %>% relocate(model)
fl.ens <- fl.ens %>% relocate(model)
################################################################################
# Connect to Zoltar

# Set this to FALSE before running; set it back to TRUE after done pulling from script
httr::set_config(config(ssl_verifypeer = FALSE))

# Make a new zoltar connection
zoltar_connection <- new_connection()

# Authenticate (this pulls from your .Renviron file)
zoltar_authenticate(zoltar_connection, 
                    Sys.getenv("Z_USERNAME"), 
                    Sys.getenv("Z_PASSWORD"))

# Check the connection
zoltar_connection

################################################################################
# Get the right data and pointers for Zoltar

# Load the zoltar projects
the_projects <- projects(zoltar_connection)

# Get the FluSight project URL
project_url <- 
  the_projects[the_projects$name == "CDC FluSight Forecast Hub", "url"]

# Use the zoltar connection and the project URL to pull down the correct models
the_models <- models(zoltar_connection, project_url)

################################################################################
# Define times of forecasts
# These are the target dates for the forecasts minus 7
times.seq <- 
  seq.Date(from = as.Date("2023-12-09"), to = as.Date("2024-01-06"), by = 7)

# Get FL and PA forecasts
fl.pa.data <- do_zoltar_query(zoltar_connection, 
                                  project_url, 
                                  query_type = "forecasts", 
                                  # this is where we pull different models
                                  # Can pass the_models$model_abbr and get all the models
                                  models = the_models$model_abbr,
                                  # Units refer to different locations
                                  units = c("12","42"), 
                                  # Different targets go here
                                  targets = c("0 wk inc flu hosp", 
                                              "1 wk inc flu hosp", 
                                              "2 wk inc flu hosp",
                                              "3 wk inc flu hosp"),
                                  timezeros = times.seq,
                                  types = c("point", "quantile"))

# Convert these to 1, 2, 3, 4 wk forecasts
fl.pa.data$target.old <- fl.pa.data$target
fl.pa.data$target <- 
  ifelse(fl.pa.data$target.old == "3 wk inc flu hosp", "4 wk inc flu hosp",
  ifelse(fl.pa.data$target.old == "2 wk inc flu hosp", "3 wk inc flu hosp",
  ifelse(fl.pa.data$target.old == "1 wk inc flu hosp", "2 wk inc flu hosp",
         "1 wk inc flu hosp")))

################################################################################
# Filter FL forecasts to be only from 2023-12-09
flor.data <- fl.pa.data %>%
  filter(unit == "12") %>%
  select(-(family:param3)) %>%
  # Group by model, timezero, and target (1, 2, or 3 week inc hosp)
  group_by(model, timezero, target) %>%
  mutate(current.group.id = cur_group_id()) %>%
  filter(timezero == "2023-12-09")

head(flor.data)

# Setup FL data to be in the right format (same as other forecasters)
# Make changes here 

fl1.card <- fl1 %>%
  mutate(cat = NA, prob = NA, sample = NA,
         season = "2023-2024") %>%
  rename(class = type) %>%
  group_by(model, timezero, target) %>%
  mutate(current.group.id = cur_group_id() + 1000) %>%
  select(model, timezero, season, unit, target, class, 
         value, cat, prob, sample, quantile, current.group.id) # %>%
  # filter(target != "1 wk inc hosp")

fl.ens.card <- fl.ens %>%
  mutate(cat = NA, prob = NA, sample = NA,
         season = "2023-2024") %>%
  rename(class = type) %>%
  group_by(model, timezero, target) %>%
  mutate(current.group.id = cur_group_id() + 2000) %>%
  select(model, timezero, season, unit, target, class, 
         value, cat, prob, sample, quantile, current.group.id) # %>%
  # filter(target != "1 wk inc hosp")

################################################################################
# Filter PA forecasts to be only from 2023-12-09
pa.data <- fl.pa.data %>%
  filter(unit == "42") %>%
  select(-(family:param3)) %>%
  # Group by model, timezero, and target (1, 2, or 3 week inc hosp)
  group_by(model, timezero, target) %>%
  mutate(current.group.id = cur_group_id()) %>%
  filter(timezero == "2023-12-09") 

# Setup PA data to be in the right format (same as other forecasters)
pa1.card <- pa1 %>%
  mutate(cat = NA, prob = NA, sample = NA,
         season = "2023-2024") %>%
  rename(class = type) %>%
  group_by(model, timezero, target) %>%
  mutate(current.group.id = cur_group_id() + 1000) %>%
  select(model, timezero, season, unit, target, class, 
         value, cat, prob, sample, quantile, current.group.id) # %>%
  # filter(target != "1 wk inc hosp")

pa.ens.card <- pa.ens %>%
  mutate(cat = NA, prob = NA, sample = NA,
         season = "2023-2024") %>%
  rename(class = type) %>%
  group_by(model, timezero, target) %>%
  mutate(current.group.id = cur_group_id() + 2000) %>%
  select(model, timezero, season, unit, target, class, 
         value, cat, prob, sample, quantile, current.group.id) # %>%
  # filter(target != "1 wk inc hosp")

################################################################################
# Bind the PA and FL forecasts from CIEMSS together to put into loop
pa.all <- rbind(pa.data, pa1.card, pa.ens.card)
fl.all <- rbind(flor.data, fl1.card, fl.ens.card)

################################################################################
# Get truth data for PA and FL forecasts
# These are the target dates for the forecasts minus 7
times.seq.truth <- 
  seq.Date(from = as.Date("2023-12-09"), to = as.Date("2023-12-30"), by = 7)

# Get truth data for FL and PA
truth.fl.pa <- do_zoltar_query(zoltar_connection, 
                              project_url, 
                              query_type = "truth", 
                              models = NULL, 
                              # Units refer to different locations
                              units = c("12","42"), 
                              # Different targets go here
                              # only 1 wk inc flu hosp works for target
                              targets = c("0 wk inc flu hosp",
                                          "1 wk inc flu hosp", 
                                          "2 wk inc flu hosp", 
                                          "3 wk inc flu hosp"),
                              timezeros = times.seq.truth, 
                              as_of = "2024-12-03 12:00:00 UTC")


# Filter data for just FL
fl.truth <- truth.fl.pa %>% 
  filter(unit == "12") %>%
  rename(actual = value) %>%
  mutate(target_end_date = timezero + 7) %>%
  select(-c(timezero,target)) %>%
  group_by(unit, target_end_date) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(-unit)

# Filter data for just PA
pa.truth <- truth.fl.pa %>% 
  filter(unit == "42") %>%
  rename(actual = value) %>%
  mutate(target_end_date = timezero + 7) %>%
  select(-c(timezero,target)) %>%
  group_by(unit, target_end_date) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(-unit)

head(fl.truth)

head(pa.truth)

################################################################################
# Florida loop cards
fl.loop.cards <- NULL

# Reminder: check 0 wk inc hosp

for(i in unique(fl.all$current.group.id)){
  
  # get current group
  current.g <- fl.all %>% 
    filter(current.group.id == i) %>%
    select(-current.group.id) %>% 
    rename(forecaster = model)
  
  # evaluate the current group
  current.eval <- evaluate_predictions(
    predictions_cards = current.g,
    truth_data = fl.truth,
    err_measures = list(wis = weighted_interval_score, 
                        ae = absolute_error,
                        coverage_80 = interval_coverage(coverage = 0.8),
                        spread = sharpness))
  
  # Add code to get the target and filter to only include the forecasts for that target
  current.target <- as.numeric(substr(current.eval$target[1],1,1))
  
  # Filter the card to only be the relevant time window
  current.eval2 <- current.eval %>%
    filter(target_end_date == timezero + 7*current.target)
  
  cat(i, "\n")
  # Bind them to each other to make a full scorecard
  fl.loop.cards <- bind_rows(fl.loop.cards, current.eval2)
  cat(nrow(fl.loop.cards), "\n")
  
}
fl.loop.cards

# Add a variable to make the CIEMSS lines dashed
fl.loop.cards <- 
  fl.loop.cards %>% 
  # mutate(Performer = ifelse(forecaster %in% c("FL1","FLEns"),"DARPA Performer","FluSight Forecaster"))
  mutate(Performer = ifelse(forecaster %in% c("FL1"),"DARPA Performer Single Model",
                            ifelse(forecaster %in% c("FLEns"), "DARPA Performer Ensemble",
                                   "FluSight Forecaster")))
################################################################################
# Pennsylvania loop cards
pa.loop.cards <- NULL

# Reminder: check 0 wk inc hosp

for(i in unique(pa.all$current.group.id)){
  
  # get current group
  current.g <- pa.all %>% 
    filter(current.group.id == i) %>%
    select(-current.group.id) %>% 
    rename(forecaster = model)
  
  # evaluate the current group
  current.eval <- evaluate_predictions(
    predictions_cards = current.g,
    truth_data = pa.truth,
    err_measures = list(wis = weighted_interval_score, 
                        ae = absolute_error,
                        coverage_80 = interval_coverage(coverage = 0.8),
                        spread = sharpness))
  
  # Add code to get the target and filter to only include the forecasts for that target
  current.target <- as.numeric(substr(current.eval$target[1],1,1))
  
  # Filter the card to only be the relevant time window
  current.eval2 <- current.eval %>%
    filter(target_end_date == timezero + 7*current.target)
  
  cat(i, "\n")
  # Bind them to each other to make a full scorecard
  pa.loop.cards <- bind_rows(pa.loop.cards, current.eval2)
  cat(nrow(pa.loop.cards), "\n")
  
}
pa.loop.cards

# Add a variable to make the CIEMSS lines dashed
pa.loop.cards <- 
  pa.loop.cards %>% 
#  mutate(Performer = ifelse(forecaster %in% c("PA1","PAEns"),"DARPA Performer","FluSight Forecaster"))
  mutate(Performer = ifelse(forecaster %in% c("PA1"),"DARPA Performer Single Model",
                            ifelse(forecaster %in% c("PAEns"), "DARPA Performer Ensemble",
                                   "FluSight Forecaster")))



################################################################################
# Make graphs: WIS loop cards

fl.loop.cards %>% 
  # Is this the right thing to do? Depends on if we are treating their targets correctly
  # filter(target_end_date < "2024-01-06") %>%
  ggplot(aes(x = target_end_date, y = wis, 
             color = forecaster, group = forecaster,
             linetype = Performer, linewidth = Performer)) +
  scale_linetype_manual(values=c("twodash", "solid", "dotted")) +
  scale_linewidth_manual(values = c(1, 1, 0.5)) +
  geom_line() + 
  geom_point() +
  ggtitle("FL Dec 2023 Forecast: WIS Comparison of inc hosp") + 
  scale_x_date(breaks = fl.truth$target_end_date) + #date_labels =  "%b %Y") +
  xlab("Target date for forecasts") + 
  ylab("Weighted Interval Score") +
  labs(caption = "Forecast and ground truth data obtained from zoltar package") +
  theme_test()


pa.loop.cards %>% 
  # Is this the right thing to do? Depends on if we are treating their targets correctly
  # filter(target_end_date < "2024-01-06") %>%
  ggplot(aes(x = target_end_date, y = wis, 
             color = forecaster, group = forecaster,
             linetype = Performer, linewidth = Performer)) +
  scale_linetype_manual(values=c("twodash", "solid", "dotted")) +
  scale_linewidth_manual(values = c(1, 1, 0.5)) +
  geom_line() + 
  geom_point() +
  ggtitle("PA Dec 2023 Forecast: WIS Comparison of inc hosp") + 
  scale_x_date(breaks = pa.truth$target_end_date) + #date_labels =  "%b %Y") +
  xlab("Target date for forecasts") + 
  ylab("Weighted Interval Score") +
  labs(caption = "Forecast and ground truth data obtained from zoltar package") +
  theme_test()


################################################################################
# Make graphs: median predictions vs real data

# Plot median predictions vs real predictions to see what happened
# Start with PA

# Get median predictions
pa.medians <- pa.all %>% 
  filter(quantile == 0.5) 

# Add variable for x-axis
pa.medians$time.end = rep(times.seq[2:5],35)

# Add a variable to aid the plotting
pa.medians <- pa.medians %>%
  mutate(Performer = ifelse(model %in% c("PA1"),"DARPA Performer Single Model",
                            ifelse(model %in% c("PAEns"), "DARPA Performer Ensemble",
                                   "FluSight Forecaster")))

# Add more variables to aid the plotting
pa.truth$model <- "Actual"
pa.truth$Performer <- "FluSight Actual"

# Note that the plot uses two data frames
# pa.medians at the top and pa.truth later on
pa.medians %>%
  ggplot(aes(x = time.end, y = value, color = model, group = model,
             linetype = Performer, linewidth = Performer)) + 
  geom_line() + 
  geom_point() +
  ggtitle("Pennsylvania Forecasts vs Real Data") + 
  scale_x_date(breaks = times.seq[2:5]) + #date_labels =  "%b %Y") +
  xlab("Target date for forecasts") + 
  ylab("1 week incident hospitalizations") +
  labs(caption = "Forecast and ground truth data obtained from zoltar package") +
  #labs(caption = "Comparing forecasts") +
  theme_test() +
  geom_line(data=pa.truth, aes(x=target_end_date, y = actual), 
            color = "purple") +
  scale_linetype_manual(values=c("twodash","solid", "solid","dotted")) + 
  scale_linewidth_manual(values = c(1, 1, 1, 0.5))

################################################################################
# Do the same for Florida

# Get median predictions
fl.medians <- fl.all %>% 
  filter(quantile == 0.5) 

# Make an x-axis variable for the plotting
fl.medians$time.end = rep(times.seq[2:5],35)

# Add variables for plotting
fl.medians <- fl.medians %>%
  mutate(Performer = ifelse(model %in% c("FL1"),"DARPA Performer Single Model",
                            ifelse(model %in% c("FLEns"), "DARPA Performer Ensemble",
                                   "FluSight Forecaster")))

# Add variables for plotting
fl.truth$model <- "Actual"
fl.truth$Performer <- "FluSight Actual"

# Make the plot

# Note that the plot uses two data frames
# pa.medians at the top and pa.truth later on
fl.medians %>%
  ggplot(aes(x = time.end, y = value, color = model, group = model,
             linetype = Performer, linewidth = Performer)) + 
  geom_line() + 
  geom_point() +
  ggtitle("Florida Forecasts vs Real Data") + 
  scale_x_date(breaks = times.seq[2:5]) + #date_labels =  "%b %Y") +
  xlab("Target date for forecasts") + 
  ylab("1 week incident hospitalizations") +
  labs(caption = "Forecast and ground truth data obtained from zoltar package") +
  #labs(caption = "Comparing forecasts") +
  theme_test() +
  geom_line(data=fl.truth, aes(x=target_end_date, y = actual), 
            color = "purple") +
  scale_linetype_manual(values=c("twodash","solid", "solid","dotted")) + 
  scale_linewidth_manual(values = c(1, 1, 1, 0.5))

################################################################################
# Write results to files
write_csv(x = fl.loop.cards, file = "Graph data files/fl-results.csv")
write_csv(x = pa.loop.cards, file = "Graph data files/pa-results.csv")
write_csv(x = fl.medians, file = "Graph data files/fl-predictions.csv")
write_csv(x = pa.medians, file = "Graph data files/pa-predictions.csv")

