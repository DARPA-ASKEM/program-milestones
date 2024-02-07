# August 23, 2023
# File to compute plots and metrics for T2

# Updated February 7th for TA4 use
################################################################################
# Set WD to folder with the RDS file
setwd()...

# Load libraries
library(evalcast)
library(covidcast)
library(tibble)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readr)
library(lubridate)
library(dplyr)

###############################################################################
# Load the RDS file
preds <- readRDS("T2-cum-deaths-predictions.RDS")

###############################################################################
# Read in our own truth data
ny.data <- read_csv("New_York_State_Statewide_COVID-19_Fatalities_by_Age_Group.csv")

# Format date to filter on it
ny.data$`Report Date` <- as.Date(ny.data$`Report Date`, format = "%m/%d/%Y")

# Filter on it
ny2.actual <- ny.data %>%
  filter(`Age Group` == "Statewide Total") %>%
  filter(`Report Date` >= "2021-05-22" & `Report Date` <= "2021-09-18" ) %>%
  select(`Fatality Count`, `Report Date`) %>%
  rename(actual = `Fatality Count`) %>%
  rename(target_end_date = `Report Date`)

# Get actual data
head(ny2.actual)

###############################################################################
# Write a loop to do errors for each (forecaster, forecast_date) pair

# Group by forecaster and forecast date
# Add group id as a variable
preds.g <- preds %>%
  group_by(forecaster, forecast_date) %>%
  mutate(group.ind = cur_group_id())

# Define loop.cards
loop.cards <- NULL

# Get each group individually and run evaluate predictions
for(i in unique(preds.g$group.ind)){

    # get current group
  current.g <- preds.g %>% 
    filter(group.ind == i) %>%
    select(-group.ind)
  
  # evaluate the current group
  current.eval <- evaluate_predictions(
    predictions_cards = current.g,
    truth_data = ny2.actual,
    err_measures = list(wis = weighted_interval_score, 
                        ae = absolute_error,
                        coverage_80 = interval_coverage(coverage = 0.8),
                        spread = sharpness))
  
  cat(i, "\n")
  # Bind them to each other to make a full scorecard
  loop.cards <- bind_rows(loop.cards, current.eval)
  cat(nrow(loop.cards), "\n")
}

# Check output: looks fine
loop.cards %>% 
  group_by(forecaster) %>%
  summarise(mean.wis = mean(wis,na.rm=T),
            num.na = sum(is.na(wis)))

# Print a few rows
head(as.data.frame(loop.cards))

# Check ahead for these
# Performer ahead is on daily basis; ensemble is on weekly basis
# But since we are looking at cumulative deaths, it should be fine either way
loop.cards %>% 
  filter(forecaster == "performer-single") %>% 
  select(ahead) %>% 
  pull() %>% 
  table()

loop.cards %>% 
  filter(forecaster == "COVIDhub-4_week_ensemble") %>% 
  select(ahead) %>% 
  pull() %>% 
  table()

# write.csv(x = loop.cards, file = "../intermediate results/DARPA Plots/t2-cum-death.csv")

###############################################################################
# This generates a PDF of the results in the current working directory
pdf("T2 Cum Death Aug 23 Comparisons.pdf")
loop.cards %>% 
  ggplot(aes(x = target_end_date, y = wis, color = forecaster, group = forecaster)) +
  geom_line() + 
  geom_point() +
  ggtitle("Timepoint 2: WIS Comparison of Cumulative Deaths") + 
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  xlab("Target date for forecasts") + 
  ylab("Weighted Interval Score") +
  labs(caption = "Ground truth data: coronavirus.health.ny.gov") +
  theme_test()

loop.cards %>%
  ggplot(aes(x = target_end_date, y = spread, color = forecaster, group = forecaster)) +
  # ßgeom_line() +
  geom_point() +
  ggtitle("Timepoint 2: Spread Comparison of Cumulative Deaths") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  xlab("Target date for forecasts") +
  ylab("Spread") +
  labs(caption = "Ground truth data: coronavirus.health.ny.gov") +
  theme_test()

loop.cards %>% 
  ggplot(aes(x = target_end_date, y = ae, color = forecaster, group = forecaster)) +
  geom_line() + 
  geom_point() + 
  ggtitle("Timepoint 2: Absolute Error Comparison of Cumulative Deaths") + 
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  xlab("Target date for forecasts") + 
  ylab("Absolute Error") +
  labs(caption = "Ground truth data: coronavirus.health.ny.gov") +
  theme_test()

# This closes the PDF
dev.off()
###############################################################################

# Look at mean WIS, AE, spread
loop.cards %>%
  group_by(forecaster) %>%
  summarise(mean.wis = mean(wis, na.rm=T),
            mean.ae = mean(ae,na.rm=T))
