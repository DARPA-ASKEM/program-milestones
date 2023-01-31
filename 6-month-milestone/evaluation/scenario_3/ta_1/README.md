## Datasets for Scenario 3

Google doc tracking methodology for finding datasets at https://docs.google.com/document/d/1Xd1ug-7yg2hxikBr6HQ8b9PBDOeN9AxuvTBa72RZkPw/edit#


### Population Data
1. US population data (2021) downloaded from US Census Bureau. Checked in as usa-2021-population-age-stratified.csv 


### Data from Google Health Data
Steps followed
1. Followed Holly's link from Slack
2. Found a number of relevant datasets at Google Health data - https://health.google.com/covid-19/open-data/raw-data
3. Parsed CSVs to extract USA specific data for cases, deaths, hospitalizations, vaccinations and
   by-age stratified data.
4. These are present in the google-health-data directory


### Rates
1. Hospitalization rate over time. Parsed using CDC's case surveillance data
   (https://data.cdc.gov/Case-Surveillance/COVID-19-Case-Surveillance-Public-Use-Data/vbim-akqf) and then processed
   using datatable to group each case into the month it was reported and whether it resulted in a
   hospitalization. Checked in as hospitalization-rate-by-month.csv  

