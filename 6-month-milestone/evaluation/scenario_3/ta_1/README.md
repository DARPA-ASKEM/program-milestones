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
2. Using same dataset, also stratified by age group and generated age-stratified-by-month-hosp-cases-deaths.csv.
   Code for both parts is in cdc-case-data/process.py
3. Hazard rates conditioned on vaccination is derived from https://data.cdc.gov/Public-Health-Surveillance/Rates-of-COVID-19-Cases-or-Deaths-by-Age-Group-and/3rge-nu2a and checked in as vaccination-hazard-rates-age-month.csv. Note: This data has also be stratified by age.

4. Hospitalization rates conditioned on vaccination. From Brian:
   "Using California data from May - Nov 2021: unvaccinated individuals had a 12.7x increased rate of hospitalization compared to vaccinated."
   This is also from https://doi.org/10.15585/mmwr.mm7104e1. Method: Queried "covid vaccine hazard rates" in Terarium


### Reinfection data

1. From Brian:  For California, median time to reinfection is 262 days for vaccinated and 277 days for unvaccinated.
Data was from May - Nov 2021. Source: https://doi.org/10.15585/mmwr.mm7104e1 Method: Put in query "covid reinfection hazard rates" in Terarium


### Sero prevalence data
1. From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9716971/pdf/main.pdf 
"From November 4, 2020, to January 27, 2021, the ratio was 2.8 (CI: 2.8–2.9) and then decreased, reaching a low point from April 21 to July 1, 2021 (1.1, CI: 0.6–1.7). These ratios increased to 2.3 (CI: 2.0–2.5) from July 1 to September 20, 2021, held steady from September 20 to December 8, 2021 (2.2, CI: 2.0–2.5), and increased again from December 8, 2021, to February 26, 2022 (3.1, CI: 3.0–3.3)."

The change ratios, ratios estimating the change in seroprevalence compared to the change in reported case prevalence, can be used as a multiplier to enhance the understanding of the infection burden represented by officially reported case rates

### Text Extractions
- `scenario3_skema_extractions.xlsx`: Parameter values and descriptions mined from relevant papers identified by TA1, which may inform TA2-TA3 
