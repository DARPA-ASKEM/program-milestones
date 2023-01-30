Directory for TA1 delivarables.

- `xdd_query_reults`: JSON files that contain the results of the information seeking queries for scenario 1. This will be the basis for the text extraction pipelines for MIT and Arizona.

## Contact matrices

1. Search result for "SIR age" from xDD leads to https://www.nature.com/articles/s41598-021-94609-3 

2. Figure 1 comes from data in reference 15, "Prem, K. & Cook, A. R. Projecting social contact matrices in 152 countries using contact surveys and demographic data. PLoS Comput. Biol. 13, 20 (2017)." (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697). Figure 1 was converted into a CSV at figure1-10.1038-s41598-021-94609-3.csv. This article was independently found as one of the top hits via a direct xDD query (https://xdd.wisc.edu/api/articles?term=contact%20matrix%20population%20country&match=true&max=10&include_highlights=true)

3. The reference included additional supporting information: 
https://doi.org/10.1371/journal.pcbi.1005697.s001 PDF explaining the data.
https://doi.org/10.1371/journal.pcbi.1005697.s002 Zip file of multiple Excel spreadsheets

4. The Excel spreadsheets include, for ~150 countries, result of contact surveys.  For survey participants, the average number of individuals in each age bin (each age bin is 5 years, up to X=16) reported having contact with. The surveys additionally broke down data by contact location (home, school, work, other). The supporting PDF suggested reweighting the different location datasets to simulate interventions (e.g., reweight to 0 for school contacts to represent school closer, reweight to 0.5 the work and other interactions to represent social distancing).

5. We have uploaded data from USA (for all locations), India (as an example country with
   multi-generational contact) and Belgium (as an example country with no multi-generational contact). These are in files belgium_all_locations_cm.csv, usa_*_cm.csv and india_all_locations_cm.csv

6. We stopped here.  Note that potentially the weights could be played with for the UK dataset to see how close one can get to the UK during-pandemic surveys.

## Data from UK

1. h2020_cm_imputed.csv is from https://github.com/jarvisc1/comix_covid-19-first_wave and is used in
the following paper https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01597-8

