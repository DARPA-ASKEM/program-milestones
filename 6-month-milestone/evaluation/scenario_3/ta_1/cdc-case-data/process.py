from datatable import dt, f, g, by, join, sort, update, fread

# File can be downloaded from https://data.cdc.gov/Case-Surveillance/COVID-19-Case-Surveillance-Public-Use-Data/vbim-akqf
# NOTE: This is a 10G file and can be slow to download

all_cases = dt.fread("./CDC_COVID_Cases.csv")

# Get month that case was reported
all_cases[:, dt.update(case_month = dt.str.slice(f.cdc_case_earliest_dt, 0, 7))]

hosp_cases = all_cases[f.hosp_yn == "Yes", :]
death_cases = all_cases[f.death_yn == "Yes", :]

# Group by case_month and count total cases, hosp, deaths
cases_by_month = all_cases[:, {"num_cases": dt.count()} , by("case_month")]
hosp_by_month = hosp_cases[:, {"num_hosp": dt.count()} , by("case_month")]
death_by_month = death_cases[:, {"num_death": dt.count()} , by("case_month")]

hosp_by_month.key = "case_month"
death_by_month.key = "case_month"

join1 = cases_by_month[:, [f.case_month, f.num_cases, g.num_hosp], join(hosp_by_month)]
aggregate_joined = join1[:, [f.case_month, f.num_cases, f.num_hosp, g.num_death], join(death_by_month)]
aggregate_joined.to_csv("./all_month_hosp_cases_deaths.csv")

# Repeat above steps but with age stratification
cases_by_age_month = all_cases[:, {"num_cases": dt.count()}, by("age_group", "case_month")]
hosp_by_age_month = hosp_cases[:, {"num_hosp": dt.count()}, by("age_group", "case_month")]
death_by_age_month = death_cases[:, {"num_death": dt.count()}, by("age_group", "case_month")]

hosp_by_age_month.key = ['age_group', 'case_month']
death_by_age_month.key = ['age_group', 'case_month']

join1 = cases_by_age_month[:, [f.age_group, f.case_month, f.num_cases, g.num_hosp], join(hosp_by_age_month)]

aggregate_joined = join1[:, [f.age_group, f.case_month, f.num_cases, f.num_hosp, g.num_death], join(death_by_age_month)]

aggregate_joined.to_csv("./age_month_hosp_cases_deaths.csv")
