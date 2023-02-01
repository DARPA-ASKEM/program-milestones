import pandas as pd

#open another dataset
df = pd.read_csv('usa-IRDVHN_age.csv')

#remove the second column which is an unnamed column
df = df.drop(columns=['Unnamed: 0'])

#get the first row of the population breakdowns
age_groups = ['N_0-9', 'N_10-19', 'N_20-29', 'N_30-39', 'N_40-49', 'N_50-59', 'N_60-69', 'N_70-79', 'N_80-89', 'N_90-99', 'N_100+']
pop_breakdown = df[age_groups].iloc[0].to_numpy()
total_pop = pop_breakdown.sum()
pop_percentages = pop_breakdown / total_pop



#create new columns for hospitalizations by age group
for i, group in enumerate(age_groups):
    d_name = 'D_' + group[2:]
    df[d_name] = (df['D'] * pop_percentages[i]).astype(int)

#do the same for deaths
for i, group in enumerate(age_groups):
    h_name = 'H_' + group[2:]
    df[h_name] = (df['H'] * pop_percentages[i]).astype(int)


#create columns for hospitalizations by vaccination status
df['H_vac'] = (df['V'] / total_pop * df['H']).astype(int)
df['H_unvac'] = (df['H'] - df['H_vac']).astype(int)

#create columns for deaths by vaccination status
df['D_vac'] = (df['V'] / total_pop * df['D']).astype(int)
df['D_unvac'] = (df['D'] - df['D_vac']).astype(int)



#create columns for hospitalizations/deaths by vaccination status and age group (e.g. H_vac_0-9, H_vac_10-19, ..., H_vac_90-99, H_vac_100+, H_unvac_0-9, H_unvac_10-19, ..., H_unvac_90-99, H_unvac_100+)
for i, group in enumerate(age_groups):
    df['H_vac_' + group[2:]] = (df['H_vac'] * pop_percentages[i]).astype(int)

for i, group in enumerate(age_groups):
    df['H_unvac_' + group[2:]] = (df['H_unvac'] * pop_percentages[i]).astype(int)

for i, group in enumerate(age_groups):
    df['D_vac_' + group[2:]] = (df['D_vac'] * pop_percentages[i]).astype(int)

for i, group in enumerate(age_groups):
    df['D_unvac_' + group[2:]] = (df['D_unvac'] * pop_percentages[i]).astype(int)



df.to_csv('usa-IRDVHN_age_HD_breakdown.csv', index=False)