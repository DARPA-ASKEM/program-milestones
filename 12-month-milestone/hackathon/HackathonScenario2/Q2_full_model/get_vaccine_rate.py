# returns vaccination rate and dose delay
import numpy as np
from datetime import date
    
def get_vaccine_rate(t): 
    day0 = date(2020, 1, 1)
    t_vac0 = (date(2020, 12, 14) - day0).days #t_vac0 = daysact(day0,'14-dec-2020')
    t_delay = (date(2021, 5, 31) - day0).days #t_delay = daysact(day0,'31-may-2021')
    t_vac_peak = (date(2021, 6, 30) - day0).days #t_vac_peak = daysact(day0,'30-jun-2021')
    
    
    # 55.2# Ontarian fully vaccinated, 68.7# at least one dose
    fraction_first_dose = 0.6
    
    peak_vac_rate = 220 * 5e-05
    
    ss_vac_rate = 10 * 5e-05
    
    # Pfizer vs AstraZeneca
    eps_PZ = 0.95
    eps_AZ = 0.05
    # eps_PZ = 0.1;
    # eps_AZ = 0.9;
    
    if (t < t_vac0): # no vaccine for the first year
        vr1 = np.array([[0],[0]])
        vr2 = np.array([[0],[0]])
    elif (t < t_delay): # 16-week dose delay
        # assume vaccination rate linearly increases from 0 on 1/1/2021 to 6/30/2021
        curr_vac_rate = peak_vac_rate * (t - t_vac0) / (t_vac_peak - t_vac0)
        # 95# Pfizer, 5# AZ
        vr1 = curr_vac_rate * np.array([[eps_PZ],[eps_AZ]])
        vr2 = np.array([[1 / (16 * 7)],[1 / (16 * 7)]])
    elif (t < t_vac_peak): # 3-week dose delay
        # assume vaccination rate linearly increases from 0 on 1/1/2021 to 6/30/2021
        curr_vac_rate = peak_vac_rate * (t - t_vac0) / (t_vac_peak - t_vac0)
        # 95# Pfizer, 5# AZ
        vr1 = curr_vac_rate * np.array([[eps_PZ],[eps_AZ]])
        vr2 = np.array([[1 / 28],[1 / (8 * 7)]])
    else:
        curr_vac_rate = (peak_vac_rate - ss_vac_rate) * np.exp(- 0.03 * (t - t_vac_peak)) + ss_vac_rate
        # 95# Pfizer, 5# AZ
        vr1 = curr_vac_rate * np.array([[eps_PZ],[eps_AZ]])
        vr2 = np.array([[1 / 28],[1 / (8 * 7)]])
    
    return vr1,vr2