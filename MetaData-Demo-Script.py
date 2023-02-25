#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re

from MyFunctions import *


#%%

#======== Load CCM_Records.xlsx =======#
CCMRecords_path = f'CCM_Records_rarice.xlsx'
CCMRecords = pd.read_excel(CCMRecords_path, sheet_name=None)

#======== Load XRF Data ===============#
# Define XRF Data as a sheet in CCM Records
XRFdata = CCMRecords['XRF_Data']

# %%

performance_data = calculate_performance(CCMRecords_path, hfr_dir, ccm_group_names, V_oc)

independent = gen_independent(performance_data['id_dict'])

# Collect all parameters in df
cols = ['ccm_tags', 'independent','Anode Loading', 'Anode_Loading_stdev', 'OER_act','VHf_100', 'VHf_4000', 'HFRavg', 'HFRavg std', 'VHFRavg']
performance_comparison_df = pd.DataFrame(list(zip(performance_data['ccm_tags'], independent, anode_loading, anode_loading_stdev, performance_data['OER_act'], performance_data['VHF_dict']['100'], performance_data['VHF_dict']['4000'], performance_data['HFRavg_dict']['HFR_avg'], performance_data['HFRavg_dict']['HFR_stdev'], performance_data['HFRavg_dict']['VHFRavg'])), columns = cols)

#%%

y_colname_dict = {
    'Anode Loading':[r'Loading (mg$_{Ir}$/cm$^2$)', (0.2,0.65)],
    'OER_act':[r'OER Activity (mA/mg$_{Ir}$)', (0, 325)],
    'HFRavg':[r'HFR (m$\Omega cm^2$)', (105, 115)],
    'VHf_100':['$V_{HFRfree}$ ($mV$)', (1430, 1600)],
    }

#%%

plot_performance(performance_comparison_df, ind_label, y_colname_dict, performance_data['project'], (4,4), odir, True)