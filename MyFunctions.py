import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy.polynomial.polynomial as poly
import re
from pathlib import Path
from collections import defaultdict


def formatting():
    '''
    Returns 3x pair color palates, 5 color ibm palate, 8 color tol palate, and 5 symbol marker palate.
    Call:
        colors_pair, colors_ibm, colors_tol = color_pallate()
    '''
    # colorblind accessible pallates
    colors_tol = ['#332288', '#117733', '#44AA99', '#88CCEE',
                  '#DDCC77', '#CC6677', '#AA4499', '#882255',
                  '#332288', '#117733', '#44AA99', '#88CCEE',
                  '#DDCC77', '#CC6677', '#AA4499', '#882255',]
    
    colors_ibm = ['#648FFF', '#785EF0','#DC267F', '#FE6100', '#FFB000']
    
    colors_pair = [['#E1BE6A', '#40B0A6'], ['#FFC20A', '#0C7BDC'], ['#E66100', '#5D3A9B']]
    
    # Unique Markers
    markers_5 = ['v', 's', 'd', 'h', 'o',
                 'v', 's', 'd', 'h', 'o',]
    markers_8 = ['v', 's', 'd', 'h', 'o', '^', 'D', 'P',
                 'v', 's', 'd', 'h', 'o', '^', 'D', 'P',]
    
    # fontsizes
    fontsizes = {
        'title': 16,
        'ticks': 12,
        'label': 16,
        'legend': 14,
        }
    
    return(colors_pair, colors_ibm, colors_tol, markers_5, markers_8, fontsizes)

def calculate_VBA(J, V, HFR, ni, nf, V_oc):
    '''
    Given a polarization curve, HFR, and V_oc input with points for a Tafel fit,
    Inputs:
        * J = current density (A/cm2)
        * V = polarization curve voltage (V)
        * HFR = High Frequency Resistance (Ohm*cm2)
        * ni = initial Tafel fitting point 
        * nf = final Tafel fitting point
        * V_oc = open circuit voltage (V)
    
    returns 
        * vba_params = (b, m, V_HFRfree,
                        mt, V_HFR, V_kin, V_oc,
                        mt_diff, HFR_diff, kin_diff, OC_diff)
        * b = Tafel intercept
        * m = Tafel slope
        * V_HFRfree = High Frequency Resistance free Voltage (V)
        * mt = mass transport voltage loss (V)
        * V_HFR = 
        * V_kin = kinetic losses (V)
        * V_oc = open circuit potential of cell
        * _diff = difference between respective losses
    '''
    # note that log10 is used throughout, 
    # loge could also be used with identical result, but must be used consistently
    log10J = np.log10(J) # log(A/cm2)
    V_HFRfree = V - J * HFR # [V]
    V_oc_arr = np.full(len(log10J), V_oc) # array of V_oc, same shape as log10J

    # tafel linear regression
    b, m = poly.polyfit(log10J.iloc[ni:nf], V_HFRfree.iloc[ni:nf], 1, full=False)

    # since my V_kin are using log10J, all η are using log10J, not ln
    V_kin = (m * log10J) + b # m*log10(x)+b
    η_kin = V_kin - V_oc_arr
    η_ohmic = J * HFR
    η_mt = V_HFRfree - V_kin

    V_HFR = V_kin + η_ohmic
    mt = V_HFR + η_mt

    # take difference between losses on vba plot (area between curves), 
    mt_diff = mt - V_HFR
    HFR_diff = V_HFR - V_kin
    kin_diff = V_kin - V_oc_arr
    OC_diff = V_oc_arr
    
    return (b, m, V_HFRfree,
            mt, V_HFR, V_kin, V_oc,
            mt_diff, HFR_diff, kin_diff, OC_diff)




def monotonic(test_list):
    '''
    monotonic(test_list) takes a test_list and splits into monotonic groupings.
    
    returns:
    res = list of grouped values
    '''
    res = [] # result
    temp = [] 
    is_up = True
    if test_list[0] > test_list[1]:
        is_up = False
    for curr, nex in zip(test_list, test_list[1:]):
        temp.append(curr)

        # checking for increasing or decreasing to change list
        # if at trough or peak, complete current group, clear temp, add nex to temp
        if (nex > curr and not is_up) or (nex <= curr and is_up):
            res.append(temp)
            temp = []

            # toggling
            is_up = not is_up

    temp.append(nex)
    res.append(temp)
    
    return res


def process_cv_data(path):
    '''
    Given an experiment is carried out in the following order: CV/other measurements/CV
    such that there is a pre and post CV for that experiment, this function
    will plot paiarwise the 10th (or other specified) scan in both pre and post CVs
    according to the CCM id as extracted from file names. 
    
    File names are extracted assuming the following format: X-###_YYYY-MM-DD
    where X is the abbreviated project associated with the CCM, ### is the CCM number, 
    and YYYY-MM-DD is the year, month, and day of that experiment. The regex search key 
    can be changed in the extract CCM IDs section.
    
    At present this script assumes CV data will be in Autolab format
    ("Column 1 (V/s)", "Scan,Index", "Time (s)", "WE(1).Potential (V)", "WE(1).Current (A)").
    
    Input:
    path - folder path of relevant CV data
    
    Returns:
    cv_df_dict - key: CCM id, values: dataframe of processed CV data
    '''
    cv_df_dict = {}
    CCM_list = []
    CCMdate_list = []
    
    for l, file in enumerate(os.listdir(path)):
        
        # Extract CCM IDs
        # search for X-### CCM id format in file name
        CCM_match = re.search('\w-\d{3}', file)
        # extract match from result
        CCM = CCM_match.group(0)
        # search for X-###_####-##-## CCM id format in file name
        CCMdate_match = re.search('\w-\d{3}_\d{4}-\d{2}-\d{2}', file)
        # extract match from result
        CCMdate = CCMdate_match.group(0)

        # create CV dict of data to plot with
        temp_df = pd.read_csv(f'{path}/{file}')
        # extract only the 10th scan
        # change this to detect final scan - generalize
        temp_df = temp_df[temp_df['Scan'] == 10]
        # extract only Potential and Current Columns
        temp_df = temp_df.loc[:, ['WE(1).Potential (V)','WE(1).Current (A)']]
        # assign temp_df to dictionary according to CCM id
        cv_df_dict[CCMdate] = temp_df
        # collect IDs
        CCM_list.append(CCM)
        CCMdate_list.append(CCMdate)
        
    
    return (CCM_list, CCMdate_list, cv_df_dict)



def calculate_performance(CCMRecords_path, hfr_dir, ccm_group_names, V_oc):
    
    '''
    Inputs:
        * CCMRecords_path = complete path of CCMRecords.xlsx
        * hfr_dir = directory of hfr files
        * ccm_group_names = list of CCM IDs as strings
        * V_oc = array of V_oc constant
    
    Returns:
        * Dictionary with keys:
    ['ccm_tags', 'id_dict', 'project', 'catalyst_loading_dict', 'ind_label', 'OER_act', 'VHF_dict', 'HFRavg_dict']
    '''
    
    
    # Clean up lists, combine into dictionaries where possible
#     data = defaultdict(list)
    
    fontsize_ttl = 16
    fontsize_ticks = 12
    fontsize_lbl = 16
    fontsize_lgd = 14

    # this should be extracted from the data
    Jsp = [0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0]
    # initialize lists / dictionaries
    CCM_list = []

    # This dictionary will hold VHF at 100mA and 4000mA
    VHF_dict = defaultdict(list)
    # This dictionary will hold HFRavg, HFRavg_stdev, and VHFRavg
    HFRavg_dict = defaultdict(list)

    OER_act = []
    catalyst_loading_dict = {}
    ccm_tags = []
    ind_label = []

    # Load Full CCMRecords.xlsx
    CCMRecords = pd.read_excel(CCMRecords_path, sheet_name=None)
    # Define XRFdata_df
    XRFdata_df = CCMRecords['XRF_Data']


    for l, file in enumerate(os.listdir(hfr_dir)):

        # OREO CCM IDs use two letter 'OR'
        if 'OREO' in hfr_dir:
            ccmid_search = '\w\w-\d{3}'

        # All other projects use single letters ('H', 'F', etc.)
        else:
            ccmid_search = '\w-\d{3}'

        # Extract CCM ID from filename    
        CCM_match = re.search(ccmid_search, file)
        # extract match from result
        CCM = CCM_match.group(0)
        CCM_list.append(CCM)

        # Extract date from filename    
        date_match = re.search('\d{4}-\d{2}-\d{2}', file)
        # extract match from result
        date = date_match.group(0)

        # Extract Project Name
        # e.g., 'OR-001_008_KS', 
        project_match = re.search(ccmid_search + '_\d{3}_\w{2}', hfr_dir)
        project = project_match.group(0)

        #============== Import Loading Data, id_dict =================#

        # Select CCMRecords.xlsx sheet based on project
        if 'OREO' in hfr_dir:
            CCMRecords_project = CCMRecords['OREO_OR_series']
        elif 'H2New' in hfr_dir:
            CCMRecords_project = CCMRecords['H2New_H_series']

        #================ Define Independent Variable label =============#
        
        # Define df filter to slice present CCM row
        sample_filter = CCMRecords_project['project_id'] == CCM
        # Extract independent_variable for each CCM
        ind_label_i = CCMRecords_project[sample_filter]['independent_variable']
        # Convert to array to extract just value
        ind_label_i = np.array(ind_label_i)[0]
        # Collect labels
        ind_label.append(ind_label_i)

        # Define id_dict from CCMRecords.xlsx
        id_dict = dict(zip(list(CCMRecords_project['project_id']) , list(CCMRecords_project['description_abbrev'])))
        
        #===============================================================#
        
        # Slice XRF df to only relevant CCM samples
        CCM_filter = XRFdata_df['ccm_realtive_id'] == CCM
        XRF_df_temp = XRFdata_df[CCM_filter]
        
        # extract XRF data: Pt (cathode)
        Cathode_load = np.array(XRF_df_temp['mean_value_cathode (mg/cm2)'])[0]
        Cathode_load_stdev = np.array(XRF_df_temp['stdev_cathode'])[0]
        # extract XRF data: IrOx (anode)
        Anode_load = np.array(XRF_df_temp['mean_value_anode (mg/cm2)'])[0]
        Anode_load_stdev = np.array(XRF_df_temp['stdev_cathode'])[0]    
        # save tuple pair to dictionary (cathode, anode)
        catalyst_loading_dict[CCM] = ((Cathode_load, Cathode_load_stdev), (Anode_load, Anode_load_stdev))

        #============= Collect Sample Group Names ====================#

        # Collect grouping tags for each sample
        res = [ele for ele in ccm_group_names if(ele in id_dict[CCM])]
        # if res is not empty... (empty lists are false)
        if res:
            ccm_tags.append(res[0])

        #=============================================================#

        temp_hfr_df = pd.read_csv(f'{hfr_dir}/{file}')

        # I want to average HFR before calculating V_HFR
        # I dont want to use HFRavg for calculating VHFRfree though, must avg AFTER calcVBA

        # define points used for tafel linear regression
    #     ni = tafel_dict[file][0]
    #     nf = tafel_dict[file][1]
        ni = 3
        nf = 10

        J = temp_hfr_df['Current Density (A/cm2)'] # current density [A/cm2]  
        V = temp_hfr_df['Potential (V)'] # Voltage feedback [V]
        HFR = temp_hfr_df['HFR (Ωcm2)'] # high frequency resistance, extracted from EIS [Ωcm2]
        V_oc_arr = np.full(len(J), V_oc) # array of V_oc, same shape as J

        # ========== HFR Mean ===============# 
        pi = 0
    #     plt.plot(J[pi:], HFR[pi:], 'o-',mfc='none', label = id_dict[CCM])
    #     plt.legend(loc=(1.05,0))
    #     plt.title(f'HFR v. J')

        # ===================================#
        
        # Collect mean HFR
        HFRavg = HFR[pi:].mean()
        HFRavg_dict['HFR_avg'].append(HFRavg*1e3)
        
        # Collect stdev of mean HFR
        HFRavg_stdev_i = np.std(HFR)
        HFRavg_dict['HFR_stdev'].append(HFRavg_stdev_i*1e3)

        # Perform VBA calculation, assign outputs
        # for 10 prev. pcs, pass constant ni,nf = 3,10
        b, m, V_HFRfree, mt, V_HFR, V_kin, V_oc, mt_diff, HFR_diff, kin_diff, OC_diff \
        = calculate_VBA(J, V, HFR, ni, nf, V_oc)
        
        #=====================================================#

        #============ Calculate Catalyst Activity from V_HFRfree = 1.45V ========#
        
        # Define the desired voltage value to interpolate current density at
        # This voltage should be well within Tafel behavior regime
        tafel_yi = 1.45 # Volts
        
        # Write out the full relationship, and solve, noting the units
        # tafel_yi = m * log10tafel_xi + b
        # We have m, b, and a y_value (voltage) to interpolate at...
        # We want to solve for the x_value (current density)...
        log10tafel_xi = (tafel_yi - b) / m
        # correct those units
        # Base e interpolated current density
        log10_Ji_interpolated = log10tafel_xi
        # We want base 10 current density
        Ji_interpolated = 10 ** log10tafel_xi
        
        # Normalize extracted current density by anode loading
        # (A/cm2)/(mgIr/cm2) * 1e3 = mA/mgIr
        OER_acti = (Ji_interpolated / Anode_load) * 1e3
        
        # Collect OER value
        OER_act.append(OER_acti)

        # extract values at specific J
        val_extract_df = pd.DataFrame(list(zip(Jsp, V, V_HFRfree, HFR)), columns = ['Jsp', 'V', 'V_HFRfree', 'HFR'])

        # VHFRfree (VHf) _X in mA/cm2
        VHf_100_value = np.array(val_extract_df[val_extract_df['Jsp'] == 0.1]['V_HFRfree'])[0]
        VHF_dict['100'].append(VHf_100_value*1e3)

        VHf_4000_value = np.array(val_extract_df[val_extract_df['Jsp'] == 4.0]['V_HFRfree'])[0]
        VHF_dict['4000'].append(VHf_4000_value*1e3)

        # HFR
        # take average of HFR
        VHFRavg = V_HFR.mean()
        # collect avg
        HFRavg_dict['VHFRavg'].append(VHFRavg)
    
    # Return unique ind_label
    ind_label = np.unique(np.array(ind_label))[0]
    
    # filter catalyst_loading_dict to only ccms in ccm_list
    catalyst_loading_dict = {k: v for k, v in catalyst_loading_dict.items() if k in CCM_list}
    # filter id_dict likewise
    id_dict = {k: v for k, v in id_dict.items() if k in CCM_list}
    
    # Collect all performance data into a single output dictionary
    return_data_values = [ccm_tags, id_dict, project, catalyst_loading_dict, ind_label, OER_act, VHF_dict, HFRavg_dict]
    return_data_keys = ['ccm_tags', 'id_dict', 'project', 'catalyst_loading_dict', 'ind_label', 'OER_act', 'VHF_dict', 'HFRavg_dict']
    
    data = dict(zip(return_data_keys, return_data_values))
    
    return(data)



def unzip_loading(catalyst_loading_dict):
    '''
    Given a loading dictionary with keys (mean, stdev),
    unpacks the values into cathode / anode mean loading / stdev
    Input:
        catalyst_loading_dict
    Return:
        cathode_loading, cathode_loading_stdev, 
        anode_loading, anode_loading_stdev
    '''
    
    # Define the cathode and anode loading data lists
    cathode_loading_data, anode_loading_data = list(zip(*catalyst_loading_dict.values()))
    # Unpack loading data into loading and stdev lists
    cathode_loading, cathode_loading_stdev = zip(*cathode_loading_data)
    anode_loading, anode_loading_stdev = zip(*anode_loading_data)
    
    return(cathode_loading, cathode_loading_stdev, anode_loading, anode_loading_stdev)


def gen_independent(id_dict):
    '''
    Generates a list of xtick labels from id_dict
    Input: 
        id_dict
    Return:
        indepdent
    '''
    
    ind_list = []
    # Adds whitespace character to ticklabel so that each tick is unique
    # This is necessary for "#X" style ticklabels
    # If there are two '#1' samples assigned to different groups, 
    # they will share a tick if the ticklabels are identical
    for i, key in enumerate(list(id_dict.values())):
        ind_list.append(f'{" "*i}{list(id_dict.values())[i][-2:]}')
        
    independent = ind_list
        
    return(independent)


def plot_performance(performance_comparison_df, ind_label, y_colname_dict, project, figsize, odir, save_plot):
    
    '''
    Input: 
        * Performance Comparison dataframe
        * ind_label = label used for figure x-axis
        * y_colname_dict = column names and ylimits in the df to be plotted
        * project = name of project
        * figsize = figure size
        * odir = folder to save plot in
        * save plot = boolean to save the plot
    '''
    
    fontsize_ticks = 12
    fontsize_lbl = 14
    fontsize_lgd = 16
    
    # Define color palatte
    colors_ibm = ['#648FFF', '#785EF0','#DC267F', '#FE6100', '#FFB000']
    
    nrows = len(y_colname_dict.keys())
    # one column per group
    ncols = 1

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)

    # Define groups
    groups = list(performance_comparison_df['ccm_tags'].unique())
    
    # iterate through sample groups
    for j, group in enumerate(groups):

        # Filter df to one group at a time
        percomp_temp = performance_comparison_df[performance_comparison_df['ccm_tags'] == group]

        for i, ax in enumerate(axes.flat):

            # define plot formatting
            yi = [*y_colname_dict][i]
            lbli = y_colname_dict[yi][0]
            ylimi = y_colname_dict[yi][1]
            curves = [percomp_temp['independent'], percomp_temp[yi], 'o-']
            labels = ''

            if 'VHf' in yi:
                curves = [percomp_temp['independent'], percomp_temp['VHf_4000'], 'o-',
                          percomp_temp['independent'], percomp_temp[yi], 'o--']
                labels = ['J=4.0 A/cm$^2$', 'J=0.1 A/cm$^2$']
                
            if 'Loading' in yi:
                x = percomp_temp['independent']
                y = percomp_temp['Anode Loading']
                error = percomp_temp['Anode_Loading_stdev']
                
                ax.fill_between(x, y-error, y+error, alpha=0.25, color=colors_ibm[j])
                # There are three things being plotted for this loop
                # To label the fill_between I simply have to pass three labels, two empty
            
            # If only one group, color by row
            if len(groups) == 1:
                clr_iteration = i
            # If multiple groups color by group (column)
            else:
                clr_iteration = j
            
            ax.plot(*curves,
                    mfc='none', color=colors_ibm[clr_iteration])
            ax.legend(labels)
            ax.set_ylim(ylimi)
            lbl_bottom = False
            rot = 0
            # if last subplot, add xlabel
            if i == (nrows-1):
                lbl_bottom = True
                rot = 0
                ax.set_xlabel(ind_label, fontsize = fontsize_lbl)
            ax.tick_params('x', labelbottom=lbl_bottom, labelsize = fontsize_ticks, rotation=rot)
            ax.grid(alpha=0.35)
            ax.set_ylabel(lbli, fontsize=fontsize_lbl)
    
#     plt.suptitle(f'{groups}', fontsize=fontsize_lbl)
    plt.xticks(ha='right')
    
    if save_plot:
        plt.savefig(f'{odir}/{project}_PerformanceComp.svg', dpi=300)
        
        
#=====================================================================#

def plot_eis_pc(hfr_dir, id_dict, fontsizes, colors, markers, figsize, odir, save_pc):
    '''
    Plots polarization curves extracted from eis measurement
    Inputs:
        * hfr_dir = folder containing hfr data files in form J, V, HFR
        HFR folder has column names: 'Current Density (A/cm2)', 'Potential (V)', 'HFR (Ωcm2)'
        * id_dict = CCM ID paired with description
        * ft_dict = fontsize dictionary with keys: label, ticks, legend, title; values: integers
        * colors = color palate as list of hexidecimal strings
        * markers = markers to be used in plotting, list of strings
        * figusize = figure size
        * odir = folder path in which to save plots
        * save_pc = Boolean to save plot or not
    '''
    # Define list of CCMs 
    CCM_list = [*id_dict.keys()]
    
    # Iterate through files in HFR folder
    for l, file in enumerate(os.listdir(hfr_dir)):
        
        temp_hfr_df = pd.read_csv(f'{hfr_dir}/{file}')

        # Extract current density, voltage, and HFR
        Ji = temp_hfr_df['Current Density (A/cm2)']
        Vi = temp_hfr_df['Potential (V)']
        HFRi = temp_hfr_df['HFR (Ωcm2)']
        # Calculate HFR free voltage
        VHFRfreei = Vi - HFRi * Ji
        
        
        # Plot Polarization Curve
        plt.figure('polcurve', figsize = figsize)
        plt.plot(Ji, Vi, 
                     marker = markers[l], ls = '-', mfc='none', 
                     color = colors[l],
                     label = f'{id_dict[CCM_list[l]]}')
        plt.plot(Ji, VHFRfreei, 
                     marker = markers[l], ls = (0, (1, 1)), mfc='none', 
                     color = colors[l])

        plt.ylabel('Voltage (V)', fontsize = fontsizes['label'])
        plt.xlabel('Current Density (A/cm$^2$)', fontsize = fontsizes['label'])
        plt.tick_params(axis='both', labelsize = fontsizes['ticks'])
        plt.grid(alpha = 0.35)
        plt.legend(loc = 'best', fontsize = fontsizes['legend'])
        plt.tight_layout()

        if save_pc:
            fname_pc_raw = 'pc_raw.svg'
            plt.savefig(f'{odir}/{fname_pc_raw}')
        
        # Plot Tafel Curve (HFR corrected polarization curve)
        plt.figure('vhfr_free', figsize = (7, 5))
        plt.semilogx(Ji, VHFRfreei, 
                     marker = markers[l], ls = (0, (1, 1)), mfc='none', 
                     color = colors[l], 
                     label = f'{id_dict[CCM_list[l]]}')

        plt.ylabel('HFR-free Voltage (V)', fontsize = fontsizes['label'])
        plt.xlabel('Current Density (A/cm$^2$)', fontsize = fontsizes['label'])
        plt.tick_params(axis='both', labelsize = fontsizes['ticks'])
        plt.grid(alpha = 0.35)
        plt.legend(loc = 'best', fontsize = fontsizes['legend'])
        plt.tight_layout()

        if save_pc:
            fname_pc_hfr_corrected = 'pc_hfr_corrected.svg'
            plt.savefig(f'{odir}/{fname_pc_hfr_corrected}')
            
# ======================================================== #

#====================================================================#
#========== Functions to Parse Gamry data into Autolab format =======#
#====================================================================#

def line_num_for_phrase_in_file(filename, phrase):
    '''
    Given a phrase, returns all lines in file where that phrase is found \n
    Inputs:
        * filename = path to the file in question
        * phrase = string to be searched for in file \n
    Returns:
        * phrase_rows = list of rows where phrase is found in file
    '''
    line_nums = []
    # Using with to open files is the pythonic idiom as it ensures the file will be properly closed when the block using the file ends. 
    # There are special characters in the DTA files which require a specific encoding to read
    with open(filename, encoding='ISO-8859-1') as file:
        # iterate through lines in the file
        for row, line in enumerate(file, 1):
            # search each line for the phrase
            if phrase in line:
                line_nums.append(row)
    return(line_nums)

#==============================================================#

def dta_to_df(file, phrase = '	Pt	'):
    '''
    Inputs:
        * file = relative path to DTA file (not directory) in question
        * skiprows = number of rows in file to skip before reading data \n
            - this should be the line above table headers \n
            - skiprows variable passed to pandas read_csv function \n
    Returns:
        * df = extracted dataframe from DTA file
    '''
    # The unique phrase to be searched for is the beginning column headers line
    # This will be used to parse the txt file and load the below table
    table_location = line_num_for_phrase_in_file(file, phrase)
    # Read table at table location
    df = pd.read_table(file, skiprows = table_location[0]-1, delimiter='\t', header=0, encoding='ISO-8859-1')
    # Drop unnamed first column
    df = df.drop(columns='Unnamed: 0')
    if 'Over' in df:
        df = df.drop(columns='Over')
    df = df.drop(labels=0, axis=0)
    df = df.reset_index(drop=True)
    inds = []
    for i in range(len(df['Pt'])):
        if df['Pt'][i]=='TABLE':
            inds.append(i)
            inds.append(i+1)
            inds.append(i+2)
    df = df.drop(df.index[inds])
    
    df = df.reset_index()
    # del df['index']
    df = df.drop(['Pt'], axis=1)
    df = df.astype(float)
    return(df)

#==============================================================#

def parse_gamryEIS(dir_EIS):
    # create list for eis filenames
    EIS_filenames = [] 
    # create list for polcurve filenames
    polcurve_filenames = [] 

    for file in os.listdir(dir_EIS):
        # Sort out the relevant files in the directory
        if 'EIS' in file:
            EIS_filenames.append(file) # add EIS filename to list
        if 'polcurve' in file:
            polcurve_filenames.append(file) # add polcurve filename to list

    # create list for EIS dfs
    EIS_dfs = []
    # create list for polcurve dfs
    polcurve_dfs = []

    # for each 
    for dat in EIS_filenames: 
        EIS_dfs.append(dta_to_df(f'{dir_EIS}/{dat}')) # load / parse EIS df, add to df list
    for dat in polcurve_filenames:
        polcurve_dfs.append(dta_to_df(f'{dir_EIS}/{dat}')) # load / parse polcurve df, add to list

    EIS_df = pd.concat(EIS_dfs) # combine all EIS dfs (single experiment) into single df
    EIS_df = EIS_df.reset_index() # reset implicit index for concatenated df
    polcurve_df = pd.concat(polcurve_dfs) # likewise for polcurve dfs
    polcurve_df = polcurve_df.reset_index() # likewise for polcurve concatenated df

    column_names_dict = {
        'Vdc':'WE(1).Potential (V)',
        'Idc':'WE(1).Current (A)',
        'level_0':'Index',
        'Freq':'Frequency (Hz)',
        'Zreal':"Z' (Ω)",
        'Zimag':"Z'' (Ω)",
        'Zmod':"Z (Ω)",
        'Zphz':'Phase (°)',
        'Time':'Time (s)'
    }
    # Autolab Columns: WE(1).Potential (V),WE(1).Current (A),Index,Frequency (Hz),Z' (Ω),-Z'' (Ω),Z (Ω),-Phase (°),Time (s)
    # 'Pt', 'Time', 'Freq', 'Zreal', 'Zimag', 'Zsig', 'Zmod', 'Zphz', 'Idc', 'Vdc', 'IERange'

    # select only necessary columns
    EIS_df = EIS_df[['Vdc', 'Idc','level_0', 'Freq', 'Zreal', 'Zimag', 'Zmod', 'Zphz', 'Time']]
    # rename columns to match Autolab format
    EIS_df = EIS_df.rename(column_names_dict, axis="columns")
    # make Z'' --> -Z''
    EIS_df["-Z'' (Ω)"] = EIS_df["Z'' (Ω)"]*(-1)
    EIS_df = EIS_df.drop("Z'' (Ω)", axis=1)

    return(EIS_df)