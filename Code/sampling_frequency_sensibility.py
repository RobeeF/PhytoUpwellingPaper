# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 12:22:19 2022

@author: rfuchs
"""


import os
import pickle
import numpy as np
import pandas as pd
import ruptures as rpt
from os.path import join
from copy import deepcopy
import matplotlib.pyplot as plt
from ruptures.costs import CostLinear

# Change the path with yours:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper')


dates_col = ['WUI_start', 'WUI_end', 'T_start', 'T_end', 'Tmax', 'Tmin',
       'window_start', 'window_end', 'ORGNANO_start',
       'ORGNANO_end', 'ORGPICOPRO_start',
       'ORGPICOPRO_end', 'REDNANO_start', 'REDNANO_end',
       'REDPICOEUK_start',
       'REDPICOEUK_end', 'REDPICOPRO_start', 'REDPICOPRO_end']

#############################################################################
# Temperature and phyto 
#############################################################################

# Import Upwelling and Downwelling dates 
with open(join('Results','date_pickles','welling_events.pickle'), 'rb') as f1:
     welling = pickle.load(f1)

with open(join('Results','date_pickles','stratify.pickle'), 'rb') as f1:
     stratif_dict = pickle.load(f1)
     
stps = pd.read_csv(join('Data', 'INSITU','stps_1H.csv'),  sep = ',', parse_dates=['date'],\
                   index_col = 'date')

wui = pd.read_csv(join('Data', 'WUI','WUI_1H.csv'), parse_dates=['date'], index_col = 'date')
wui_index = 'marseille'

# Import also the T data
wui_sst = wui.join(stps, how="inner")
wui_sst = wui_sst[['marseille', 'T']]


#=======================================================
# Reference values
#=======================================================

path = join('Results', 'responses', 'abundance' + '.csv')
ref_abundance = pd.read_csv(path, parse_dates = dates_col)   
path = join('Results', 'responses', 'biomass' + '.csv')
ref_biomass = pd.read_csv(path, parse_dates = dates_col)   


#==============================================================
# Rupture detection for abundances, biomass for different lags
#==============================================================

# Number of breakpoints to look for rupture detection
n_bkps = 2  
pfg_entities = ['abundance', 'biomass']
pfg_colors = dict(zip(['Orgnano', 'Orgpicopro', 'Rednano', 'Redpicoeuk', 'Redpicopro'],\
                      ['#9467bd', '#1f77b4', '#2ca02c', '#d62728', 'orange']))

lag_res = {'abundance': {}, 'biomass': {}} 
lag_res['abundance']['2'] = ref_abundance
lag_res['biomass']['2'] = ref_biomass

for pfg_entity in pfg_entities:
    for lag in [2, 3, 4, 5]:        
        
        event_df = pd.read_csv(join('Results','event_T_anomalies_suited.csv'),\
                               parse_dates=['WUI_start', 'WUI_end',\
                               'T_start', 'T_end', 'Tmax', 'Tmin',\
                               'window_start', 'window_end'])
            
        #====================================
        # Data importation 
        #====================================
    
        phyto = pd.read_csv(join('Data','INSITU','PFG_' + pfg_entity + '.csv'),\
                            parse_dates = ['date'], index_col = 'date')    
        pfg_list = phyto.columns
        
        event_df['is_enough_phyto_data'] = True   
        # Create the columns to store the dates of the pfg reactions
        for pfg in pfg_list:
            event_df[pfg + '_start'] = pd.NaT
            event_df[pfg + '_end'] = pd.NaT
            event_df[pfg + '_phyto_react_before_T'] = np.nan
            event_df[pfg + 'median_before'] = np.nan
            event_df[pfg + 'median_during'] = np.nan
            event_df[pfg + 'median_after'] = np.nan
             
        #====================================
        # Iterate through the events
        #====================================
        
        for event_idx, event in event_df.iterrows():
            event_data = phyto.loc[event['window_start']:event['window_end']]
            event_data = event_data.loc[::lag]
            
            # If one phase do not have enough data, then it is not considered
            if (pd.isna(event_data).mean() >= 0.3).all():
                event_df.loc[event_idx, 'is_enough_phyto_data'] = False
                continue
            
            if len(event_data) <= 2:
                #print('No data for event:', event['window_start'], event['window_end'])
                continue
        
            # Interpolate the data
            event_data = event_data.interpolate(method = 'polynomial', order = 1)
    
            #====================================
            # Iterate through the PFGs
            #====================================
        
            for pfg in pfg_list:
                
                # Rupture detection (pfg_after_Tdrop useless to clear)
                pfg_data = event_data[[pfg]].dropna()
                
                if len(pfg_data) == 0:
                    #print(event['window_start'], ': no data for', pfg)
                    event_df.loc[event_idx, pfg + '_phyto_react_before_T'] = True
                    continue
    
                #====================================
                # Choose the kernel 
                #====================================
    
                pfg_data = pfg_data.join(wui_sst['T'])
                pfg_data['intercept'] = 1
                algo = rpt.Binseg(custom_cost = CostLinear, jump = 1).fit(pfg_data.values)
    
                result = algo.predict(n_bkps = n_bkps)
                result[-1] = result[-1] - 1 # Odd pattern of the package
                result_index =  pfg_data.index[result] # the date index 
                
                # Characterize the behaviour
                if len(result) != n_bkps + 1:
                    print('Lag:', lag, 'Entity', pfg_entity, 'pfg', pfg, 'I break')
                    continue
                    
                #assert len(result) == n_bkps + 1
                event_df.loc[event_idx, pfg + '_start'] = result_index[0]
                event_df.loc[event_idx, pfg + '_end'] = result_index[1]
                
                # Compute the mean of the pfg before and after
                median_before = pfg_data.loc[pfg_data.index < result_index[0], pfg].median()
                median_during = pfg_data.loc[(pfg_data.index >= result_index[0]) \
                                       & (pfg_data.index < result_index[1]), pfg].median()
                median_after = pfg_data.loc[pfg_data.index >= result_index[1], pfg].median()
        
                # Store them
                event_df.loc[event_idx, pfg + 'median_before'] = median_before
                event_df.loc[event_idx, pfg + 'median_during'] = median_during
                event_df.loc[event_idx, pfg + 'median_after'] = median_after
                
                # Check that the change does not occur before the event
                if result_index[0] < event['Tmax']:
                    event_df.loc[event_idx, pfg + '_phyto_react_before_T'] = True
                else:
                    event_df.loc[event_idx, pfg + '_phyto_react_before_T'] = False
            
        #========================================
        # Delete events which have no phyto data 
        #========================================
    
        event_df = event_df[event_df['is_enough_phyto_data']].reset_index(drop = True)
        lag_res[pfg_entity][str(2 + lag * 2)] = deepcopy(event_df)


#==============================================================
# Compute stat diff with lags
#==============================================================

responses = {'abundance': {}, 'biomass': {}} 
for pfg_entity in pfg_entities:
    for lagH, event_df in lag_res[pfg_entity].items():
        
        #***************************************
        # Data importation 
        #***************************************
        #path = join('Results', 'responses', pfg_entity + '.csv')
        #event_df = pd.read_csv(path, parse_dates = dates_col)     
        
        #***************************************
        # Store the quantities of interest
        #***************************************
        
        response = pd.DataFrame(columns = ['reaction_delay', 'reaction_duration', 'reaction_magnitude',\
                                                       'relaxation_magnitude', 'pfg', 'WUI_start']) 
        for pfg in pfg_list:
            # Trigger length
            reaction_delay = (event_df[pfg + '_start'] - event_df['Tmax']).apply(lambda x: x.total_seconds() / (3600 * 24))
            reaction_delay = np.where(event_df[pfg + '_phyto_react_before_T'], np.nan, reaction_delay)
        
            # Trigger time
            reaction_duration =  (event_df[pfg + '_end'] - event_df[pfg + '_start']).apply(lambda x: x.total_seconds() / (3600 * 24))
                
            # Reaction_strength
            reac_magn =  ((event_df[pfg + 'median_during'] / event_df[pfg + 'median_before']) - 1) * 100
            relax_magn =  ((event_df[pfg + 'median_after'] / event_df[pfg + 'median_during']) - 1) * 100
        
            df = pd.DataFrame(data = np.stack([reaction_delay, reaction_duration, reac_magn, relax_magn]).T,\
                              columns = ['reaction_delay', 'reaction_duration', 'reaction_magnitude', 'relaxation_magnitude'])
            
            nan_indices = event_df[pfg + '_phyto_react_before_T'].fillna(True)            
            df.loc[nan_indices] = np.nan
            
            df['pfg'] = pfg
            df['WUI_start'] = event_df['WUI_start']
            
            response = pd.concat([response, df])
                
        response = response[~response.isna().any(1)]
        response['pfg'] = response['pfg'].apply(lambda x: str.capitalize(x))
    
        responses[pfg_entity][lagH] = deepcopy(response)


#==============================================================
# Compute the difference
#==============================================================

diffs = pd.DataFrame()
for pfg_entity in pfg_entities:
    ref = responses[pfg_entity]['2'].set_index(['pfg', 'WUI_start'])
    for lagH, resp in responses[pfg_entity].items():
        #if lagH == '0H':
            #continue
        resp = resp.set_index(['pfg', 'WUI_start'])
        
        missing_indices = set(ref.index) - set(resp.index) 
        if missing_indices:
            missing_data = pd.DataFrame(index = list(missing_indices))
            resp = pd.concat([resp, missing_data])
            
        resp = resp.loc[ref.index]
        
        diff = pd.DataFrame((np.abs((ref - resp) / ref)).mean()).T
        diff['lag'] = lagH
        diff['pfg_entity'] = pfg_entity
        diffs = pd.concat([diffs, diff])

diffs.columns = ['Reaction delay', 'Reaction duration', 'Reaction magnitude',
       'Relaxation magnitude', 'lag', 'pfg_entity']


for pfg_entity, diff_data in diffs.groupby(['pfg_entity']):
    fig, ax = plt.subplots(figsize = (7, 7))
    diff_data = diff_data.set_index('lag').iloc[:,:4]  * 100 # In percentage
    ax.plot(diff_data)
    ax.set_ylabel('Variation wrt. 2H-data estimates (%)')
    ax.set_xlabel('Data frequency (H)')

    plt.legend(diff_data.columns)
    plt.tight_layout()
    fig_path = join('Figures', 'sensitivity', 'sensitivity_' + pfg_entity + '.png')
    plt.savefig(fig_path)
    plt.show()
        
    