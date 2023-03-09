# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 15:33:00 2021

@author: rfuchs
"""

import os
import pickle
import numpy as np
import pandas as pd
import ruptures as rpt
from os.path import join

from ruptures.costs import CostLinear

from scipy.signal import filtfilt, butter

# Change the path with yours:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper/Code')
from utilities import split_sequences, merge_close_events, is_stratified

os.chdir('../Data')

#===================================================================
# Extract the stratified periods
#===================================================================

# Files importation
stps = pd.read_csv(join('INSITU','stps_1H.csv'),  sep = ',')
stps['date'] = pd.to_datetime(stps['date'], format = "%Y-%m-%d %H:%M:%S")
stps = stps.set_index('date')
del(stps['Salinity'])

meanT = stps['T'].mean() # Mean Temperature over the period

# Filter the series
b, a = butter(1, 0.0007)
stps['filter'] = filtfilt(b, a, x=stps['T'], method='gust')

# Indentify the the regime shift dates
from_strat_to_unstrat = np.where((stps['filter'] > meanT).astype(int).diff() == -1)[0]
from_unstrat_to_strat = np.where((stps['filter'] > meanT).astype(int).diff() == 1)[0] 

strat_shift = np.hstack([from_strat_to_unstrat, from_unstrat_to_strat]) 
strat_shift.sort()

strat_shift_dates = [stps.iloc[strat_event].name for strat_event in strat_shift]
strat_shift_dates.sort()

strat_shift_dates.insert(0, stps.iloc[0].name)
strat_shift_dates.append(stps.iloc[-1].name)


end_date = pd.to_datetime('2021-11-20')
strat_dict = dict()
strat_dict['Stratified'] = [[strat_shift_dates[0], strat_shift_dates[1]],\
                            [strat_shift_dates[2], strat_shift_dates[3]],\
                            [strat_shift_dates[4], strat_shift_dates[5]]]
    
strat_dict['Unstratified'] = [[strat_shift_dates[1], strat_shift_dates[2]],\
                              [strat_shift_dates[3], strat_shift_dates[4]],
                              [strat_shift_dates[5], end_date]]

with open(join('..','Results','date_pickles','stratify.pickle'), 'wb') as f1:
    pickle.dump(strat_dict, f1)

#****************************
# Extract Upwelling- Downwelling events
#****************************
# Import Meteorological data
wui = pd.read_csv(join('WUI','WUI_1H.csv'), parse_dates = ['date'], index_col = 'date')

dates = wui[(wui.index >= pd.to_datetime('2019-09-18 14:00:00')) & \
                    (wui.index < pd.to_datetime('2021-11-01 00:00:00'))]

#=======================================
# Upwelling
#=======================================  
valid_dates = dates[dates['marseille'] >= 0]['marseille']
event_dates_list = split_sequences(valid_dates)
    
sequences_mean_wui = [np.mean(wui.loc[event[0]:event[-1]]['marseille'])\
                  for event in event_dates_list] # To check

q_mean_up = 0.433 # From Odic et al.

upwelling_events = [[event_dates_list[i][0], event_dates_list[i][-1]] for i in range(len(event_dates_list)) if
                    sequences_mean_wui[i] >= q_mean_up]

# Merge the events
welling = merge_close_events(upwelling_events)
nb_events = len(welling)

with open(join('..', 'Results', 'date_pickles','welling_events.pickle'), 'wb') as f1:
    pickle.dump(welling, f1)  

#############################################################################
# Temperature and phyto 
#############################################################################

# Import Upwelling and Downwelling dates 
with open(join('..','Results','date_pickles','welling_events.pickle'), 'rb') as f1:
     welling = pickle.load(f1)

with open(join('..','results','date_pickles','stratify.pickle'), 'rb') as f1:
     stratif_dict = pickle.load(f1)
     
stps = pd.read_csv(join('INSITU','stps_1H.csv'),  sep = ',', parse_dates=['date'],\
                   index_col = 'date')

wui = pd.read_csv(join('WUI','WUI_1H.csv'), parse_dates=['date'], index_col = 'date')
wui_index = 'marseille'

# Import also the T data
wui_sst = wui.join(stps, how="inner")
wui_sst = wui_sst[['marseille', 'T']]

#===================================================================
# Extract Upwellings in stratified period
#===================================================================
    
# Compute T anomaly
b, a = butter(1, 1 / (24 * 15))
wui_sst['T_anomaly'] = wui_sst['T'] - filtfilt(b, a, x=wui_sst['T'], method='gust')
wui_sst['T_filted'] = filtfilt(b, a, x=wui_sst['T'], method='gust')


neg_anomalies =  wui_sst[wui_sst['T_anomaly'] <= 0]['T_anomaly']
neg_anomalies = split_sequences(neg_anomalies)

pos_anomalies =  wui_sst[wui_sst['T_anomaly'] >= 0]['T_anomaly']
pos_anomalies = split_sequences(pos_anomalies)

event_df = pd.DataFrame()

margin = '8H'

for event_idx, event in enumerate(welling):
    is_strat = is_stratified(event[0], stratif_dict)
    
    # The T anomalies cannot start before the WUI event or 
    # more than 24h after the end of the WUI event
    anomalies_ev = [[ano[0], ano[-1]] for ano in neg_anomalies if
                        event[0] <= ano[0] <= event[1] + pd.Timedelta('24H')]

    # Delete artifact anomalies
    anomalies_ev = [ano for ano in anomalies_ev if ano[-1] - ano[0] >= pd.Timedelta(margin)]
    
    if len(anomalies_ev) > 0:
        anomalies_ev = merge_close_events(anomalies_ev)

    else:
        continue
           
    data = wui_sst[(wui_sst.index  >= event[0] - pd.Timedelta('2D')) & \
                   (wui_sst.index  <= max(anomalies_ev[-1][-1], event[-1]))]
        
    Tmax = data[(data.index  >= event[0]) & (data.index  <= anomalies_ev[0][0])]['T'].idxmax()
    Tmin = data[(data.index  >= anomalies_ev[0][0]) & (data.index  <= anomalies_ev[-1][-1])]['T'].idxmin()
        
    event_df = event_df.append({'WUI_start': event[0], 'WUI_end': event[-1],\
                                'T_start': anomalies_ev[0][0], 'T_end': anomalies_ev[-1][-1],\
                                'nb_T_anomalies': int(len(anomalies_ev)),\
                                'is_stratified': is_strat, \
                                'Tmax': Tmax, 'Tmin': Tmin},\
                                ignore_index = True) # Merge the anomalies. If several could be interesting to save them
        
        
event_df = event_df[['WUI_start', 'WUI_end', 'T_start', 'T_end', 'Tmax', 'Tmin',\
                     'nb_T_anomalies', 'is_stratified']]
event_df = event_df.sort_values(by=['WUI_start'])
event_df['is_stratified'] = event_df['is_stratified'].astype(bool)
event_df = event_df.reset_index(drop = True) 
event_df.to_csv(join('..','Results','event_T_anomalies.csv'), index = False)     


#=======================================
# Select only suited events
#=======================================
event_df = pd.read_csv('../Results/event_T_anomalies.csv', parse_dates=['WUI_start', 'WUI_end',\
                                                  'T_start', 'T_end', 'Tmax', 'Tmin'])

# Identify events where the ocean was not in a "normal" state before the event
event_df['was_normal_start'] = pd.NaT
event_df['back_to_normal_ok'] = pd.NaT

for event_idx, event in event_df.iterrows():
    
    # Look for the end of the previous event
    prev_event_end = event_df.iloc[event_idx - 1]['T_end']\
            if event_idx > 0 else pd.to_datetime('2019-09-17 14:00:00')
       
    # Look for the start of the next event
    next_event_start = event_df.iloc[event_idx + 1]['Tmax']\
            if event_idx + 1 < len(event_df) else pd.to_datetime('2021-11-01')
    
    event_df.loc[event_idx, 'was_normal_start'] = event['Tmax'] - prev_event_end >= pd.Timedelta('1d')
    event_df.loc[event_idx, 'back_to_normal_ok'] = next_event_start - event['T_end'] >= pd.Timedelta('1d')
    
# Keep well suited events:
event_df = event_df[event_df['is_stratified'].astype(bool) &\
                        (event_df['T_end'] - event_df['T_start'] > pd.Timedelta('1d')) &\
                        event_df['back_to_normal_ok'] &\
                        event_df['was_normal_start']].reset_index(drop = True)
    
event_df['window_start'] = event_df['Tmax']  - pd.Timedelta('1d')
event_df['window_end'] = event_df['T_end']  + pd.Timedelta('1d')

# Delete the columns used to filter: lighter dataset
event_df = event_df[['WUI_start', 'WUI_end', 'T_start', 'T_end', 'Tmax', 'Tmin',
        'window_start', 'window_end']]

event_df.to_csv(join('..','Results','event_T_anomalies_suited.csv'), index = False)     


#======================================================
# Rupture detection for abundances, biovolume, biomass
#=======================================================

# Number of breakpoints to look for rupture detection
n_bkps = 2  
pfg_entities = ['abundance', 'biomass']

event_df = pd.read_csv(join('..','Results','event_T_anomalies_suited.csv'),\
                       parse_dates=['WUI_start', 'WUI_end',\
                       'T_start', 'T_end', 'Tmax', 'Tmin',\
                       'window_start', 'window_end'])
    
pfg_colors = dict(zip(['OraNano', 'OraPicoProk', 'RedNano', 'RedPico', 'RedPicoProk'],\
                      ['#9467bd', '#1f77b4', '#2ca02c', '#d62728', 'orange']))


for pfg_entity in pfg_entities:
    
    #====================================
    # Data importation 
    #====================================

    phyto = pd.read_csv(join('INSITU','PFG_' + pfg_entity + '.csv'),\
                        parse_dates = ['date'], index_col = 'date')    
    pfg_list = phyto.columns
    
    event_df['is_enough_phyto_data'] = True   
    # Create the columns to store the dates of the pfg reactions
    for pfg in pfg_list:
        event_df[pfg + '_start'] = pd.NaT
        event_df[pfg + '_end'] = pd.NaT
        event_df[pfg + '_phyto_react_before_T'] = np.nan
        event_df[pfg + '_median_before'] = np.nan
        event_df[pfg + '_median_during'] = np.nan
        event_df[pfg + '_median_after'] = np.nan
         
    #====================================
    # Iterate through the events
    #====================================
    
    for event_idx, event in event_df.iterrows():
        event_data = phyto.loc[event['window_start']:event['window_end']]
    
        # If one phase do not have enough data, then it is not considered
        if (pd.isna(event_data).mean() >= 0.3).all():
            event_df.loc[event_idx, 'is_enough_phyto_data'] = False
            continue
        
        if len(event_data) <= 2:
            print('No data for event:', event['window_start'], event['window_end'])
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
                print(event['window_start'], ': no data for', pfg)
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
            assert len(result) == n_bkps + 1
            event_df.loc[event_idx, pfg + '_start'] = result_index[0]
            event_df.loc[event_idx, pfg + '_end'] = result_index[1]
            
            # Compute the mean of the pfg before and after
            median_before = pfg_data.loc[pfg_data.index < result_index[0], pfg].median()
            median_during = pfg_data.loc[(pfg_data.index >= result_index[0]) \
                                   & (pfg_data.index < result_index[1]), pfg].median()
            median_after = pfg_data.loc[pfg_data.index >= result_index[1], pfg].median()
    
            # Store them
            event_df.loc[event_idx, pfg + '_median_before'] = median_before
            event_df.loc[event_idx, pfg + '_median_during'] = median_during
            event_df.loc[event_idx, pfg + '_median_after'] = median_after
            
            # Check that the change does not occur before the event
            if result_index[0] < event['Tmax']: #+ pd.Timedelta('2H') * 3:
                event_df.loc[event_idx, pfg + '_phyto_react_before_T'] = True
            else:
                event_df.loc[event_idx, pfg + '_phyto_react_before_T'] = False
        
    #========================================
    # Delete events which have no phyto data 
    #========================================

    event_df = event_df[event_df['is_enough_phyto_data']].reset_index(drop = True)
    event_df.to_csv(join('..','Results','responses', pfg_entity + '.csv'), index = False)     
    