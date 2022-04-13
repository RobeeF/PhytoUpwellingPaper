# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 10:04:30 2021

@author: chloe
"""

import os
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
from os.path import join
from copy import deepcopy
import matplotlib.pyplot as plt
from scipy.stats import kruskal

# Change with your path:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper/')

sn = pd.read_csv(join('Data','INSITU','sn.csv'),\
                 sep = ',', parse_dates = ['date'], index_col=['date'])
sn.columns = ['nitrites', 'phosphates', 'nitrates', 'Ammonium']

# Import the events temporal boundaries
event_df = pd.read_csv(join('Results','event_T_anomalies_suited.csv'),\
       parse_dates=['WUI_start', 'WUI_end', 'T_start',\
       'T_end', 'Tmax', 'Tmin', 'window_start', 'window_end']) 
   
# Import the stratification dates
with open(join('Results', 'date_pickles', 'stratify.pickle'), 'rb') as f1:
     stratif_dict = pickle.load(f1)
     
strat = stratif_dict['Stratified']
unstrat = stratif_dict['Unstratified']

stps = pd.read_csv(join('Data', 'INSITU','stps_1H.csv'),\
                   parse_dates = ['date'], index_col = 'date')

# Import Meteorological data
wui = pd.read_csv(join('Data', 'WUI', 'WUI_1H.csv'), parse_dates = ['date'], index_col = 'date')
wui.columns = ['WUI_calanque', 'WUI_marseille', 'WUI_cotebleue']
wui_index = 'WUI_marseille'

sn['N/P'] = sn[['nitrites', 'nitrates', 'Ammonium']].sum(1) / sn['phosphates']

#=====================================================================
# Plot the SN and color the WUI events
#=====================================================================

nb_subplots = sn.shape[1]
fig, axs = plt.subplots(nb_subplots, 1, figsize=(12, 5), sharex = True)

for s_idx, s in enumerate(sn.columns):
    axs[s_idx].scatter(sn.index, sn[s], s= 3, label = s)
    axs[s_idx].set_ylabel(s + '\n ($\mu{M}$)')

    axs[-1].set_ylabel(s)
    axs[-1].axhline(16, linestyle = "--", color = 'red') 

# Display the selected upwelling events
T_changes = []
for event_idx, event in event_df.iterrows():
                
    for ax_nb in range(nb_subplots):   
        # Plot WUI ranges
        color = 'dodgerblue' #if event['is_expected'] else 'red'
        axs[ax_nb].axvspan(event['WUI_start'], event['WUI_end'],\
                           facecolor= color, alpha=0.3)

fig.tight_layout()
plt.savefig(join('Figures','all_nutrients.png'))
plt.show()

#=========================================================================
# Get the indices of the sn in stratified wo up, with up and unstratified
#=========================================================================

#*****************
# Get the indices
#*****************
begin_col = 'WUI_start'
end_col = 'T_end'

nb_sn_measures = len(sn)
unstrat_indices = pd.DataFrame(index = sn.index,\
                             data = np.any([(sn.index >= s[0]) & (sn.index < s[1])\
                                            for s in unstrat], 0), columns = ['valid'])

strat_indices = ~unstrat_indices    
stratNonUp_indices = deepcopy(strat_indices)
stratUp_indices = pd.DataFrame(index = sn.index, data = False, columns = ['valid'])
    
for idx, event in event_df.iterrows():
    stratNonUp_indices.loc[:,'valid'] = np.where(pd.Series(stratNonUp_indices.index).between(event['window_start'],\
                                           event['window_end']),\
                                           np.full(nb_sn_measures, False), stratNonUp_indices['valid'])
        
    stratUp_indices.loc[:,'valid'] = np.where(pd.Series(stratUp_indices.index).between(event[begin_col],\
                                           event[end_col]),\
                                           np.full(nb_sn_measures, True), stratUp_indices['valid'])
 
#*****************
# Get the values
#*****************

unstr = sn.loc[unstrat_indices['valid']].replace([np.inf, -np.inf], np.nan).dropna()
str_ = sn.loc[strat_indices['valid']].replace([np.inf, -np.inf], np.nan).dropna()
str_up = sn.loc[stratUp_indices['valid']].replace([np.inf, -np.inf], np.nan).dropna()
str_nonup = sn.loc[stratNonUp_indices['valid']].replace([np.inf, -np.inf], np.nan).dropna()

#pd.set_option('display.float_format', '{:.3g}'.format)
median = pd.concat([unstr.median(), str_up.median(), str_nonup.median()], axis = 1).round(2)
q1 =  pd.concat([unstr.quantile(0.25), str_up.quantile(0.25), str_nonup.quantile(0.25)], axis = 1).round(2)
q3 =  pd.concat([unstr.quantile(0.75), str_up.quantile(0.75), str_nonup.quantile(0.75)], axis = 1).round(2)
iqr = pd.DataFrame(q3.values - q1.values, index = median.index).round(2)

stats = (median.astype(str) + " (" + iqr.astype(str) + ')')
stats.columns = ['Unstratified', 'Stratified Upwelling', 'Stratified NonUpwelling']

stats.to_csv(join('Results','sn_stratification.csv'))

#************************
# Test differences in N:P
#************************

for s in sn.columns:
    small_sample = str_nonup[s] if str_nonup[s].mean() <= str_up[s].mean() else str_up[s]
    big_sample = str_nonup[s] if str_nonup[s].mean() > str_up[s].mean()  else str_up[s]

    print(s, kruskal(small_sample, big_sample).pvalue)
    print('-----------------------------------')

#=====================================================================
# Plot the SN and the WUI for all events
#=====================================================================

n_subplots = 4
colors = sns.color_palette("Set3").as_hex()

for idx, event in event_df.iterrows():
    
    start = event['window_start'] - pd.Timedelta('3D')
    end = event['window_end'] + pd.Timedelta('3D')
    
    # Test if there are some nutrients information
    sn_data = sn[(sn.index >= start)\
                & (sn.index <= end)]
        
    fig, axs = plt.subplots(4, 1, figsize=(7, 5), sharex=True, dpi = 300)
    
    axs[0].plot(wui[(wui.index >= start)\
                 & (wui.index <= end)][wui_index],\
                label = wui_index)
    axs[0].set_ylabel('(m3/s/m)')
    axs[0].tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False)
    
    axs[0].axvspan(start, event['Tmax'], facecolor=colors[0], alpha=0.2)
    axs[0].axvspan(event['Tmax'], event['T_end'], facecolor=colors[2], alpha=0.2)
    axs[0].axvspan(event['T_end'], end, facecolor=colors[3], alpha=0.2)
    
    cols = ['red', 'green']
    for col_idx, nut in enumerate(['nitrites', 'phosphates']):
            
        axs[1].plot(sn_data[nut],\
                    label = nut, marker = '+', linestyle = 'dashed',\
                    color = cols[col_idx])
            
    axs[1].legend(loc = 'upper right', fontsize = 5)
    axs[1].grid(True)
    axs[1].tick_params(axis='x', rotation = 30, labelsize = 10) 

    for col_idx, nut in enumerate(['nitrates', 'Ammonium']):
        axs[2].plot(sn_data[nut],\
                    label = nut, marker = '+', linestyle = 'dashed')

    axs[2].legend(loc = 'upper right', fontsize = 5)

    axs[3].plot(sn_data['N/P'],\
                label = 'N/P', marker = '+', linestyle = 'dashed',\
                color = 'purple')

    axs[3].legend(loc = 'upper right', fontsize = 5)
    
    for axis_nb in range(n_subplots):
        axs[axis_nb].axvspan(start, event['Tmax'], facecolor=colors[0], alpha=0.2)
        axs[axis_nb].axvspan(event['Tmax'], event['T_end'], facecolor=colors[2], alpha=0.2)
        axs[axis_nb].axvspan(event['T_end'], end, facecolor=colors[3], alpha=0.2)

    # Axis legend handling
    for ax_nb in range(n_subplots -1):
        axs[ax_nb].tick_params(
                    axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom=False,      # ticks along the bottom edge are off
                    top=False,         # ticks along the top edge are off
                    labelbottom=False) 

    for ax_nb in range(1, n_subplots - 1):
        axs[ax_nb].set_ylabel('$\mu$M')
        axs[ax_nb].grid(True)


    axs[0].set_title('WUI and nutrient variations')
    axs[-1].set_xlabel('Date')
    axs[-1].tick_params(axis='x', rotation = 30, labelsize = 10)

    fig.tight_layout()
    fig_root = join('Figures', 'individual_events', 'SN')
    plt.savefig(join(fig_root, str(event[0])[:10] + '_' + str(event[1])[:10] + '.png'))
    plt.show()
