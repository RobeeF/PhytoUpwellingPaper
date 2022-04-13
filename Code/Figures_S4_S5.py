# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 14:30:12 2021

@author: rfuchs
"""

import os
import pandas as pd
import seaborn as sns
from os.path import join
import matplotlib.pyplot as plt

# Change with your path:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper/')

# Choose a quantity to track: biomass or abundance
entity_tracked = 'abundance'

#==============================================================================
# Data Importation
#==============================================================================
    
units = {'loss': '($d^{-1}$)', 'mu': '($d^{-1}$)', 'biomass': '($mgC.mL^{-1}$)',\
         'abundance': '($cells.m{L}^{-1}$)', 'biovolume': '($\mu^3$)',\
         'NPP': '($mgC.m^{-3}d^{-1}$)'}
    
freq = {'loss': '1H', 'mu': '1H', 'biomass': '2H',\
         'abundance': '2H', 'biovolume': '2H',\
         'NPP': '1H'}
    
conv = {entity: 1 for entity in list(units.keys())}
conv['abundance'] = 10 ** 3
conv['biomass'] = 10 ** -9


phyto = pd.read_csv(join('Data', 'INSITU', 'PFG_' + entity_tracked + '.csv'),\
                    parse_dates = ['date'], index_col = 'date') 

# Import Meteorological data

wui = pd.read_csv(join('Data', 'WUI', 'WUI_2H.csv'), parse_dates = ['date'])
wui.set_index('date', inplace = True)
wui.columns = ['WUI_calanque', 'WUI_marseille', 'WUI_cotebleue']

# Format the dataset
phyto_wui = phyto.join(wui)


stps = pd.read_csv(join('Data', 'INSITU', 'stps_' + freq[entity_tracked] + '.csv'),\
                   parse_dates = ['date'], index_col='date')

phyto = phyto[phyto.index < pd.to_datetime('2021-10-31')]
stps = stps[stps.index < pd.to_datetime('2021-10-31')]
wui = wui.loc[pd.to_datetime('2019-09-18'):pd.to_datetime('2021-10-31')]

# Import the events temporal boundaries

event_df = pd.read_csv(join('Results', 'event_T_anomalies_suited.csv'),\
       parse_dates=['WUI_start', 'WUI_end', 'T_start',\
       'T_end', 'Tmax', 'Tmin']) 

# Focus on one event
eventOfInterest = event_df.loc[0]

#==============================================================================
# Plotting utility: PFGs vs Temp. vs WUI
#==============================================================================

# Import also the T data
phyto_wui_T = phyto_wui.join(stps, how="inner")
phyto_wui_T.index = pd.to_datetime(phyto_wui_T.index)

phyto_cols = phyto.columns
             
# construct cmap
pfg_colors = dict(zip(phyto.columns, ['#9467bd', '#1f77b4', '#2ca02c', '#d62728', 'orange']))
colors = sns.color_palette("Set3").as_hex()

nb_subplots = 3

# Plotting options depending on the quantity to plot
if entity_tracked == 'biomass':
    axis1 = ["REDPICOPRO"]
    axis1_bis = ['ORGPICOPRO']
    axis2 = ["REDPICOEUK"]
    axis2_bis = ['REDNANO', "ORGNANO"]
elif entity_tracked == 'abundance':
    axis1 = ["REDPICOPRO"]
    axis1_bis = ["ORGPICOPRO"]
    axis2 = ["REDPICOEUK", "REDNANO"] 
    axis2_bis = ["ORGNANO"]
elif entity_tracked in ['mu', 'NPP']:
    axis1 = ['REDPICROPRO', 'REDPICOEUK']
    axis1_bis = ['ORGPICOPRO']
    axis2 = ['REDNANO']  
    axis2_bis = ["ORGNANO"]
else:
    raise ValueError('Please enter a legal entity to track')

pfg_fontsize = 8
fig, axs = plt.subplots(3, 1, figsize=(8, 5), sharex=True, dpi = 300)

#=================================
# First axis
#=================================

axs[0].grid(True)
axs[0].plot(wui['WUI_marseille'])
axs[0].tick_params(axis='x', rotation = 30, labelsize = 10)
axs[0].tick_params(axis='y', labelcolor = 'tab:blue')
axs[0].set_ylabel('WUDI (m3/s/m)', color = 'tab:blue')

ax0_bis = axs[0].twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:orange'
ax0_bis.set_ylabel('Temperature (°C)', color=color)  # we already handled the x-label with ax1
ax0_bis.plot(stps['T'], color=color)
ax0_bis.tick_params(axis='y', labelcolor=color)

#=================================
# First axis
#=================================

for idx, pfg in enumerate(axis1):
    axs[1].plot(phyto[pfg] * conv[entity_tracked], label = pfg, color = pfg_colors[pfg]) 

axs[1].grid(True)
axs[1].legend(loc = 'upper left', fontsize = pfg_fontsize)
axs[1].tick_params(axis='y', labelcolor =  pfg_colors[axis1[0]])
axs[1].set_ylabel(units[entity_tracked], color = pfg_colors[axis1[-1]])

ax1_bis = axs[1].twinx()  # instantiate a second axes that shares the same x-axis

for idx, pfg in enumerate(axis1_bis):
    ax1_bis.plot(phyto[pfg] * conv[entity_tracked], label = pfg,\
                     color = pfg_colors[pfg]) 
            
ax1_bis.tick_params(axis='y')
ax1_bis.legend(loc = 'upper right', fontsize = pfg_fontsize)
ax1_bis.tick_params(axis='y', labelcolor = pfg_colors[axis1_bis[0]])
ax1_bis.set_ylabel(units[entity_tracked], color = pfg_colors[axis1_bis[-1]])

#=================================
# Second axis
#=================================

for idx, pfg in enumerate(axis2):
    axs[2].plot(phyto[pfg] * conv[entity_tracked], label = pfg, color = pfg_colors[pfg])
    
axs[2].legend(loc = 'upper left', fontsize = pfg_fontsize)
axs[2].grid(True)
axs[2].tick_params(axis='y', labelcolor = pfg_colors[axis2[0]])
axs[2].set_ylabel(units[entity_tracked], color = pfg_colors[axis2[-1]])

ax2_bis = axs[2].twinx()  # instantiate a second axes that shares the same x-axis
    
for idx, pfg in enumerate(axis2_bis):
    ax2_bis.plot(phyto[pfg] * conv[entity_tracked], label = pfg, color = pfg_colors[pfg])
    
        
ax2_bis.legend(loc = 'upper right', fontsize = pfg_fontsize) 
ax2_bis.tick_params(axis='y', labelcolor = pfg_colors[axis2_bis[0]])
ax2_bis.set_ylabel(units[entity_tracked], color = pfg_colors[axis2_bis[-1]])

# Display the selected upwelling events
for ax_nb in range(nb_subplots):
    
    # Plot all the events
    for event_idx, event in event_df.iterrows():
        axs[ax_nb].axvspan(event['WUI_start'], event['WUI_end'],\
                           facecolor= 'dodgerblue', alpha=0.3)

    # Create a Rectangle patch on the illustrative event
    axs[ax_nb].axvline(eventOfInterest['WUI_start'], color = 'darkblue')    
    axs[ax_nb].axvline(eventOfInterest['WUI_end'], color = 'darkblue')
        
    
# Axis legend handling
for ax_nb in range(nb_subplots - 1):
    axs[ax_nb].tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False,
                rotation = 45)   
    

fig.tight_layout()
plt.savefig(join('Figures', 'all_PFGS_' + entity_tracked + '.png'))
plt.show()