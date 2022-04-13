# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 14:30:12 2021

@author: rfuchs
"""

import os
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
from os.path import join
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Change with your path:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper/')

# Choose a quantity to track: biomass or abundance, mu
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
conv['abundance'] = 10 ** -3
conv['biomass'] = 10 ** -9

phyto = pd.read_csv(join('Data', 'INSITU', 'PFG_' + entity_tracked + '.csv'),\
                    parse_dates = ['date'], index_col = 'date') 

# Import Meteorological data
wui = pd.read_csv(join('Data', 'WUI', 'WUI_2H.csv'), parse_dates = ['date'], index_col = 'date')
wui.columns = ['WUI_calanque', 'WUI_marseille', 'WUI_cotebleue']

# Format the dataset
phyto_wui = phyto.join(wui)

stps = pd.read_csv(join('Data', 'INSITU', 'stps_' + freq[entity_tracked] + '.csv'),\
                   parse_dates = ['date'], index_col='date')

phyto = phyto[phyto.index < pd.to_datetime('2021-11-01')]
phyto_wui = phyto_wui[phyto_wui.index < pd.to_datetime('2021-11-01')]

# Import the events temporal boundaries
event_df = pd.read_csv(join('Results', 'event_T_anomalies_suited.csv'),\
       parse_dates=['WUI_start', 'WUI_end', 'T_start',\
       'T_end', 'Tmax', 'Tmin', 'window_start', 'window_end']) 
    
# Import the stratification dates
with open(join('Results', 'date_pickles', 'stratify.pickle'), 'rb') as f1:
     stratif_dict = pickle.load(f1)
strat = stratif_dict['Stratified']

# Focus on one event
eventOfInterest = event_df.loc[0]

blooms = pd.read_csv(join('Results','blooms', entity_tracked  + '50.csv'),\
                     index_col = 'PFG') 

#==============================================================================
# Plotting utility: PFGs vs Temp. vs WUI
#==============================================================================

# Import also the T data

phyto_wui_T = phyto_wui.join(stps, how="inner")
phyto_wui_T.index = pd.to_datetime(phyto_wui_T.index)

phyto_cols = phyto.columns
             
# construct cmap
c = ['#9467bd', '#1f77b4', '#2ca02c', '#d62728', 'orange']
pfg_colors = dict(zip(phyto.columns, c))
colors = sns.color_palette("Set3").as_hex()

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

fig, axs = plt.subplots(3, 3, figsize=(12, 5), gridspec_kw = {'width_ratios':[1, 2, 2]},\
                        dpi = 300)
nb_subplots = len(axs)
years = [2019, 2020, 2021]

#=================================
# First axis
#=================================


# Store the max and min of original and bis y axes 
max_ax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]]['WUI_marseille'].max()\
          for i, year in enumerate(years)]
min_ax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]]['WUI_marseille'].min()\
          for i, year in enumerate(years)]

max_bisax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]]['T'].max()\
             for i, year in enumerate(years)]
min_bisax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]]['T'].min()\
             for i, year in enumerate(years)]


for i, year in enumerate(years):
    
    year_wui_T = phyto_wui_T.loc[strat[i][0]:strat[i][1]]
    year_stps = stps.loc[strat[i][0]:strat[i][1]]
    year_wui = wui.loc[strat[i][0]:strat[i][1]]
    
    axs[0][i].grid(True)    
    axs[0][i].plot(year_wui.index, year_wui['WUI_marseille'])    
    axs[0][i].tick_params(axis='y', labelcolor = 'tab:blue')
    
    axs[0][i].set_ylim([np.min(min_ax) * 0.9, np.max(max_ax) * 1.1])
    
    if i != 0:
        axs[0][i].set_yticklabels([])
    
    ax0_bis = axs[0][i].twinx()  # instantiate a second axes that shares the same x-axis
    
    color = 'tab:orange'
    ax0_bis.plot(year_stps.index, year_stps['T'], color=color)
    ax0_bis.tick_params(axis='y', labelcolor=color)
            
    ax0_bis.set_ylim([np.min(min_bisax) * 0.9, np.max(max_bisax) * 1.1])
    
    if i != 2:
        ax0_bis.set_yticklabels([])
    else:
        ax0_bis.tick_params(axis='y', labelcolor = 'tab:orange')


# Set the labels on the corner plots
axs[0][0].set_ylabel('WUDI (m3/s/m)', color = 'tab:blue')  
ax0_bis.set_ylabel('Temperature (°C)', color=color)  # we already handled the x-label with ax1

#=================================
# Second axis
#=================================

# Store the max and min of original and bis y axes 
max_ax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]][pfg].max()\
          for i, year in enumerate(years) for pfg in axis1]
max_ax += [blooms.loc[pfg].max() for pfg in axis1]

min_ax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]][pfg].min()\
          for i, year in enumerate(years) for pfg in axis1]

max_bisax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]][pfg].max()\
          for i, year in enumerate(years) for pfg in axis1_bis]
max_bisax += [blooms.loc[pfg].max() for pfg in axis1_bis]

min_bisax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]][pfg].min()\
          for i, year in enumerate(years) for pfg in axis1_bis]

for i, year in enumerate(years):
    year_wui_T = phyto_wui_T.loc[strat[i][0]:strat[i][1]]
        
    for idx, pfg in enumerate(axis1):
        axs[1][i].plot(year_wui_T[pfg] * conv[entity_tracked], label = pfg,\
                       color = pfg_colors[pfg])#, zorder=1) 
        axs[1][i].axhline(blooms.loc[pfg, str(year)] * conv[entity_tracked],\
                        color = pfg_colors[pfg], linestyle = '--')#, zorder=2)
          
    axs[1][i].grid(True)
    axs[1][i].tick_params(axis='y', labelcolor =  pfg_colors[axis1[0]])

    # Hide ticks
    axs[1][i].set_ylim([np.min(min_ax) * 0.9 * conv[entity_tracked],\
                        np.max(max_ax) * 1.1 * conv[entity_tracked]])
    
    if i != 0:
        axs[1][i].set_yticklabels([])

    ax1_bis = axs[1][i].twinx()  # instantiate a second axes that shares the same x-axis
    
    for idx, pfg in enumerate(axis1_bis):
        ax1_bis.plot(year_wui_T[pfg] * conv[entity_tracked], label = pfg,\
                         color = pfg_colors[pfg])#, zorder=1) 
        # The unstratified max level
        ax1_bis.axhline(blooms.loc[pfg, str(year)] * conv[entity_tracked],\
                        color = pfg_colors[pfg], linestyle = '--')#, zorder=2)
                
    ax1_bis.tick_params(axis='y')

    # Hide ticks
    ax1_bis.set_ylim([np.min(min_bisax) * 0.9 * conv[entity_tracked],\
                      np.max(max_bisax) * 1.1 * conv[entity_tracked]])
    
    if i != 2:
        ax1_bis.set_yticklabels([])
    else:
        ax1_bis.tick_params(axis='y', labelcolor = pfg_colors[axis1_bis[0]])
    
# Set the labels on the corner plots
axs[1][0].set_ylabel(units[entity_tracked], color = pfg_colors[axis1[-1]])
ax1_bis.set_ylabel(units[entity_tracked], color = pfg_colors[axis1_bis[-1]])


#=================================
# Third axis
#=================================

# Store the max and min of original and bis y axes 
max_ax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]][pfg].max()\
          for i, year in enumerate(years) for pfg in axis2]
max_ax += [blooms.loc[pfg].max() for pfg in axis2]

min_ax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]][pfg].min()\
          for i, year in enumerate(years) for pfg in axis2]

max_bisax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]][pfg].max()\
          for i, year in enumerate(years) for pfg in axis2_bis]
max_bisax += [blooms.loc[pfg].max() for pfg in axis2_bis]

min_bisax = [phyto_wui_T.loc[strat[i][0]:strat[i][1]][pfg].min()\
          for i, year in enumerate(years) for pfg in axis2_bis]

for i, year in enumerate(years):
    year_wui_T = phyto_wui_T.loc[strat[i][0]:strat[i][1]]

    for idx, pfg in enumerate(axis2):
        axs[2][i].plot(year_wui_T[pfg] * conv[entity_tracked], label = pfg,\
                       color = pfg_colors[pfg])#, zorder=1)
        axs[2][i].axhline(blooms.loc[pfg, str(year)] * conv[entity_tracked],\
                        color = pfg_colors[pfg], linestyle = '--')#, zorder=2)
        
        
    axs[2][i].grid(True)
    axs[2][i].tick_params(axis='y', labelcolor = pfg_colors[axis2[0]])

    # Hide ticks
    axs[2][i].set_ylim([np.min(min_ax) * 0.9 * conv[entity_tracked],\
                        np.max(max_ax) * 1.1 * conv[entity_tracked]])
    
    if i != 0:
        axs[2][i].set_yticklabels([])

    ax2_bis = axs[2][i].twinx()  # instantiate a second axes that shares the same x-axis
        
    for idx, pfg in enumerate(axis2_bis):
        ax2_bis.plot(year_wui_T[pfg] * conv[entity_tracked], label = pfg,\
                     color = pfg_colors[pfg])#, zorder=1)
        # The unstratified max level
        ax2_bis.axhline(blooms.loc[pfg, str(year)] * conv[entity_tracked],\
                        color = pfg_colors[pfg], linestyle = '--')#, zorder=2)
        
    ax2_bis.tick_params(axis='y', labelcolor = pfg_colors[axis2_bis[0]])

    # Hide ticks
    ax2_bis.set_ylim([np.min(min_bisax) * 0.9 * conv[entity_tracked],\
                      np.max(max_bisax) * 1.1 * conv[entity_tracked]])
    
    if i != 2:
        ax2_bis.set_yticklabels([])
    else:
        ax2_bis.tick_params(axis='y', labelcolor = pfg_colors[axis2_bis[0]])


# Set the labels on the corner plots
axs[2][0].set_ylabel(units[entity_tracked], color = pfg_colors[axis2[-1]])
ax2_bis.set_ylabel(units[entity_tracked], color = pfg_colors[axis2_bis[-1]])


# Display the selected upwelling events
for i, year in enumerate(years):
    for ax_nb in range(nb_subplots):
        event_df_year = event_df[event_df['window_start'].dt.year == year]
        for event_idx, event in event_df_year.iterrows():
            axs[ax_nb][i].axvspan(event['WUI_start'], event['WUI_end'],\
                               facecolor= 'dodgerblue', alpha=0.3)#, zorder=3)
                       
# Create a Rectangle patch
for ax_nb in range(nb_subplots):
    axs[ax_nb][0].axvline(eventOfInterest['WUI_start'], color = 'darkblue')#, zorder=4)
    axs[ax_nb][0].axvline(eventOfInterest['WUI_end'], color = 'darkblue')#, zorder=4)
    
# Hide the x labels except on the last row
for i, year in enumerate(years):
    for ax_nb in range(nb_subplots - 1):
        axs[ax_nb][i].tick_params(
                    axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom=False,      # ticks along the bottom edge are off
                    top=False,         # ticks along the top edge are off
                    labelbottom=False)   

for i, year in enumerate(years):
    axs[nb_subplots - 1][i].tick_params(
                        axis='x',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        rotation = 45)  

# Create the legends at the bottom of the image
ccs = [Line2D([None],[None], color = c[i]) for i in range(5)]

fig.legend(ccs, [p.capitalize() for p in phyto.columns],\
           loc="lower center", bbox_to_anchor=(0.5, -0.05), ncol = 5)

fig.tight_layout()
plt.savefig(join('Figures', 'stratify_all_PFGS_' + entity_tracked + '.png'))
plt.show()