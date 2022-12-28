# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 16:53:44 2021

@author: rfuchs
"""

import os
import pandas as pd
import seaborn as sns
from os.path import join
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Change with your path:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper/')

pfg_entities = ['biomass', 'abundance', 'mu']

#==============================================================================
# Data Importation
#==============================================================================

units = {'loss': '($h^{-1}$)', 'mu': '($h^{-1}$)', 'biomass': '($mgC.m{L}^{-1}$)',\
         'abundance': '($cells.m{L}^{-1}$)', 'biovolume': '($\mu^3$)',\
         'NPP': '($mgC.m^{-3}h^{-1}$)'}
        
conv = {entity: 1 for entity in list(units.keys())}
conv['abundance'] = 10 ** 3
conv['biomass'] = 10 ** -9

daily_estimation = False 
if daily_estimation:
    units['loss'] = '($d^{-1}$)'
    units['mu'] = '($d^{-1}$)'
    units['NPP'] = '($mgC.m^{-3}d^{-1}$)'
    
# Create the legends at the bottom of the image
c = ['#9467bd', '#1f77b4', '#2ca02c', '#d62728', 'orange']
ccs = [Line2D([None],[None], color = c[i]) for i in range(5)]

# construct cmap
pfg_colors = dict(zip(['OraNano', 'OraPicoProk', 'RedNano', 'RedPico',\
                       'RedPicoProk'], c))
colors = sns.color_palette("Set3").as_hex()

for entity_tracked in pfg_entities:
    root = join('Data','INSITU') if entity_tracked in ['biomass', 'abundance'] else 'Results'
    path = join(root, 'PFG_' + entity_tracked + '.csv')
    phyto = pd.read_csv(path, parse_dates = ['date']) 
            
    # Import phytoplantkon data
    if entity_tracked == 'biomass':
        freq = '2H'
    elif entity_tracked == 'abundance':
        freq = '2H'
    
    elif entity_tracked == 'biovolume':
        freq = '2H'
        
    elif entity_tracked in ['mu', 'loss', 'NPP']:
                        
        if daily_estimation:
            freq = '24H'
            phyto = phyto.groupby(pd.Grouper(key = 'date', freq = freq)).sum() * (24 / 25)
            phyto = phyto.reset_index()
        else:
            freq = '1H'
    else:
        raise RuntimeError('Please enter a legal entity to track')
    
    phyto.set_index('date', inplace = True)
            
    # Import Meteorological data
    wui = pd.read_csv('Data/WUI/WUI_1H.csv', parse_dates = ['date'])
    wui.set_index('date', inplace = True)
    wui.columns = ['WUI_calanque', 'WUI_marseille', 'WUI_cotebleue']
    
    # Format the dataset
    phyto_wui = phyto.join(wui)
    
    stps = pd.read_csv('Data/INSITU/stps_' + freq + '.csv', parse_dates=['date'],\
                       index_col=['date'])
    
    # Import the rupture points
    event_df = pd.read_csv(join('Results', 'responses', entity_tracked + '.csv'),\
                            parse_dates=['WUI_start', 'WUI_end', 'T_start',\
                                         'T_end', 'Tmax', 'Tmin']) 

    #==============================================================================
    # Plotting utility: PFGs vs Temp. vs WUI
    #==============================================================================
            
    # Import also the T data
    phyto_wui_T = phyto_wui.join(stps, how="inner")        
    phyto_cols = phyto.columns
                 
    nb_subplots = 3
    
    # Plotting options depending on the quantity to plot
    # Plotting options depending on the quantity to plot
    if entity_tracked == 'biomass':
        axis1 = ["RedPicoProk"]
        axis1_bis = ['OraPicoProk']
        axis2 = ['RedPico']
        axis2_bis = ["RedNano", "OraNano"]
    elif entity_tracked == 'biovolume':
        axis1 = ["RedPicoProk", "OraPicoProk"]
        axis1_bis = ['RedPico']
        axis2 = ["RedNano"]
        axis2_bis = ["OraNano"]
    elif entity_tracked == 'abundance':
        axis1 = ["RedPicoProk"]
        axis1_bis = ["OraPicoProk"]
        axis2 = ["RedPico", "RedNano"] 
        axis2_bis = ["OraNano"]
    elif entity_tracked in ['mu', 'NPP']:
        axis1 = ['RedPico']
        axis1_bis = ['OraPicoProk', 'RedPico']
        axis2 = ['RedNano']  
        axis2_bis = ["OraNano"]
    else:
        raise ValueError('Please enter a legal entity to track')
    
    
    for event_idx, event in event_df.iterrows():        
        event_data = phyto_wui_T.loc[(phyto_wui_T.index >= event['window_start']) &\
                                     (phyto_wui_T.index < event['window_end'])]
            
        if len(event_data) <= 2:
            print('event', event['window_start'], 'to', event['window_end'], 'has no data')
            continue
        
        #=================================
        # First axis
        #=================================
        fig, axs = plt.subplots(3, 1, figsize=(7, 5), sharex=True, dpi = 300)
        axs[0].grid(True)

        axs[0].plot(event_data['WUI_marseille'])
        axs[0].set_xlabel('Date')
        axs[0].tick_params(axis='y', labelcolor = 'tab:blue')
        axs[0].set_ylabel('WUDI (m3/s/m)', color = 'tab:blue')
        
        ax0_bis = axs[0].twinx()  # instantiate a second axes that shares the same x-axis
    
        color = 'tab:orange'
        ax0_bis.set_ylabel('Temperature (°C)', color=color)  # we already handled the x-label with ax1
        ax0_bis.plot(event_data['T'], color=color)
        ax0_bis.tick_params(axis='y', labelcolor=color)
        
        #=================================
        # Second axis
        #=================================

        axs[1].grid(True)
        for idx, pfg in enumerate(axis1):
            axs[1].plot(event_data[pfg] * conv[entity_tracked], label = pfg, color = pfg_colors[pfg])
            
            if not(pd.isnull(event[pfg + '_start']) or pd.isnull(event[pfg + '_end'])):
                axs[1].axvline(pd.to_datetime(event[pfg + '_start']),\
                               color = pfg_colors[pfg], linestyle='dashed')
                axs[1].axvline(pd.to_datetime(event[pfg + '_end']),\
                               color = pfg_colors[pfg], linestyle='dashed')
                          
        axs[1].set_ylabel(units[entity_tracked], color = pfg_colors[axis1[-1]])
        axs[1].tick_params(axis='y', labelcolor =  pfg_colors[axis1[0]])
        
        ax1_bis = axs[1].twinx()  # instantiate a second axes that shares the same x-axis

        for idx, pfg in enumerate(axis1_bis):
            ax1_bis.plot(event_data[pfg] * conv[entity_tracked], label = pfg,\
                             color = pfg_colors[pfg]) 
                    
            if not(pd.isnull(event[pfg + '_start']) or pd.isnull(event[pfg + '_end'])):
                ax1_bis.axvline(pd.to_datetime(event[pfg + '_start']),\
                               color = pfg_colors[pfg], linestyle='dashed')
                ax1_bis.axvline(pd.to_datetime(event[pfg + '_end']),\
                               color = pfg_colors[pfg], linestyle='dashed')

        ax1_bis.tick_params(axis='y')
        ax1_bis.tick_params(axis='y', labelcolor = pfg_colors[axis1_bis[0]])
        ax1_bis.set_ylabel(units[entity_tracked], color = pfg_colors[axis1_bis[-1]])

        #=================================
        # Third axis
        #=================================
        
        for idx, pfg in enumerate(axis2):
            axs[2].plot(event_data[pfg] * conv[entity_tracked], label = pfg, color = pfg_colors[pfg])
            
            if not(pd.isnull(event[pfg + '_start']) or pd.isnull(event[pfg + '_end'])):
                axs[2].axvline(pd.to_datetime(event[pfg + '_start']),\
                               color = pfg_colors[pfg], linestyle='dashed')
                axs[2].axvline(pd.to_datetime(event[pfg + '_end']),\
                               color = pfg_colors[pfg], linestyle='dashed')
        
        axs[2].grid(True)
        axs[2].tick_params(axis='x', rotation = 30, labelsize = 10)
        axs[2].tick_params(axis='y', labelcolor = pfg_colors[axis2[0]])
        axs[2].set_ylabel(units[entity_tracked], color = pfg_colors[axis2[-1]])

        ax2_bis = axs[2].twinx()  # instantiate a second axes that shares the same x-axis
            
        for idx, pfg in enumerate(axis2_bis):
            ax2_bis.plot(event_data[pfg] * conv[entity_tracked], label = pfg, color = pfg_colors[pfg])
            
            if not(pd.isnull(event[pfg + '_start']) or pd.isnull(event[pfg + '_end'])):
                ax2_bis.axvline(pd.to_datetime(event[pfg + '_start']),\
                               color = pfg_colors[pfg], linestyle='dashed')
                ax2_bis.axvline(pd.to_datetime(event[pfg + '_end']),\
                               color = pfg_colors[pfg], linestyle='dashed')
                
        ax2_bis.tick_params(axis='y', labelcolor = pfg_colors[axis2_bis[0]])
        ax2_bis.set_ylabel(units[entity_tracked], color = pfg_colors[axis2_bis[-1]])
        
        #=================================
        # Final axes formatting
        #=================================
        
        # Display the color spans for the physical phases
        for ax_nb in range(nb_subplots):
            axs[ax_nb].axvspan(event['window_start'], event['Tmax'], facecolor=colors[0], alpha=0.3)
            axs[ax_nb].axvspan(event['Tmax'], event['T_end'], facecolor=colors[2], alpha=0.3)
            axs[ax_nb].axvspan(event['T_end'], event['window_end'], facecolor=colors[3], alpha=0.3)
        
        # Axis legend handling
        for ax_nb in range(nb_subplots - 1):
            axs[ax_nb].tick_params(
                        axis='x',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        bottom=False,      # ticks along the bottom edge are off
                        top=False,         # ticks along the top edge are off
                        labelbottom=False,
                        rotation = 45)   
            
        fig.legend(ccs, [c for c in phyto.columns],\
                   loc="lower center", bbox_to_anchor=(0.5, -0.05),\
                   ncol = 5, fontsize = 9)
        fig.tight_layout()
                  
        fig_path = os.path.join('Figures', 'individual_events',\
                                entity_tracked, str(event['WUI_start'].date())\
                                + '.png')
        plt.savefig(fig_path)
        plt.show()