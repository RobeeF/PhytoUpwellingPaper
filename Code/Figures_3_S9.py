# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:30:50 2022

@author: rfuchs
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
from os.path import join
from copy import deepcopy
from matplotlib import pyplot as plt

# Change with your path:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper/')

pfg_colors = dict(zip(['Orgnano', 'Orgpicopro', 'Rednano', 'Redpicoeuk', 'Redpicopro'],\
                      ['#9467bd', '#1f77b4', '#2ca02c', '#d62728', 'orange']))

dates_col = ['WUI_start', 'WUI_end', 'T_start', 'T_end', 'Tmax', 'Tmin',
       'window_start', 'window_end', 'ORGNANO_start',
       'ORGNANO_end', 'ORGPICOPRO_start',
       'ORGPICOPRO_end', 'REDNANO_start', 'REDNANO_end',
       'REDPICOEUK_start',
       'REDPICOEUK_end', 'REDPICOPRO_start', 'REDPICOPRO_end']

ylims = {'reaction_delay': [0,9], 'reaction_duration': [0,9], 'reaction_magnitude': [-100,430],\
         'relaxation_magnitude': [-100,130]}

pfg_entities = ['abundance', 'biomass']

responses = {}
for pfg_entity in pfg_entities:

    #***************************************
    # Data importation 
    #***************************************
    path = join('Results', 'responses', pfg_entity + '.csv')
    event_df = pd.read_csv(path, parse_dates = dates_col)     
    pfg_list = [pfg.upper() for pfg in pfg_colors.keys()]
    
    #***************************************
    # Store the quantities of interest
    #***************************************
    
    response = pd.DataFrame(columns = ['reaction_delay', 'reaction_duration', 'reaction_magnitude',\
                                                   'relaxation_magnitude', 'pfg', 'WUI_start']) 
    for pfg in pfg_list:
        # Trigger length
        reaction_delay =  (event_df[pfg + '_start'] - event_df['Tmax']).apply(lambda x: x.total_seconds() / (3600 * 24))
        reaction_delay = np.where(event_df[pfg + '_phyto_react_before_T'], np.nan, reaction_delay)
    
        # Trigger time
        reaction_duration =  (event_df[pfg + '_end'] - event_df[pfg + '_start']).apply(lambda x: x.total_seconds() / (3600 * 24))
            
        # Reaction_strength
        reac_magn =  ((event_df[pfg + 'median_during'] / event_df[pfg + 'median_before']) - 1) * 100
        relax_magn =  ((event_df[pfg + 'median_after'] / event_df[pfg + 'median_during']) - 1) * 100
    
        df = pd.DataFrame(data = np.stack([reaction_delay, reaction_duration, reac_magn, relax_magn]).T,\
                          columns = ['reaction_delay', 'reaction_duration', 'reaction_magnitude', 'relaxation_magnitude'])
        df.loc[event_df[pfg + '_phyto_react_before_T']] = np.nan
        
        df['pfg'] = pfg
        df['WUI_start'] = event_df['WUI_start']
        
        response = response.append(df)
            
        
    response = response[~response.isna().any(1)]
    response['pfg'] = response['pfg'].apply(lambda x: str.capitalize(x))

    responses[pfg_entity] = deepcopy(response)
     
    
#========================================
# Draw the boxplots
#========================================

for pfg_entity in pfg_entities:
    response = responses[pfg_entity]
    
    for entity_tracked in ['reaction_delay', 'reaction_duration', 'reaction_magnitude', 'relaxation_magnitude']:
        
        # Add units and lines
        if entity_tracked in ['reaction_delay', 'reaction_duration']:
            unit = 'days'
        else:
            unit = '% variation'
            
        ax = sns.boxplot(x = "pfg", y = entity_tracked,
                         data=response, \
                         showfliers = False, palette = pfg_colors)
            
        # Add the number of points on which each boxplot is based upon
        medians = response.groupby(['pfg'])[entity_tracked].median().values
        nobs = response['pfg'].value_counts().values
        nobs = [str(x) for x in nobs.tolist()]
        nobs = ["n: " + i for i in nobs]
         
        # Add it to the plot
        pos = range(len(nobs))
        for tick,label in zip(pos,ax.get_xticklabels()):
            ax.text(pos[tick],
                    medians[tick] + np.abs(medians[tick]) * 0.08,
                    nobs[tick],
                    horizontalalignment='center',
                    size='x-small',
                    color='black',
                    weight='semibold')
            
        # Display settings
        ax.tick_params('x', rotation = 45)
        ax.set_ylabel(entity_tracked + ' (' + unit + ')') 
        ax.set_title('PFG ' + entity_tracked + ' (' + pfg_entity + ')' + ' during an Upwelling')
        ax.set_ylim(ylims[entity_tracked])
        
        plt.tight_layout()
        fig_path = join('Figures', 'boxplots', pfg_entity, entity_tracked + '.png')
        plt.savefig(fig_path)
        plt.show()
    
#==============================================
# Focusing on Reaction vs relaxation magnitude
#==============================================

fig, axs = plt.subplots(5, 2, figsize = (15, 25))
for entity_idx, pfg_entity in enumerate(pfg_entities):
    response = responses[pfg_entity].set_index('pfg')
    
    for pfg_idx, pfg in enumerate(pfg_colors.keys()):
        axs[pfg_idx, entity_idx].scatter(response.loc[pfg, 'reaction_magnitude'],\
                                    response.loc[pfg, 'relaxation_magnitude'],\
                                    c = pfg_colors[pfg], label = pfg)

        axs[pfg_idx, entity_idx].set_title(pfg)
        axs[pfg_idx, entity_idx].set_xlabel('Reaction magnitude (%)')
        axs[pfg_idx, entity_idx].set_ylabel('Relaxation magnitude (%)')
        
plt.tight_layout()
fig_path = join('Figures', 'Reaction_relaxation.png')
plt.savefig(fig_path)
plt.show()
    
