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

pfg_colors = dict(zip(['OraNano', 'OraPicoProk', 'RedNano', 'RedPico', 'RedPicoProk'],\
                      ['#9467bd', '#1f77b4', '#2ca02c', '#d62728', 'orange']))

dates_col = ['WUI_start', 'WUI_end', 'T_start', 'T_end', 'Tmax', 'Tmin',
       'window_start', 'window_end', 'OraNano_start',
       'OraNano_end', 'OraPicoProk_start',
       'OraPicoProk_end', 'RedNano_start', 'RedNano_end',
       'RedPico_start',
       'RedPico_end', 'RedPicoProk_start', 'RedPicoProk_end']

ylims = {'reaction_delay': [0,9], 'reaction_duration': [0,9], 'reaction_magnitude': [-100,430],\
         'relaxation_magnitude': [-100,130]}

pfg_entities = ['abundance', 'biomass']
pfg_list = [pfg for pfg in pfg_colors.keys()]

responses = {}
for pfg_entity in pfg_entities:

    #***************************************
    # Data importation 
    #***************************************
    path = join('Results', 'responses', pfg_entity + '.csv')
    event_df = pd.read_csv(path, parse_dates = dates_col)     
    
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
        reac_magn =  ((event_df[pfg + '_median_during'] / event_df[pfg + '_median_before']) - 1) * 100
        relax_magn =  ((event_df[pfg + '_median_after'] / event_df[pfg + '_median_during']) - 1) * 100
    
        df = pd.DataFrame(data = np.stack([reaction_delay, reaction_duration, reac_magn, relax_magn]).T,\
                          columns = ['reaction_delay', 'reaction_duration', 'reaction_magnitude', 'relaxation_magnitude'])
        df.loc[event_df[pfg + '_phyto_react_before_T']] = np.nan
        
        df['pfg'] = pfg
        df['WUI_start'] = event_df['WUI_start']
        
        response = pd.concat([response, df])
            
    response = response[~response.isna().any(1)]
    response['pfg'] = response['pfg']#.apply(lambda x: str.capitalize(x))

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
    

#========================================
# Summary table
#========================================    

path = join('Results', 'responses', 'abundance' + '.csv')
abundance_df = pd.read_csv(path, parse_dates = dates_col, index_col=['window_start', 'window_end'])   
path = join('Results', 'responses', 'biomass' + '.csv')
biomass_df = pd.read_csv(path, parse_dates = dates_col, index_col=['window_start', 'window_end'])   

# Set unindentified ruptures to Nan
for pfg in pfg_list:
    col_to_nans = [pfg + '_median_before', pfg + '_median_during', pfg + '_median_after']
    abundance_df.loc[abundance_df[pfg + '_phyto_react_before_T'], col_to_nans] = np.nan
    biomass_df.loc[biomass_df[pfg + '_phyto_react_before_T'], col_to_nans] = np.nan

cols_to_keep = ['OraNano_median_before',
'OraNano_median_during', 'OraNano_median_after', 
'OraPicoProk_median_before', 'OraPicoProk_median_during', 
'OraPicoProk_median_after',  'RedNano_median_before', 
'RedNano_median_during', 'RedNano_median_after', 
'RedPico_median_before', 'RedPico_median_during', 
'RedPico_median_after', 'RedPicoProk_median_before', 
'RedPicoProk_median_during', 'RedPicoProk_median_after']

cols_aliases = ['OraNano pre-reaction',
'OraNano reaction', 'OraNano relaxation', 
'OraPicoProk pre-reaction', 'OraPicoProk reaction', 
'OraPicoProk relaxation',  'RedNano pre-reaction', 
'RedNano reaction', 'RedNano relaxation', 
'RedPico pre-reaction', 'RedPico reaction', 
'RedPico relaxation', 'RedPicoProk pre-reaction', 
'RedPicoProk reaction', 'RedPicoProk relaxation']

# Round and convert to the right unit
abundance_df = (abundance_df[cols_to_keep] * 10 ** 3).round(2) 
abundance_df.columns = cols_aliases
abundance_df.to_csv(join('Results', 'abundances_all_events.csv'))

# !!! Scientific notation ?
biomass_df = biomass_df[cols_to_keep] * 10 ** -9
biomass_df.columns = cols_aliases
biomass_df.round(8).to_csv(join('Results', 'biomass_all_events.csv'))

