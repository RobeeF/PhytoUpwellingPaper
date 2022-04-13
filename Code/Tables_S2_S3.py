# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 14:24:16 2022

@author: rfuchs
"""

import os
import pickle
import numpy as np
import pandas as pd
from os.path import join
from copy import deepcopy

# Change with your path:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper/')

# Choose a quantity to track: biomass or abundance
entity_tracked = 'biomass'

#==============================================================================
# Data Importation
#==============================================================================

phyto = pd.read_csv(join('Data','INSITU','PFG_' + entity_tracked + '.csv'),\
                    parse_dates = ['date'], index_col = 'date') 
    
# Import the events temporal boundaries
event_df = pd.read_csv(join('Results','event_T_anomalies_suited.csv'),\
       parse_dates=['WUI_start', 'WUI_end', 'T_start',\
       'T_end', 'Tmax', 'Tmin', 'window_start', 'window_end']) 
   
# Import the stratification dates
with open('Results/date_pickles/stratify.pickle', 'rb') as f1:
     stratif_dict = pickle.load(f1)
     
strat = stratif_dict['Stratified']
unstrat = stratif_dict['Unstratified']

units = {'loss': '($d^{-1}$)', 'mu': '($d^{-1}$)', 'biomass': '($mgC.mL^{-1}$)',\
         'abundance': '($cells.m{L}^{-1}$)', 'biovolume': '($\mu^3$)',\
         'NPP': '($mgC.m^{-3}d^{-1}$)'}
    
freq = {'loss': '1H', 'mu': '1H', 'biomass': '2H',\
         'abundance': '2H', 'biovolume': '2H',\
         'NPP': '1H'}
    
conv = {entity: 1 for entity in list(units.keys())}
conv['abundance'] = 10 ** 3
conv['biomass'] = 10 ** -9

pfg_list = ['ORGNANO', 'ORGPICOPRO', 'REDNANO', 'REDPICOEUK', 'REDPICOPRO']
        
#=========================================================================
# Get the indices of the sn in stratified wo up, with up and unstratified
#=========================================================================

#*****************
# Get the indices
#*****************
begin_col = 'Tmax'
end_col = 'T_end'

nb_measures = len(phyto)
unstrat_indices = np.any([(phyto.index >= s[0]) & (phyto.index < s[1]) for s in unstrat], 0)
unstrat_indices = pd.DataFrame(data = unstrat_indices, columns = ['valid'], index = phyto.index)

strat_indices = ~unstrat_indices    
stratNonUp_indices = deepcopy(strat_indices)
stratUp_indices = pd.DataFrame(index = phyto.index, data = False, columns = ['valid'])
    
for idx, event in event_df.iterrows():
    stratNonUp_indices.loc[:,'valid'] = np.where(pd.Series(stratNonUp_indices.index).between(event['window_start'],\
                                           event['window_end']),\
                                           np.full(nb_measures, False), stratNonUp_indices['valid'])
        
    stratUp_indices.loc[:,'valid'] = np.where(pd.Series(stratUp_indices.index).between(event[begin_col],\
                                           event[end_col]),\
                                           np.full(nb_measures, True), stratUp_indices['valid'])
 
phyto['Stratify_period'] = np.where(unstrat_indices, 'Unstratified', np.nan)
phyto['Stratify_period'] = np.where(stratNonUp_indices['valid'].values, np.full(nb_measures, 'Stratified (Non-events)'),\
                                    phyto['Stratify_period'])
phyto['Stratify_period'] = np.where(stratUp_indices['valid'].values, np.full(nb_measures, 'Stratified (events)'),\
                                    phyto['Stratify_period'])

phyto = phyto[phyto['Stratify_period'].isin(['Stratified (events)', 'Stratified (Non-events)', 'Unstratified'])]  

phyto[pfg_list] *= conv[entity_tracked]

pd.set_option('display.float_format', '{:.3g}'.format)
median = phyto.groupby('Stratify_period').median()#.round(2)
q3 = phyto.groupby('Stratify_period').quantile(0.75)
q1 = phyto.groupby('Stratify_period').quantile(0.25)
iqr = (q3 - q1)#.round(2)

stats = (median.astype(str) + " (" + iqr.astype(str) + ')').T

stats = stats[['Unstratified', 'Stratified (events)', 'Stratified (Non-events)']]
stats.index = stats.index.str.capitalize()

stats.to_csv(join('Results', entity_tracked + '_stratification.csv'))

