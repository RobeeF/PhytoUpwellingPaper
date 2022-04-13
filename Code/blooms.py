# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 11:44:54 2022

@author: rfuchs
"""

import os
import pickle
import numpy as np
import pandas as pd
from os.path import join
from scipy.signal import butter, filtfilt, find_peaks

# Change with your path:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper/')

#==============================================================================
# Data Importation
#==============================================================================

# Choose a quantity to track: biomass, abundance, mu, NPP, biovolume
units = {'loss': '($d^{-1}$)', 'mu': '($d^{-1}$)', 'biomass': '($fgC.\mu{L}^{-1}$)',\
         'abundance': '($cells.\mu{L}^{-1}$)', 'biovolume': '($\mu^3$)',\
         'NPP': '($mgC.m^{-3}d^{-1}$)'}
freq = {'loss': '1H', 'mu': '1H', 'biomass': '2H',\
         'abundance': '2H', 'biovolume': '2H',\
         'NPP': '1H'}

phyto = pd.read_csv(join('Data','INSITU','PFG_biomass.csv'), parse_dates = ['date'],\
                    index_col = 'date') 
    
# Import the stratification dates
with open(join('Results','date_pickles','stratify.pickle'), 'rb') as f1:
     stratif_dict = pickle.load(f1)
    
stratif_dict['Unstratified'] = stratif_dict['Unstratified'][:-1]

#==============================================================================
# Determining the bloom dates
#==============================================================================

blooms_lims = {}

for unstrat_period in stratif_dict['Unstratified']:
    # Fetch the raw data
    year = str(unstrat_period[1].year)    
    phyto_unstrat = phyto.loc[unstrat_period[0]: unstrat_period[1]]
    phyto_year = phyto.loc[pd.to_datetime(year + '-01-01'):  pd.to_datetime(year + '-12-31')]
    
    # Filter the series
    filter_ = 0.05
    b, a = butter(1, filter_)    
    tb = pd.DataFrame(phyto_unstrat.sum(1), columns = ['biomass'])
    tb['filter'] = filtfilt(b, a, x=tb['biomass'], method='gust')
    
    # Find the maxima
    threshold = phyto_year.sum(1).median(0) * 1.05
    tb['up'] = tb['filter'] > threshold
    peaks = find_peaks(tb['filter'], prominence = 50000, distance = 50)[0]
    peak_indices = tb.index[peaks]
    
    # Compute the pre and post bloom indices
    pre_indices = []
    for peak_idx in peak_indices:
        pre_idx = tb.loc[:peak_idx, 'up'].sort_index(ascending = False).idxmin()
        pre_indices.append(pre_idx)

    post_indices = []
    for peak_idx in peak_indices:
        post_idx = tb.loc[peak_idx:, 'up'].sort_index().idxmin()
        post_indices.append(post_idx)
        
    # Wrap up the results
    nb_blooms = len(peak_indices)
    indices = []
    for bloom_nb in range(nb_blooms):
        indices.append([pre_indices[bloom_nb], post_indices[bloom_nb]])
    blooms_lims[year] = indices

# Determining the references periods
blooms_refs = {}
for year in ['2020', '2021']:
    blooms_refs[year] = [[start - pd.Timedelta('1W'), start] for start, end in blooms_lims[year]]
   
# Store the results
with open(join('Results','date_pickles','blooms_lims.pickle'), 'wb') as f1:
    pickle.dump(blooms_lims, f1)
    
with open(join('Results','date_pickles','blooms_refs.pickle'), 'wb') as f1:
    pickle.dump(blooms_refs, f1)
    
#==============================================================================
# Taking the 90% of the bloom distribution
#==============================================================================

entity_tracked = 'abundance'
phyto = pd.read_csv(join('Data','INSITU','PFG_' + entity_tracked + '.csv'),\
                    parse_dates = ['date'], index_col = 'date') 

thr = 0.5
blooms = pd.DataFrame(index = phyto.columns, columns = ['2020', '2021'])

for unstrat_period in stratif_dict['Unstratified']:
    year = str(unstrat_period[1].year)
    bloom_year = blooms_lims[year]
    
    phyto_blooms = pd.DataFrame()
    for bloom_dates in bloom_year:
        phyto_bloom = phyto.loc[bloom_dates[0]: bloom_dates[1]]
        phyto_blooms = phyto_blooms.append(phyto_bloom)
    
    phyto_unstrat = phyto.loc[unstrat_period[0]: unstrat_period[1]]
    phyto_year = phyto.loc[pd.to_datetime(year + '-01-01'):  pd.to_datetime(year + '-12-31')]

    for pfg in phyto_unstrat.columns:
        data = phyto_blooms[pfg]        
   
        blooms_y = blooms_lims[year]
        refs = blooms_refs[year]
        refs2 = phyto_year.median(0)

        blooms.loc[pfg, year] = data.quantile([thr]).values[0]
        

blooms['2019'] = np.nan
blooms.index.name = 'PFG'
blooms.to_csv(join('Results', 'blooms', entity_tracked + str(int(thr * 100)) + '.csv'))

    
