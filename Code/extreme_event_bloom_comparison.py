# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 13:47:16 2021

@author: rfuchs
"""

import os
import pickle
import pandas as pd
from numpy import trapz
from os.path import join

# Change with your path:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper/')

from Code.utilities import is_stratified
pfgs = ['ORGNANO', 'ORGPICOPRO', 'REDNANO', 'REDPICOEUK', 'REDPICOPRO']

#==============================================================
# Import data
#==============================================================

# Use mean data rather than median data
reactions = pd.read_csv('Results/responses/biomass.csv',\
            parse_dates = ['WUI_start', 'WUI_end', 'T_start', 'T_end',\
                   'Tmax', 'Tmin', 'window_start', 'window_end',\
                   'ORGNANO_start', 'ORGNANO_end', 'ORGPICOPRO_start',
                   'ORGPICOPRO_end', 'REDNANO_start', 'REDNANO_end',\
                   'REDPICOEUK_start', 'REDPICOEUK_end',\
                   'REDPICOPRO_start', 'REDPICOPRO_end'])


biomass = pd.read_csv('Data/INSITU/PFG_biomass.csv', parse_dates = ['date'],\
            index_col='date')
   
# Import the stratification dates
pickle_dir = join('Results', 'date_pickles')
with open(join(pickle_dir,'stratify.pickle'), 'rb') as f1:
     stratif_dict = pickle.load(f1)

# Import the upwelling dates
with open(join(pickle_dir,'welling_events.pickle'), 'rb') as f1:
     welling_dict = pickle.load(f1)
  
# Keep stratfied events only
welling_dict = [w for w in welling_dict if is_stratified(w[0], stratif_dict)]

# Import Upwelling and Downwelling dates 
with open(join(pickle_dir,'blooms_lims.pickle'), 'rb') as f1:
     blooms_lims = pickle.load(f1)

with open(join(pickle_dir,'blooms_refs.pickle'), 'rb') as f1:
     blooms_refs = pickle.load(f1)

#==============================================================
# Compute the cumulated biomass increase due to events
#==============================================================

dbiomass_dt = pd.DataFrame(index = reactions['window_start'], columns = pfgs)
dbiomass = pd.DataFrame(index = reactions['window_start'], columns = pfgs)

for idx, event in reactions.iterrows():
    for pfg in pfgs:
        
        if event[pfg + '_phyto_react_before_T']:
            continue
        
        # Time deltas
        dt1 = (event[pfg + '_end'] - event[pfg + '_start']).total_seconds() / 3600
        dt_ref = (event[pfg + '_start'] - event['window_start']).total_seconds() / 3600
        
        # Compute the integrals of the reference and the pfg reaction
        pfg_data = biomass.loc[event[pfg + '_start']: event[pfg + '_end'], pfg].interpolate(method = 'polynomial', order = 1)
        b_pfg = trapz(pfg_data.dropna()) 
        
        ref_data = biomass.loc[event['window_start']: event[pfg + '_start'], pfg].interpolate(method = 'polynomial', order = 1)
        b_ref = trapz(ref_data.dropna())
        
        dbiomass.loc[event['window_start'], pfg] = (b_pfg - b_ref) 

        b_ref = dt1 * b_ref / dt_ref # Make it comparable
        dbiomass_dt.loc[event['window_start'], pfg] = (b_pfg - b_ref) / dt1
        
    
#==============================================================
# Compute the biomass change due to blooms 
#==============================================================

dblooms_dt = {'2020': [], '2021': []}
dblooms = {'2020': [], '2021': []}

for year in ['2020', '2021']:
    refs = blooms_refs[year]
    blooms_y = blooms_lims[year]
    for bloom_idx, bloom in enumerate(blooms_y):
        
        # Compute the dts
        dt_bloom = (bloom[1] - bloom[0]).total_seconds() / 3600
        dt_ref = (refs[bloom_idx][1] - refs[bloom_idx][0]).total_seconds() / 3600

        # Compute the area under the curve and substract the ref 
        b_dt = biomass.loc[bloom[0]:bloom[1]].interpolate(method = 'polynomial', order = 1)
        b_dt = b_dt.apply(lambda x: trapz(x.dropna())) 
        
        b_ref = biomass.loc[refs[bloom_idx][0]:refs[bloom_idx][1]].interpolate(method = 'polynomial', order = 1)
        b_ref = b_ref.apply(lambda x: trapz(x.dropna()))
        dblooms[year].append((b_dt - b_ref).sum())

        b_ref = dt_bloom * b_ref / dt_ref
        dblooms_dt[year].append((b_dt - b_ref).sum() / dt_bloom)
 
#==============================================================
# Results
#==============================================================


b_ratio = dbiomass.loc[pd.to_datetime('2019-11-06 17:00:00')].sum() / dblooms['2020'][0]
b_ratio = (100 * b_ratio).round(1)
print('The total increase in biomass during the biggest event was of', b_ratio,\
      '% with respect to the corresponding year spring bloom')
    

b_dt_ratio = dbiomass_dt.loc[pd.to_datetime('2020-06-18 01:00:00')].sum() / dblooms_dt['2020'][0]
b_dt_ratio = (100 * b_dt_ratio).round(1)
print('The per unit of time increase in biomass during the biggest event was of', b_dt_ratio,\
      '% with respect to the corresponding year spring bloom')
    

# The 2020-09-23 19:00:00 event presented an even higher increase but presented too many NAs
a = dbiomass_dt.loc[pd.to_datetime('2020-01-01 01:00:00'):pd.to_datetime('2021-01-01 01:00:00')].sum(1) / dblooms_dt['2020'][0]
b = dbiomass_dt.loc[pd.to_datetime('2021-01-01 01:00:00'):pd.to_datetime('2022-01-01 01:00:00')].sum(1) / dblooms_dt['2021'][0]
a.max()
