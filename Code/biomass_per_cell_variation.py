# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 14:43:01 2022

@author: rfuchs
"""

import os
import pandas as pd
from os.path import join

pfgs = ['ORGNANO', 'ORGPICOPRO', 'REDNANO', 'REDPICOEUK', 'REDPICOPRO']

dates_col = ['WUI_start', 'WUI_end', 'T_start', 'T_end', 'Tmax', 'Tmin',
       'window_start', 'window_end', 'ORGNANO_start',
       'ORGNANO_end', 'ORGPICOPRO_start',
       'ORGPICOPRO_end', 'REDNANO_start', 'REDNANO_end',
       'REDPICOEUK_start',
       'REDPICOEUK_end', 'REDPICOPRO_start', 'REDPICOPRO_end']

os.chdir(r'C:\Users\rfuchs\Documents\GitHub\PhytoUpwelling_paper')

## Total Biomass and Biomass per cell
biomass = pd.read_csv('Data/INSITU/PFG_biomass.csv', parse_dates=['date'], index_col = ['date'])
abundance = pd.read_csv('Data/INSITU/PFG_abundance.csv', parse_dates=['date'], index_col = ['date'])
biomass_per_cell = biomass.div(abundance) 

path = join('Results', 'responses', 'biomass' + '.csv')
event_df = pd.read_csv(path, parse_dates = dates_col)   


vars_during_b = pd.DataFrame()
vars_after_b = pd.DataFrame()

vars_during_bC = pd.DataFrame()
vars_after_bC = pd.DataFrame()

vars_during_a = pd.DataFrame()
vars_after_a = pd.DataFrame()

for event_idx, event in event_df.iterrows():
    biomCell_data = biomass_per_cell.loc[event['window_start']:event['window_end']]
    biom_data = biomass.loc[event['window_start']:event['window_end']]
    abund_data = abundance.loc[event['window_start']:event['window_end']]

    var_during_pfgs_b = []
    var_after_pfgs_b = []
    
    var_during_pfgs_bC = []
    var_after_pfgs_bC = []
    
    var_during_pfgs_a = []
    var_after_pfgs_a = []
    
    for pfg in pfgs:
        biomCell_pfg_data = biomCell_data[[pfg]].dropna()
        biom_pfg_data = biom_data[[pfg]].dropna()
        abund_pfg_data = abund_data[[pfg]].dropna()

        #==================================================
        # Biomass
        #==================================================
        
        median_before = biom_pfg_data.loc[biom_pfg_data.index < event[pfg + '_start'], pfg].mean()
        median_during = biom_pfg_data.loc[(biom_pfg_data.index >= event[pfg + '_start']) \
                               & (biom_pfg_data.index < event[pfg + '_end']), pfg].mean()
        median_after = biom_pfg_data.loc[biom_pfg_data.index >= event[pfg + '_end'], pfg].mean()
        
        var_during_biom = (median_during - median_before) / median_before
        var_after_biom = (median_after - median_during) / median_during 
        
        #==================================================
        # Biomass/cell
        #==================================================

        median_before = biomCell_pfg_data.loc[biomCell_pfg_data.index < event[pfg + '_start'], pfg].mean()
        median_during = biomCell_pfg_data.loc[(biomCell_pfg_data.index >= event[pfg + '_start']) \
                               & (biomCell_pfg_data.index < event[pfg + '_end']), pfg].mean()
        median_after = biomCell_pfg_data.loc[biomCell_pfg_data.index >= event[pfg + '_end'], pfg].mean()
        
        var_during_biomCell = (median_during - median_before) / median_before
        var_after_biomCell = (median_after - median_during) / median_during
                
        #==================================================
        # Abundance
        #==================================================
        median_before = abund_pfg_data.loc[abund_pfg_data.index < event[pfg + '_start'], pfg].mean()
        median_during = abund_pfg_data.loc[(abund_pfg_data.index >= event[pfg + '_start']) \
                               & (abund_pfg_data.index < event[pfg + '_end']), pfg].mean()
        median_after = abund_pfg_data.loc[abund_pfg_data.index >= event[pfg + '_end'], pfg].mean()
        
        var_during_abund = (median_during - median_before) / median_before
        var_after_abund = (median_after - median_during) / median_during

        #==================================================
        # Wrap up
        #==================================================
          
        var_during_pfgs_b.append(var_during_biom)
        var_after_pfgs_b.append(var_after_biom)
        
        var_during_pfgs_bC.append(var_during_biomCell)
        var_after_pfgs_bC.append(var_after_biomCell)
        
        var_during_pfgs_a.append(var_during_abund)
        var_after_pfgs_a.append(var_after_abund)
        
        
    vars_during_b = pd.concat([vars_during_b, pd.DataFrame(var_during_pfgs_b).T])
    vars_after_b = pd.concat([vars_after_b, pd.DataFrame(var_after_pfgs_b).T])
    
    vars_during_bC = pd.concat([vars_during_bC, pd.DataFrame(var_during_pfgs_bC).T])
    vars_after_bC = pd.concat([vars_after_bC, pd.DataFrame(var_after_pfgs_bC).T])
    
    vars_during_a = pd.concat([vars_during_a, pd.DataFrame(var_during_pfgs_a).T])
    vars_after_a = pd.concat([vars_after_a, pd.DataFrame(var_after_pfgs_a).T])

vars_during_b.columns = pfgs
vars_after_b.columns = pfgs

vars_during_bC.columns = pfgs
vars_after_bC.columns = pfgs

vars_during_a.columns = pfgs
vars_after_a.columns = pfgs

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

cols = ['During: Biomass vs. Biomass/cell', 'After: Biomass vs. Biomass/cell',\
        'During: Biomass vs. Abundance', 'After: Biomass vs. Abundance']
corrs = pd.DataFrame(index = pfgs, columns = cols)
pvals = pd.DataFrame(index = pfgs, columns = cols)
 

for pfg in pfgs:
    data = pd.DataFrame([vars_during_b[pfg].values, vars_during_bC[pfg].values]).T.dropna()
    plt.scatter(data[0], data[1])
    r = spearmanr(data[0], data[1])
    plt.title('BC During' + pfg + ' corr: ' + str(r.correlation.round(3)) + ', p-value: ' + str(r.pvalue.round(3)))
    plt.ylabel('Biomass/cell variation during')
    plt.xlabel('Biomass variation during')
    plt.show()
    
    corrs.loc[pfg, 'During: Biomass vs. Biomass/cell'] = r.correlation.round(3)
    pvals.loc[pfg, 'During: Biomass vs. Biomass/cell'] = r.pvalue.round(3)
    
    data = pd.DataFrame([vars_after_b[pfg].values, vars_after_bC[pfg].values]).T.dropna()
    plt.scatter(data[0], data[1])
    r = spearmanr(data[0], data[1])
    plt.title('BC After' + pfg + 'corr: ' + str(r.correlation.round(3)) + ', p-value: ' + str(r.pvalue.round(3)))
    plt.ylabel('Biomass/cell variation during')
    plt.xlabel('Biomass variation during')
    plt.show()
    
    corrs.loc[pfg, 'After: Biomass vs. Biomass/cell'] = r.correlation.round(3)
    pvals.loc[pfg, 'After: Biomass vs. Biomass/cell'] = r.pvalue.round(3)

    
    data = pd.DataFrame([vars_during_b[pfg].values, vars_during_a[pfg].values]).T.dropna()
    plt.scatter(data[0], data[1])
    r = spearmanr(data[0], data[1])
    plt.title('A During' + pfg + ' corr: ' + str(r.correlation.round(3)) + ', p-value: ' + str(r.pvalue.round(3)))
    plt.ylabel('Biomass/cell variation during')
    plt.xlabel('Biomass variation during')
    plt.show()
    
    corrs.loc[pfg, 'During: Biomass vs. Abundance'] = r.correlation.round(3)
    pvals.loc[pfg, 'During: Biomass vs. Abundance'] = r.pvalue.round(3)
    
    data = pd.DataFrame([vars_after_b[pfg].values, vars_after_a[pfg].values]).T.dropna()
    plt.scatter(data[0], data[1])
    r = spearmanr(data[0], data[1])
    plt.title('A After' + pfg + 'corr: ' + str(r.correlation.round(3)) + ', p-value: ' + str(r.pvalue.round(3)))
    plt.ylabel('Biomass/cell variation during')
    plt.xlabel('Biomass variation during')
    plt.show()
    
    corrs.loc[pfg, 'After: Biomass vs. Abundance'] = r.correlation.round(3)
    pvals.loc[pfg, 'After: Biomass vs. Abundance'] = r.pvalue.round(3)


corrs = np.where(pvals > 0.10, np.nan, corrs)
corrs = pd.DataFrame(corrs, index = pvals.index, columns = pvals.columns)


corrs.to_csv(join('Results', 'corrs_biomass_abundance_cell.csv'))


vars_during.median(0)
vars_after.median(0)

decrease_during = pd.DataFrame((vars_during < 0.5).mean(0))
decrease_during.columns = ['decrease_during']
increase_after = pd.DataFrame((vars_after > 0.5).mean(0))
increase_after.columns = ['increase_after']
decrease_during.join(increase_after)
