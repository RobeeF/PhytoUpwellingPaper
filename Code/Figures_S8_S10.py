# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 16:23:54 2022

@author: rfuchs
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 14:24:16 2022

@author: rfuchs
"""

import os
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
from os.path import join
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from scipy.stats import spearmanr

# Choose a quantity to track: biomass, abundance
entity_tracked = 'abundance'

# Change with your path:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper/')

#==============================================================================
# Data Importation
#==============================================================================

path = join('Data', 'INSITU' , 'PFG_' + entity_tracked + '.csv')
phyto = pd.read_csv(path, parse_dates = ['date'],\
                    index_col = 'date') 
  
phyto = phyto.sort_index()
gr =  pd.read_csv(join('Results','PFG_mu.csv'), parse_dates = ['date'],\
                    index_col = 'date') 
loss =  pd.read_csv(join('Results','PFG_loss.csv'), parse_dates = ['date'],\
                    index_col = 'date') 
   

phyto = phyto.sort_index()
loss = loss.sort_index()
gr = gr.sort_index()

pfg_list = gr.columns
         
# Import the events temporal boundaries
dates_col = ['WUI_start', 'WUI_end', 'T_start', 'T_end', 'Tmax', 'Tmin',
       'window_start', 'window_end', 'ORGNANO_start',
       'ORGNANO_end', 'ORGPICOPRO_start',
       'ORGPICOPRO_end', 'REDNANO_start', 'REDNANO_end',
       'REDPICOEUK_start',
       'REDPICOEUK_end', 'REDPICOPRO_start', 'REDPICOPRO_end']

event_df = pd.read_csv(join('Results', 'responses', entity_tracked + '.csv'),\
                       parse_dates=dates_col)     
    
# Import the stratification dates
with open(join('Results', 'date_pickles', 'stratify.pickle'), 'rb') as f1:
     stratif_dict = pickle.load(f1)
     
strat = stratif_dict['Stratified']
unstrat = stratif_dict['Unstratified']

stps = pd.read_csv(join('Data', 'INSITU','stps_1H.csv'),\
                   parse_dates = ['date'], index_col = 'date')

pfg_colors = dict(zip(['Orgnano', 'Orgpicopro', 'Rednano', 'Redpicoeuk', 'Redpicopro'],\
                      ['#9467bd', '#1f77b4', '#2ca02c', '#d62728', 'orange']))


#=========================================================================
# Differenciating between reaction and relaxation 
#=========================================================================

gr_T = pd.DataFrame(index = pfg_list.str.capitalize(),\
                    columns = ['Pre-reaction', 'Reaction', 'Relaxation'])
pvals_T = pd.DataFrame(index = pfg_list.str.capitalize(),\
                    columns = ['Pre-reaction', 'Reaction', 'Relaxation'])

gr_ent = pd.DataFrame(index = pfg_list.str.capitalize(),\
                    columns = ['Pre-reaction', 'Reaction', 'Relaxation'])
pvals_ent = pd.DataFrame(index = pfg_list.str.capitalize(),\
                    columns = ['Pre-reaction', 'Reaction', 'Relaxation'])

for pfg in pfg_list:
    print(pfg)
    
    gr_entity_T = gr[[pfg]].join(phyto[pfg], how = 'inner', lsuffix = '_growth rate',\
                               rsuffix = '_' + entity_tracked).join(stps['T'], how = 'inner')
    gr_entity_T['phase'] = np.nan
    
    for event_idx, event in event_df.iterrows():
        if event[pfg + '_phyto_react_before_T']:
            continue
        
        gr_entity_T.loc[event['window_start']:event[pfg + '_start'], 'phase'] = 'Pre-reaction'
        gr_entity_T.loc[event[pfg + '_start']:event[pfg + '_end'], 'phase'] = 'Reaction'
        gr_entity_T.loc[event[pfg + '_end']:event['window_end'], 'phase'] = 'Relaxation'
        
   
    gr_entity_T = gr_entity_T.dropna()
    
    
    #*********************************
    # Growth rate vs Temperature
    #*********************************
    
    for phase, group in gr_entity_T.groupby('phase'):
        
        sp_test = spearmanr(group[pfg + '_growth rate'], group['T'])
        gr_T.loc[pfg.capitalize(), phase] = sp_test.correlation
        pvals_T.loc[pfg.capitalize(), phase] = sp_test.pvalue     
        
     
    #*********************************
    # Growth rate vs Entity
    #*********************************
    
    for phase, group in gr_entity_T.groupby('phase'):
        
        sp_test = spearmanr(group[pfg + '_growth rate'], group[pfg + '_' + entity_tracked])
        gr_ent.loc[pfg.capitalize(), phase] = sp_test.correlation
        pvals_ent.loc[pfg.capitalize(), phase] = sp_test.pvalue     
        

#=====================================================================
# Compare the variation in growth rates between the different phases
#=====================================================================

median_gr_phase = pd.DataFrame(index = pfg_list.str.capitalize(),\
                               columns = ['Pre-reaction', 'Reaction', 'Relaxation'])

for pfg in pfg_list:
    gg = pd.DataFrame()

    for event_idx, event in event_df.iterrows():
        if event[pfg + '_phyto_react_before_T']:
            continue
        
        prereac = gr.loc[event['window_start']:event[pfg + '_start'], pfg].median() 
        react = gr.loc[event[pfg + '_start']:event[pfg + '_end'], pfg].median() 
        relax = gr.loc[event[pfg + '_end']:event['window_end'], pfg].median()
        
        gg = gg.append({'window_start': event['window_start'], 'Pre-reaction': prereac,\
                   'Reaction': react, 'Relaxation': relax}, ignore_index = True)
    
    # Clean the dataset
    gg = gg.set_index('window_start').replace([np.inf,\
                                -np.inf, pd.NaT], np.nan).dropna(how = 'all')
    
    median_gr_phase.loc[pfg.capitalize()] = gg.mean()
        
    # Confirm with a Kruskal
    all_ = gg[['Pre-reaction', 'Reaction', 'Relaxation']].dropna(how = 'any')
    kw3 = kruskal(all_['Pre-reaction'], all_['Reaction'], all_['Relaxation']).pvalue
 
    print(pfg)
    #print(kw)
    #print(kw2)
    print(kw3)
    print(gg.mean())
    print('-------------------------')
    
    
#========================================================
# Loss vs mu 
#=========================================================

corr = pd.DataFrame(index = pfg_list.str.capitalize(),\
                    columns = ['Pre-reaction', 'Reaction', 'Relaxation'])
pvals = pd.DataFrame(index = pfg_list.str.capitalize(),\
                    columns = ['Pre-reaction', 'Reaction', 'Relaxation'])

for pfg in pfg_list:
    print(pfg)
    mu_loss = gr[[pfg]].join(loss[[pfg]], lsuffix = '_gr', rsuffix = '_loss')    
    mu_loss['phase'] = np.nan
    
    for event_idx, event in event_df.iterrows():
        if event[pfg + '_phyto_react_before_T']:
            continue
        
        mu_loss.loc[event['window_start']:event[pfg + '_start'], 'phase'] = 'Pre-reaction'
        mu_loss.loc[event[pfg + '_start']:event[pfg + '_end'], 'phase'] = 'Reaction'
        mu_loss.loc[event[pfg + '_end']:event['window_end'], 'phase'] = 'Relaxation'
        
    mu_loss = mu_loss.dropna()
    
    for phase, group in mu_loss.groupby('phase'):
        sp_test = spearmanr(group[pfg + '_gr'], group[pfg + '_loss'])
        corr.loc[pfg.capitalize(), phase] = sp_test.correlation
        pvals.loc[pfg.capitalize(), phase] = sp_test.pvalue

#===============================
# Heatmaps
#===============================

sns.heatmap(24 * median_gr_phase.replace([np.nan], 0),\
            annot = True, vmax = 1.3)
plt.title(entity_tracked.capitalize())
plt.savefig(join('Figures', entity_tracked + '_division_per_phase.png'))
plt.show()


sns.heatmap(corr.replace([np.nan], 0), mask = pvals.replace([np.nan], 1) > 0.05,\
            annot = True, vmin = -0.3, vmax = 0.6)
plt.title(entity_tracked.capitalize())
plt.savefig(join('Figures', entity_tracked + '_division_vs_loss.png'))
plt.show()


