# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:27:13 2021

@author: rfuchs
"""

import os
import numpy as np
import pandas as pd
from os.path import join

# Change the path with yours:
os.chdir('C:/Users/rfuchs/Documents/GitHub/PhytoUpwelling_paper')

#==========================================================
# Per day growth rate
#==========================================================

ssPop_params_dir = 'C:/Users/rfuchs/Documents/These/Oceano/Upwellings/results/growth_rates/'

#pfgs = ['ORGNANO', 'ORGPICOPRO', 'REDNANO', 'REDPICOEUK', 'REDPICOPRO']
pfgs = ['OraNano','OraPicoProk','RedNano','RedPico','RedPicoProk']

mu_all = pd.DataFrame()
for idx, pfg in enumerate(pfgs):
    # npp in mgC 
    mu_path = join('Results', 'growth_rates', 'mu', pfg + '.csv') 
    mu = pd.read_csv(mu_path, parse_dates = ['date'])
    mu.columns = ['date', pfg]
    
    if idx == 0:
        mu_all = pd.concat([mu_all, mu], axis = 1)
    else:
        mu_all = pd.concat([mu_all, mu[pfg]], axis = 1)

mu_all.columns = ['date'] + pfgs
mu_all.to_csv(join('Results', 'PFG_mu.csv'), index = False)
    
#==========================================================
# Per day loss rate
#==========================================================

loss_all = pd.DataFrame()
for idx, pfg in enumerate(pfgs):
    loss = pd.read_csv(ssPop_params_dir + 'loss/' + pfg + '.csv', parse_dates = ['date'])
    loss.columns = ['date', pfg]
    
    if idx == 0:
        loss_all = pd.concat([loss_all, loss], axis = 1)
    else:
        loss_all = pd.concat([loss_all, loss[pfg]], axis = 1)
 
# Write the hourly loss rate
loss_all.to_csv(join('Results', 'PFG_loss.csv'), index = False)


#==========================================================
# Mimick rupture files for mu and loss (no detection associated to these quantities)
#==========================================================

# Import the rupture points
template = pd.read_csv(join('Results', 'responses', 'biomass.csv'),\
                        parse_dates=['WUI_start', 'WUI_end', 'T_start',\
                                     'T_end', 'Tmax', 'Tmin']) 

template.iloc[:, 8:] = np.nan
template.to_csv(join('Results', 'responses', 'mu.csv'), index = False)
template.to_csv(join('Results', 'responses', 'loss.csv'), index = False)
