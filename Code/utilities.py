# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 16:53:39 2021

@author: rfuchs
"""

import numpy as np
import pandas as pd
from copy import deepcopy

#======================================================
# Events handling utilities
#======================================================

def extract_contig_events(dates):
    '''
    Extract upwelling/ downwelling consistent events from a list of dates

    ----------
    dates : list of datetime64
        The dates during which events occur

    Returns
    -------
    merged_events : list of list
        A list of events (each event is characterised by its beginning and end date)

    '''
    
    valid_dates = dates[(dates.index >= pd.to_datetime('2019-09-18 14:00:00')) & \
                        (dates.index <= pd.to_datetime('2020-12-31 23:00:00'))]
    
    contig_dates = valid_dates.index.to_series().diff() == pd.Timedelta('1H')
    nb_events = (contig_dates.to_numpy() == False).sum()
    contig_dates = zip(contig_dates.index, contig_dates.to_numpy())
    
    event_dates_list = [[] for event in range(nb_events)]
    
    # Isolate each event in several DataFrames
    event_counter = -1
    for date, bool_date in contig_dates:
        if bool_date == False:
            event_counter += 1
    
        event_dates_list[event_counter].append(date)
        
    
    # Add hours before and after the events as context: 
    for idx, event in enumerate(event_dates_list):
      # Define the date range to keep
      start = event[0] - 30 * pd.Timedelta('2H')
      end = event[-1] + 50 * pd.Timedelta('2H')
      end = end if end <= pd.to_datetime('2020-12-31 23:00:00') else pd.to_datetime('2020-12-31 23:00:00')   
      
      event_dates_list[idx] = [start, end]
          
    # Merge the events separated with less than 24 hours
    merged_events = [[] for event in range(nb_events)]
    merged_events[0] = deepcopy(event_dates_list[0])
    event_counter = 0
    
    for event in event_dates_list[1:]:
        current_event_end = merged_events[event_counter][-1]
        next_event_start = event[0]
    
        # If the two events are too close they are merged:
        if current_event_end >= next_event_start - pd.Timedelta('24H'):
            merged_events[event_counter] = [merged_events[event_counter][0], event[-1]]
                                                         
        # Otherwise it is stored as such as the next event
        else:
            merged_events[event_counter + 1] = event
            event_counter +=1
    
    # Cleaning it up
    merged_events = [event for event in merged_events if len(event) > 0]
    nb_events = len(merged_events)
    
    return merged_events


def is_stratified(date, stratif_dict):
    '''
    Determine if the date is in a water column stratified period or not 

    Parameters
    ----------
    date : pd.datetime
        The date to test for stratification
    stratif_dict : dict
        The dict that contains the dates of the different stratified periods.

    Returns
    -------
    is_stratif : Boolean
        True if the date belongs to a stratified period, False otherwise

    '''
    is_stratif = np.any([np.all((date >= pd.to_datetime(stratif_dates[0])) &\
    (date <= pd.to_datetime(stratif_dates[1]))) \
    for stratif_dates in stratif_dict['Stratified']])
    
    return is_stratif
    
                    

def split_sequences(valid_dates):
    '''
    Create a list of lists from a list using Nan as cutting value

    Parameters
    ----------
    valid_dates : list
        The list to split.

    Returns
    -------
    event_dates_list : list of lists
        The desired splitted lists.
    '''
    
    contig_dates = valid_dates.index.to_series().diff() == pd.Timedelta('1H')
    
    nb_events = (contig_dates.to_numpy() == False).sum()
    contig_dates = zip(contig_dates.index, contig_dates.to_numpy())
    
    event_dates_list = [[] for event in range(nb_events)]
    
    # Isolate each event in several DataFrames
    event_counter = -1
    for date, bool_date in contig_dates:
        if bool_date == False:
            event_counter += 1
    
        event_dates_list[event_counter].append(date)
        
    return event_dates_list

def merge_close_events(events):
    '''
    Merge the upwelling events separated with less than 24hours

    Parameters
    ----------
    events : list of lists
        The list of dates for each event.

    Returns
    -------
    merged_events : list of lists
        The list where too close events were merged.

    '''
    
    nb_events = len(events)
    merged_events = [[] for event in range(nb_events)]
    merged_events[0] = [events[0][0], events[0][-1]]
    event_counter = 0
        
    for event in events[1:]:
        current_event_end = merged_events[event_counter][-1]
        next_event_start = event[0]
    
        # If the two events are too close they are merged:
        if current_event_end >= next_event_start - pd.Timedelta('24H'):
            merged_events[event_counter] = [merged_events[event_counter][0], event[-1]]
                                                         
        # Otherwise it is stored as such as the next event
        else:
            merged_events[event_counter + 1] = [event[0], event[-1]]
            event_counter +=1
        
    # Cleaning it up
    merged_events = [event for event in merged_events if len(event) > 0]
        
    return merged_events