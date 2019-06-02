# Before running script use Excel to remove any "--" 

import os
import pandas as pd
import csv
import numpy as np
import scipy as sp
import scipy.stats
import datetime
import math
 
#---------------------------#
#     Assign Variables      #
#---------------------------#

file_in = 'camk_filtered_raw.csv' 
dir_in = 'C:\\Users\\haley\\Dropbox\\Projects\\camk2a\\Raw Data\\Filtered'
#dir_in = 'D:\\Dropbox\\Projects\\camk2a'
os.chdir(dir_in)
now = datetime.datetime.now()
date_string = now.strftime('%y%m%d %H.%M')
#dir_out = dir_in + '\\Calculated ' + date_string
dir_out = dir_in + '\\Calculated'
data = pd.read_csv(file_in)
exp_name = 'camk'

splitter = 20
connection = 0.1
strong = 1.0




#---------------------------#
#     Define Functions      #
#---------------------------#

# Saves data to a new subfolder 
def save_data (data, file_out, dir_out):
    try:
        os.stat(dir_out)
    except:
        os.mkdir(dir_out)

    # Write data to file
    os.chdir(dir_out)
    data.to_csv(file_out, index = False)



# clean the data by converting all negative CC to 0 since by definition, hyperpolarizing one cell cannot lead to depolarization in another
# Validated in Origin
def get_clean(data, factor1, factor2, factor3, factor4, file_out, dir_out):
    data[factor1] = np.where(data[factor1] > 0, 0, data[factor1])
    data[factor2] = np.where(data[factor2] > 0, 0, data[factor2])
    data[factor3] = np.where(data[factor3] > 0, 0, data[factor3])
    data[factor4] = np.where(data[factor4] > 0, 0, data[factor4])
    save_data (data, file_out, dir_out)
    return data

# Calculates Rin from CC data (not from VSteps to be consistent with xfer)
# Validated in Origin
def get_Rin(data_Rin, file_out, dir_out):
    data_Rin['rin1'] = np.nan
    data_Rin['rin2'] = np.nan
    data_Rin['rin1_nmda'] = np.nan
    data_Rin['rin2_nmda'] = np.nan

    for index, row in data_Rin.iterrows():	
        if row['istim1'] != 0 and math.isnan(row['istim1']) == False:
            data_Rin.loc[index:index:, 'rin1'] = row['vstim1'] / row['istim1'] * 1000
        if row['istim2'] != 0 and math.isnan(row['istim2']) == False:
            data_Rin.loc[index:index:, 'rin2'] = row['vstim2'] / row['istim2'] * 1000
        if row['istim1_nmda'] != 0 and math.isnan(row['vstim1_nmda']) == False:
            data_Rin.loc[index:index:, 'rin1_nmda'] = row['vstim1_nmda'] / row['istim1_nmda'] * 1000    
        if row['istim2_nmda'] != 0 and math.isnan(row['istim2_nmda']) == False:
            data_Rin.loc[index:index:, 'rin2_nmda'] = row['vstim2_nmda'] / row['istim2_nmda'] * 1000
        
        
    # Save data including intrinsic
    save_data (data_Rin, file_out, dir_out)
    return data_Rin

# Combine 1>2 and 2>1
# Validated in Origin
def get_connections (data_connect, file_out, dir_out):
    connections = pd.DataFrame({'id':[],'strain':[],'distance':[],'vm1':[],'vm2':[],
        'vrec1':[],'vstim1':[],'istim1':[],'vm1_nmda':[],'vm2_nmda':[],'vrec2_nmda':[],'vstim1_nmda':[],'istim1_nmda':[],
        'rin1':[],'rin2':[],'rin1_nmda':[],'rin2_nmda':[]}) 
    for index, row in data_connect.iterrows():
        connections = connections.append({'id':row['id'],'strain':row['strain'],'distance':row['distance'],
            'vrec1':row['vrec1'],'vstim1':row['vstim1'],'istim1':row['istim1'],
            'vrec1_nmda':row['vrec1_nmda'],'vstim1_nmda':row['vstim1_nmda'],'istim1_nmda':row['istim1_nmda'],
            'rin1':row['rin1'],'rin2':row['rin2'],'rin1_nmda':row['rin1_nmda'],'rin2_nmda':row['rin2_nmda'],
            'vm1':row['vm1'],'vm1_nmda':row['vm1_nmda']}, ignore_index = True) 
        connections = connections.append({'id':row['id'],'strain':row['strain'],'distance':row['distance'],
            'vrec1':row['vrec2'],'vstim1':row['vstim2'],'istim1':row['istim2'],
            'vrec1_nmda':row['vrec2_nmda'],'vstim1_nmda':row['vstim2_nmda'],'istim1_nmda':row['istim2_nmda'],
            'rin1':row['rin2'],'rin2':row['rin1'],'rin1_nmda':row['rin2_nmda'],'rin2_nmda':row['rin1_nmda'], 'vm1':row['vm2'],'vm1_nmda':row['vm2_nmda']}, ignore_index = True)
        connections = connections[['id','strain','distance','vrec1','vstim1','istim1',
            'vrec1_nmda','vstim1_nmda','istim1_nmda','rin1','rin2','rin1_nmda','rin2_nmda', 'vm1','vm2','vm1_nmda', 'vm2_nmda']] 
    save_data (connections, file_out, dir_out)
    return connections

# Get CC
# Validated in Origin
def get_cc(data_cc):
    data_cc['cc'] = np.nan
    data_cc['cc_nmda'] = np.nan
    
    # Checking for NaN in vrec rather than vstim because vrec appears first in the calculation
    # Otherwise you would get a type error on the next line because vrec is called first
    for index, row in data_cc.iterrows():
        if row['vstim1'] != 0 and math.isnan(row['vstim1']) == False: 
            data_cc.loc[index:index:, 'cc'] = row['vrec1'] / row['vstim1'] * 100
        if row['vstim1_nmda'] != 0 and math.isnan(row['vrec1_nmda']) == False:
            data_cc.loc[index:index:, 'cc_nmda'] = row['vrec1_nmda'] / row['vstim1_nmda'] * 100
    return data_cc

# Get xfer
# Validated in Origin
def get_xfer(data_xfer):
    data_xfer['xfer'] = np.nan
    data_xfer['xfer_nmda'] = np.nan
    for index, row in data_xfer.iterrows():	
        if row['istim1'] != 0 and math.isnan(row['vrec1']) == False:
            data_xfer.loc[index:index:, 'xfer'] = row['vrec1'] / row['istim1'] * 1000
        if row['istim1_nmda'] != 0 and math.isnan(row['vrec1_nmda']) == False:
            data_xfer.loc[index:index:, 'xfer_nmda'] = row['vrec1_nmda'] / row['istim1_nmda'] * 1000
    return data_xfer

# Get gj
# Validated in Origin
def get_gj(data_gj):
    data_gj['gj'] = np.nan
    data_gj['gj_nmda'] = np.nan
    for index, row in data_gj.iterrows():	
        if row['xfer'] != 0:
            data_gj.loc[index:index:, 'gj'] = 1/(((row['rin1'] * row['rin2'])-(row['xfer']*row['xfer']))/row['xfer']) * 1000000
        if row['xfer_nmda'] != 0:
            data_gj.loc[index:index:, 'gj_nmda'] = 1/(((row['rin1_nmda'] * row['rin2_nmda'])-(row['xfer_nmda']*row['xfer_nmda']))/row['xfer_nmda']) * 1000000
    return data_gj

# Get_change in gj, cc, and Rin
# Validated in Origin
def get_changes(data_change):
    data_change['cc_change'] = np.nan
    data_change['cc_pchange'] = np.nan
    data_change['rin_change'] = np.nan
    data_change['rin_pchange'] = np.nan
    data_change['gj_change'] = np.nan
    data_change['gj_pchange'] = np.nan
    data_change['vm_change'] = np.nan
    for index, row in data_change.iterrows():
        if row['cc'] != 0:
            data_change.loc[index:index:, 'cc_change'] = row['cc_nmda']-row['cc']
            data_change.loc[index:index:,'cc_pchange'] = ((row['cc_nmda']-row['cc'])/row['cc']) * 100
        if row['rin1'] != 0:
            data_change.loc[index:index:,'rin_change'] = row['rin1_nmda']-row['rin1']
            data_change.loc[index:index:,'rin_pchange'] = ((row['rin1_nmda']-row['rin1'])/row['rin1']) * 100
        if row['gj'] != 0:
            data_change.loc[index:index:,'gj_change'] = row['gj_nmda']-row['gj']
            data_change.loc[index:index:,'gj_pchange'] = ((row['gj_nmda']-row['gj'])/row['gj']) * 100
        data_change.loc[index:index:,'vm_change'] = row['vm1_nmda']-row['vm1']
    return data_change

def get_bin (data):
    data_bin = data.copy()
    data_bin['bins'] = np.nan
    for index, row in data_bin.iterrows():
        if row['distance'] >= 0 and row['distance'] < splitter:
            data_bin.loc[index:index:, 'bins'] = splitter/2
        elif row['distance'] >= splitter and row['distance'] < 2*splitter:
            data_bin.loc[index:index:, 'bins'] = 2*splitter-(splitter/2)
        elif row['distance'] >= 2*splitter and row['distance'] < 3*splitter:
            data_bin.loc[index:index:, 'bins'] = 3*splitter - (2*splitter/2)
        elif row['distance'] >= 3*splitter and row['distance'] < 4*splitter:
            data_bin.loc[index:index:, 'bins'] = 4*splitter - (3*splitter/2)
        elif row['distance'] >= 4*splitter and row['distance'] < 5*splitter:
            data_bin.loc[index:index:, 'bins'] = 5*splitter - (4*splitter/2)
        elif row['distance'] >= 5*splitter and row['distance'] < 6*splitter:
            data_bin.loc[index:index:, 'bins'] = 6*splitter - (5*splitter/2)
    return data_bin


def sort_connected (data, file_out, dir_out):
    data_connected = data.copy()
    for index, row in data_connected.iterrows():
        if row['cc'] < connection:
            data_connected.drop(index, inplace = True)
    save_data (data_connected, file_out, dir_out)

def sort_strong (data, file_out, dir_out):
    data_strong = data.copy()
    for index, row in data_strong.iterrows():
        if row['cc'] >= strong and row['cc_nmda'] >= 0:
            pass
        else:
            data_strong.drop(index, inplace = True)
    save_data (data_strong, file_out, dir_out)

def sort_weak (data, file_out, dir_out):
    data_weak = data.copy()
    for index, row in data_weak.iterrows():
        if row['cc'] < strong and row['cc'] >=0.1 and row['cc_nmda'] >= 0:
            pass
        else:
            data_weak.drop(index, inplace = True)
    save_data (data_weak, file_out, dir_out)

def sort_strengthened (data, file_out, dir_out):
    data_strengthened = data.copy()
    for index, row in data_strengthened.iterrows():
        if row['cc'] >= 0.1 and row['cc_nmda'] >= 0.1 and row['cc_pchange'] > 0:
            pass
        else:
            data_strengthened.drop(index, inplace = True)
    save_data (data_strengthened, file_out, dir_out)

def sort_weakened (data, file_out, dir_out):
    data_weakened = data.copy()
    for index, row in data_weakened.iterrows():
        if row['cc'] >= 0.1 and row['cc_pchange'] < 0:
            pass
        else:
            data_weakened.drop(index, inplace = True)
    save_data (data_weakened, file_out, dir_out)

def sort_NMDA (data, file_out, dir_out):
    data_NMDA = data.copy()
    for index, row in data_NMDA.iterrows():
        if row['cc'] < 0.1 and row['cc_nmda'] > 0.1 and row['cc_pchange'] > 0:
            pass
        else:
            data_NMDA.drop(index, inplace = True)
    save_data (data_NMDA, file_out, dir_out)

def sort_strength_created (data, file_out, dir_out):
    data_strength_created = data.copy()
    for index, row in data_strength_created.iterrows():
        if row['cc_nmda'] >= 0.1 and row['cc_pchange'] > 0:
            pass
        else:
            data_strength_created.drop(index, inplace = True)
    save_data (data_strength_created, file_out, dir_out)

#---------------------------#
#     Analyze the data      #
#---------------------------#


file_out = date_string + ' clean.csv'
# clean connections (remove negative % CC Values)
data = get_clean (data, 'vrec1', 'vrec2_nmda', 'vrec2', 'vrec2_nmda', file_out, dir_out)

# Get Rin
data_Rin = data.copy()      # Keep from overwriting original data
data_Rin = get_Rin(data_Rin, file_out, dir_out)
data_connect = data_Rin.copy()
data_connect = get_connections(data_connect, file_out, dir_out)

# Get cc, xfer, gj
file_out = exp_name + '_attempted.csv'
data_calc = data_connect.copy()
data_calc = get_cc(data_calc)
data_calc = get_xfer(data_calc)
data_calc = get_gj(data_calc)
data_calc = get_changes(data_calc)
data_calc = get_bin(data_calc)
save_data (data_calc, file_out, dir_out)

# Sorting
sort_connected (data_calc, exp_name + '_connected.csv', dir_out)
sort_strong (data_calc, exp_name + '_strong.csv', dir_out)
sort_weak (data_calc, exp_name + '_weak.csv', dir_out)
sort_strengthened (data_calc, exp_name + '_strengthened.csv', dir_out)
sort_strength_created (data_calc, exp_name + '_strengthened + Created.csv', dir_out)
sort_weakened (data_calc, exp_name + '_weakened.csv', dir_out)
sort_NMDA (data_calc, exp_name + '_nmda_created.csv', dir_out)

