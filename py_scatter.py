# If the axes plot weird, make sure all data in the column is numeric (no '--' denoting missing #data, etc)
# Replace repetitive code with lists and for statements
# Save Anova stats to file
# Reformat significance for paurwise and for main effects

import matplotlib.pyplot as plt  #For plotting/graphing
import pandas as pd              # For working with spreadsheets
import os                        # For working with files and directories 
import csv                       # For working with spreadsheets
import datetime                  # For timestamping your work
import scipy as sp               # For working with stats/advanced math functions
import numpy as np               # For working with advanced math functions
from scipy import stats          # For descriptive stats
import math                      # for isnan() function
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
from statsmodels.stats.libqsturng import psturng



# Input File variables
#dir_in = 'D:\\Dropbox\\Projects\\camk2a\wt_ko_cre_dunnet'
dir_in = 'C:\\Users\\haley\\Dropbox\\Projects\\camk2a\\Raw Data\\Filtered\\Calculated\\Attempted' # Directory on another computer
dir_stats = dir_in + '\\Stats'
dir_plots = dir_in + '\\Plots'
file_in = 'camk_attempted.csv'               # Group1 Raw data for plotting and stats
file_dunnett = 'camk_dunnett.csv'
groups = ['wt', 'ko', 'cre']                 # For naming files and deg of freedom calc

# Directory functions
os.chdir(dir_in)                            # Change the current directory to the one with all of your csv files
now = datetime.datetime.now()               # Gets the current date and time
date_string = now.strftime('%y%m%d %H.%M')  # formats the current date and time
#dir_out = dir_in + '\\Caculated ' + date_string       # Creates a timestamped name for a new folder for your output files


confidence = 0.95
repeated = 'n'                    # Repeated measures ANOVA? 'y' or 'n'      
exp_name = 'camk2a'
measures = ['cc', 'cc_nmda', 'cc_change', 'cc_pchange', 
            'gj', 'gj_nmda', 'gj_change', 'gj_pchange', 
            'rin1', 'rin1_nmda', 'rin_change', 'rin_pchange',  
            'vm1', 'vm_change']
posthoc = 'dunnett'                  # 'dunnett' or 'tukey'

#---------------------------------Function Definitions-------------------------------------#
# Load worksheets into pandas dataframes (df). 1 Worksheet = 1 dataframe
def get_frames (file_in, file_dunnett, groups):
    df_all = pd.read_csv (file_in)
    df_group1 = df_all[df_all.strain == groups[0]]
    df_group2 = df_all[df_all.strain == groups[1]]
    df_group3 = df_all[df_all.strain == groups[2]]
    df_dunnett = pd.read_csv(file_dunnett)
    df_groups = [df_group1, df_group2, df_group3, df_all, df_dunnett]
    
    
    return df_groups

# Save the dataframe to a csv
def save_csv (df_out, file_out, dir_stats, ignore):
    try:
        os.stat(dir_stats)
    except:
        os.mkdir(dir_stats)

    # Write data to file
    os.chdir(dir_stats)
    df_out.to_csv(file_out, index = ignore)
    

# Gets descriptive stats for one group
def get_desc (df_data, factor):
    parent_dir = os.getcwd()
    df_desc = df_data.copy()
    n = df_desc.groupby(factor).count()
    #print (n[1,0])
    avg = df_desc.groupby(factor).mean()
    sd = df_desc.groupby(factor).std()
    se = sd/np.sqrt(n.astype('int'))
    added = df_desc.groupby(factor).sum()
    minimum = df_desc.groupby(factor).min()
    maximum = df_desc.groupby(factor).max()
    quartile25 = df_desc.groupby(factor).quantile(q = 0.25, 
                axis = 0, numeric_only = True, interpolation = 'linear')
    quartile75 = df_desc.groupby(factor).quantile(q = 0.75, 
                axis = 0, numeric_only = True, interpolation = 'linear')
    median = df_desc.groupby(factor).median()
    conf = se * sp.stats.t._ppf((1+confidence)/2., n.astype('int')-1)
    conf_5 = avg - conf
    conf_95 = avg + conf
    
    n = n.reset_index()
    avg = avg.reset_index()
    sd = sd.reset_index()
    added = added.reset_index()
    minimum = minimum.reset_index()
    maximum = maximum.reset_index()
    quartile25 = quartile25.reset_index()
    quartile75 = quartile25.reset_index()
    median = median.reset_index()
    
    # Make new dataframe for stats
    df_return = pd.concat([n,avg,sd,se,added,minimum,maximum,quartile25,quartile75,median,conf_5,conf_95], sort = False, ignore_index=False)
    df_return = df_return.transpose()
    df_return.columns = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
      
    return df_return

def get_x_adj (stars):
    if stars == '****':
        x_adj = 0
    elif stars == '***':
        x_adj = 0
    elif stars == '**':
        x_adj = -0.1
    elif stars == '*':
        x_adj = -0.06
    elif stars == 'n.s.':
        x_adj = 0
    return x_adj


def get_stars (p):
    if p < 0.00001:
        stars = '****'
    elif p < 0.001:
        stars = '***'
    elif p < 0.01:
        stars = '**'
    elif p < 0.05:
        stars = '*'
    elif p >= 0.05:
        stars = 'n.s.'
    else:
        print (p)
    return stars

def get_anova (df_groups, measure, repeated, groups, dir_stats, exp_name, posthoc):
    
    group1_meas = df_groups[0][measure].dropna()
    group2_meas = df_groups[1][measure].dropna()
    group3_meas = df_groups[2][measure].dropna()
    
    # One-way Anova with strain as the independent factor
    f_value, p_anova = stats.f_oneway(group1_meas, group2_meas, group3_meas)
    stars_anova = get_stars(p_anova)

    # Degrees of Freedom 
    freedom1 = len(groups) - 1
    total_n = len(group1_meas) + len(group2_meas) + len(group3_meas)
    freedom2 = total_n - len(groups)
    freedom = [freedom1, freedom2]

    df_pairs = get_tukey (df_groups, measure, groups)
    
    for index, row in df_pairs.iterrows():
        stars = get_stars(row['p_value'])
        df_pairs.loc[index, 'significance'] = stars    
    
    anova_text =  'One-Way ANOVA: F(' + str(freedom[0]) + ',' + str(freedom[1]) + ') = ' + str(round(f_value,4)) + ', p = ' + str(round(p_anova, 4))
    
    if posthoc == 'tukey':
        df_pairs = get_tukey(df_groups, measure, groups)
    elif posthoc == 'dunnett':
        df_pairs = get_dunnett(df_groups, measure, groups)
    df_out = df_pairs.copy()
    df_out['anova'] = anova_text
    file_anova = exp_name + '_' + measure + '_' + posthoc + '.csv'
    save_csv (df_out, file_anova, dir_stats, False)
    
    return f_value, p_anova, freedom, stars_anova, df_pairs 

def get_tukey (df_groups, measure, groups):

    # Tukey posthoc analysis
    # See https://jpktd.blogspot.com/2013/03/multiple-comparison-and-tukey-hsd-or_25.html
    # And https://code.google.com/archive/p/qsturng-py/
    # And https://stackoverflow.com/questions/48200699/how-can-i-get-p-values-of-each-group-comparison-when-applying-the-tukey-s-hones
    # q, res_table, std_pairs, etc can be found from print(dir(result)) which will list all possible calculations
    df_tukey = df_groups[3][np.isfinite(df_groups[3][measure])]
    #print(df_tukey)
    mc = MultiComparison(df_tukey[measure], df_tukey['strain'])
    #result = pairwise_tukeyhsd(mc.data,mc.groups,0.05)
    result = mc.tukeyhsd()
    p = psturng(np.abs(result.meandiffs/result.std_pairs), len(result.groupsunique), result.df_total) 
    df_pairs = pd.DataFrame({'group1': [result._results_table[1][0], result._results_table[2][0], result._results_table[3][0]],
                             'group2': [result._results_table[1][1], result._results_table[2][1], result._results_table[3][1]],
                             'p_value': [np.around(p[0], 4), np.around(p[1], 4), np.around(p[2], 4)]})
    
    for index, row in df_pairs.iterrows():
        stars = get_stars(row['p_value'])
        df_pairs.loc[index, 'significance'] = stars   

    return df_pairs 

def get_dunnett (df_groups, measure, groups):
    
    # Dunnet posthoc analysis
    # Performed in camk2_dunnett.r and imported here
    # Manually added "strain" as name of first column because it was faster than programming it 
    
    df_dun = df_groups[4]
    for index, row in df_dun.iterrows():
        if row['measure'] == measure and row['strain'] == 'wt-ko':
            group1_wt = 'wt'
            p_val_wt = row['pval']
        if row['measure'] == measure and row['strain'] == 'cre-ko':
            group1_cre = 'cre'
            p_val_cre = row['pval']
    df_pairs = pd.DataFrame({'group1': [group1_wt, group1_cre],
                             'group2': ['ko', 'ko'],
                             'p_value': [np.around(p_val_wt,4), np.around(p_val_cre,4)]})
    

    
    for index, row in df_pairs.iterrows():
        stars = get_stars(row['p_value'])
        df_pairs.loc[index, 'significance'] = stars    
    
    return df_pairs 
    

def get_measure (df_groups, df_stats1, df_stats2, df_stats3, measure, repeated, exp_name, dir_stats, dir_plots, groups, posthoc):
    df_stats1 = df_stats1.transpose()
    df_stats2 = df_stats2.transpose()
    df_stats3 = df_stats3.transpose()
    
    y_column_name = measure                             # Header name for your y-axis column in your raw data
    save_fig = exp_name + '_' + measure +'.png'         # Name of the output figure. Jpg and TIFF cause trouble. Stick to png or pdf
    
    # Determines axis labels and data to display
    if measure == 'cc':
        col_num = 17
        y_hline = 0.1   # Add a horizontal line at a specific y value
        y_label = 'Coupling Coefficient (%) \n '# Label for the y axis
    elif measure == 'cc_nmda':
        col_num = 18
        y_hline = 0
        y_label = 'Coupling Coefficient (%) \n with NMDA \n '
    elif measure == 'cc_change':
        col_num = 23
        y_hline = 0
        y_label = '\u0394 Coupling Coefficient \n with NMDA \n '
    elif measure == 'cc_pchange':
        col_num = 24
        y_hline = 0
        y_label = '% \u0394 Coupling Coefficient (%) \n with NMDA \n'
    elif measure == 'gj':
        col_num = 21
        y_hline = 0
        y_label = 'Junctional Conductance (pS) \n'
    elif measure == 'gj_nmda':
        col_num = 22
        y_hline = 0
        y_label = 'Junctional Conductance (pS) \n with NMDA \n'
    elif measure == 'gj_change':
        col_num = 27
        y_hline = 0
        y_label = '\u0394 Junctional Conductance (pS) \n with NMDA \n'
    elif measure == 'gj_pchange':
        col_num = 28
        y_hline = 0
        y_label = '% \u0394 Junctional Conductance \n with NMDA \n'
    elif measure == 'rin1':
        col_num = 9
        y_hline = 0
        y_label = 'Input Resistance (M\u03A9) \n'
    elif measure == 'rin1_nmda':
        col_num = 11
        y_hline = 0
        y_label = 'Input Resistance (M\u03A9) \n with NMDA \n'
    elif measure == 'rin_change':
        col_num = 25
        y_hline = 0
        y_label = '\u0394 Input Resistance (M\u03A9)  \n with NMDA \n'
    elif measure == 'rin_pchange':
        col_num = 26
        y_hline = 0
        y_label = '% \u0394 Input Resistance (M\u03A9) \n with NMDA \n'
    elif measure == 'vm1':
        col_num = 13
        y_hline = -40
        y_label = 'Resting Membrane Potential (mV) \n'
    elif measure == 'vm_change':
        col_num = 29
        y_hline = 0
        y_label = '\u0394 Resting Membrane Potential (mV) \n with NMDA \n'
            
    
    # Assemble a descriptive stats spreadsheet for 1 type of measurement
    df_stats = pd.DataFrame({'strain':[df_stats1.iloc[1,0],df_stats2.iloc[1,0],df_stats3.iloc[1,0]],
                            'n':[df_stats1.iloc[0,col_num], df_stats2.iloc[0,col_num], df_stats3.iloc[0,col_num]], 
                            'mean':[df_stats1.iloc[1,col_num], df_stats2.iloc[1,col_num], df_stats3.iloc[1,col_num]],
                            'se':[df_stats1.iloc[3,col_num], df_stats2.iloc[3,col_num], df_stats3.iloc[3,col_num]],
                            'sd':[df_stats1.iloc[2,col_num], df_stats2.iloc[2,col_num], df_stats3.iloc[2,col_num]],
                            'sum':[df_stats1.iloc[4,col_num], df_stats2.iloc[4,col_num], df_stats3.iloc[4,col_num]],
                            'min':[df_stats1.iloc[5,col_num], df_stats2.iloc[5,col_num], df_stats3.iloc[5,col_num]],
                            'max':[df_stats1.iloc[6,col_num], df_stats2.iloc[6,col_num], df_stats3.iloc[6,col_num]],
                            '25_quartile':[df_stats1.iloc[7,col_num], df_stats2.iloc[7,col_num], df_stats3.iloc[7,col_num]],
                            '75_quartile':[df_stats1.iloc[8,col_num], df_stats2.iloc[8,col_num], df_stats3.iloc[8,col_num]],
                            'median':[df_stats1.iloc[9,col_num], df_stats2.iloc[9,col_num], df_stats3.iloc[9,col_num]],
                            '5_ci':[df_stats1.iloc[10,col_num], df_stats2.iloc[10,col_num], df_stats3.iloc[10,col_num]],
                            '95_ci':[df_stats1.iloc[11,col_num], df_stats2.iloc[11,col_num], df_stats3.iloc[11,col_num]]})
    measure_out = exp_name + '_' + measure + '_analyzed.csv'
    save_csv (df_stats, measure_out, dir_stats, False) # False = Keep row labels
    
    # Run stats on two groups using independent t-test
    f_value, p_values, freedom, stars_anova, df_pairs = get_anova (df_groups, measure, repeated, groups, dir_stats, exp_name, posthoc)
        
    # Plot the data with the significance
    get_plot (df_groups, df_stats, y_column_name, save_fig, y_label, stars_anova, df_pairs, y_hline, dir_plots, measure)
    
    return col_num

def get_plot (df_groups, df_stats, y_column_name, save_fig, y_label, stars_anova, df_pairs, y_hline, dir_plots, measure):
    
     # Graph variables
    x_column_name = 'strain'            # Header name for your x-axis column in your raw data
    mean_column_name = 'mean'           # Header name for you mean colum in df_stats
    se_column_name = 'se'               # Header name for your standard error column 
    x_ticklabel1 = 'WT'                 # Replaces tick labels on the x-axis for the first scatter group  
    x_ticklabel2 = 'CaMKII\u03B1 KO'                 # Replaces tick labels on the x-axis for the second scatter group 
    x_ticklabel3 = 'ptf1a-Cre'          # Replaces tick labels on the x-axis for the third scatter group 
    x_label = ''                        # label for the x axis
    x_fontsize = 14                     # Font size of the x label
    y_fontsize = 14                     # Font size of the y label
    tick_fontsize = 12                  # Font size of the x and y axis tick labels 
    plot_colors = ['black', 'dimgray', 'teal', 'blueviolet']   # Sequence of colors for the plots in your data: Mean, scatter 1, scatter2
    plot_markers = ['_', 's','o','v']   # Marker shape list: solid line (mean), square (scatter1), circle (scatter2)
    scatter_marker_size = 4             # Marker size for the scatter plots
    mean_marker_size = 50               # Increases the size of the mean marker to be a longer solid line
    x_min = -0.5                        # x axis minimum 
    x_max = 2.5                         # x axis maxumim
    y_hline_marker = '--'               # Horizontal line will be dashed 
    y_hline_color = 'silver'            # Horizontal line will be a light gray color
    y_hline_width = 1                   # Thickness of the horizontal line 
    x_ticklabel_angle = 45              # Rotates angle of the the x-axis tick labels  
    plot_dpi = 80                       # Resolution of the image displayed in the jupyter notebook
    save_dpi = 300                      # Resolution of the saved image (publication quality)

    # Autoscale: get range of both groups
    y_min = df_stats['min'].min()
    y_max =  df_stats['max'].max()
    
    # Autoscale: Pad the y-axis min and max by 25% of the range
    y_range = abs(y_min) + y_max 
    scale_pad = 0.25 * (y_range)
    y_min = y_min - scale_pad
    y_max = y_max + scale_pad
    
    # Assign data to x and y variables
    x_ticklabel1 = x_ticklabel1 + '\n(' + str(int(df_stats.iloc[0,1])) + ')'
    x_ticklabel2 = x_ticklabel2 + '\n(' + str(int(df_stats.iloc[1,1])) + ')'
    x_ticklabel3 = x_ticklabel3 + '\n(' + str(int(df_stats.iloc[2,1])) + ')'
    x_data1 = df_groups[0][x_column_name]  # Imports the x axis data for the first scatter group 
    y_data1 = df_groups[0][y_column_name]  # Imports the y axis data for the first scatter group
    x_data2 = df_groups[1][x_column_name]  # Imports the x axis data for the second scatter group
    y_data2 = df_groups[1][y_column_name]  # Imports the y axis data for the second scatter group
    x_data3 = df_groups[2][x_column_name]  # Imports the x axis data for the second scatter group
    y_data3 = df_groups[2][y_column_name]  # Imports the y axis data for the second scatter group
    x_means = df_stats[x_column_name]  # Imports the x axis data for the means 
    y_means = df_stats[mean_column_name]  # Imports the y axis data for the means
    y_err = df_stats[se_column_name]    # Imports the standard error data for the means
      
    # Generate the Figure
    fig, ax = plt.subplots(figsize=(4, 3), dpi=plot_dpi) # subplots rather than plt will be used to fine-tune the output

    # Format axes aesthetics
    ax.tick_params(labelsize = tick_fontsize)             # Increase the size of tick labels
    ax.set_xlabel(x_label, fontsize = x_fontsize)         # Increase the size of X-axis label
    ax.set_ylabel(y_label, fontsize = y_fontsize)         # Increase the size of the y-axis label
    ax.set_prop_cycle(color = plot_colors, marker = plot_markers) # Set the sequence of colors and marker shapes

    # Remove the right and top bounding lines
    ax.spines['top'].set_visible(False) 
    ax.spines['right'].set_visible(False)

    # Plot the data
    ax.plot(x_means, y_means, linestyle = 'none', markersize = mean_marker_size)    # plot first mean
    ax.errorbar(x_means, y_means, yerr = y_err, linestyle = 'none', elinewidth = 2, capsize = 5, color = 'black' ) # plot error bars
    ax.plot(x_data1, y_data1, linestyle = 'none', markersize = scatter_marker_size) # plot first scatter dataset
    ax.plot(x_data2, y_data2, linestyle = 'none', markersize = scatter_marker_size) # plot second scatter dataset
    ax.plot(x_data3, y_data3, linestyle = 'none', markersize = scatter_marker_size) # plot second scatter dataset

    # Set the axes 
    ax.axis ([x_min, x_max, y_min, y_max]) # Set axes limits [xmin, xmax, ymin, ymax]

    # Add ANOVA significance to top left of graph
    if stars_anova != 'n.s.':
        sig = stars_anova + ' ' + 'Genotype '
        ax.text(-0.4, y_max, sig, size = 10, color = 'black')

    # Add tukey significance
    # Compared to wt ***   compared to cre ###
    df_temp = df_groups[1].groupby('strain')
    y_minmax = df_temp[measure].max()
    x_range = x_max - x_min
    x_star = x_range/3
    for index, row in df_pairs.iterrows():   
        if str(row['significance']) != 'n.s.' and str(row['group1']) == 'cre' or str(row['group2']) == 'wt':  
            x_adj = get_x_adj (row['significance'])
            y_adj = y_minmax + abs(0.2* y_minmax)
            row['significance'] = row['significance'].replace('*', '#')
            ax.text(x_star +  x_adj, y_adj, row['significance'], size = 12, color = 'black', fontweight = 'bold')
        if str(row['significance']) != 'n.s.' and str(row['group1']) == 'wt' or str(row['group2']) == 'wt':
            x_adj = get_x_adj (row['significance'])
            y_adj = y_minmax + abs(0.5* y_minmax)
            ax.text(x_star + x_adj, y_adj, row['significance'], size = 16, color = 'black', fontweight = 'bold')
                    

    # Draw the horizontal line (line y-value is determined in get_measure())
    ax.axhline(y_hline, 0, 1 , color = y_hline_color, lw = y_hline_width, linestyle = y_hline_marker) # Add a horizontal line
      
    # Change the labels of the x axis (optional)
    ax.set_xticklabels ([x_ticklabel1, x_ticklabel2, x_ticklabel3], rotation = x_ticklabel_angle)

    # Plot the graph in the Jupyter Notebook
    #plt.show()

    # Save the figure to a new subfolder if the subfolder does not already exist
    # Checks to see if the output folder exists
    try:                         
        os.stat(dir_plots)
    except:
        os.mkdir(dir_plots)
    
    # Change the current directory to the output directory
    os.chdir(dir_plots)

    # Write the figure to file
    fig.savefig(save_fig, dpi=save_dpi, facecolor='w', edgecolor='w',  
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        metadata=None)

    # Reset the working directory to your input directory
    os.chdir(dir_in)



#----------------------------Execution--------------------------------------

df_groups = get_frames (file_in, file_dunnett, groups)

# Gets Descriptive Stats 
df_stats1 = get_desc (df_groups[0], 'strain')
df_stats2 = get_desc (df_groups[1], 'strain')
df_stats3 = get_desc (df_groups[2], 'strain')
save_group1 = groups[0] + '_desc'+ '.csv'
save_group2 = groups[1] + '_desc'+ '.csv'
save_group3 = groups[2] + '_desc'+ '.csv'
save_csv (df_stats1, save_group1, dir_stats, True) # Ignore index = true (no row labels)
save_csv (df_stats2, save_group2, dir_stats, True) # Ignore indec - True (no row labels)
save_csv (df_stats3, save_group3, dir_stats, True) # Ignore indec - True (no row labels)


# Gets plots for each measurement
for measure in measures:
    col_num = get_measure(df_groups, df_stats1, df_stats2, df_stats3, 
                          measure, repeated, exp_name, dir_stats, 
                          dir_plots, groups, posthoc)
    

