import astropy
import gzip
import glob
import numpy as np
import math
import xlrd
import re 
import shutil
import matplotlib.pyplot as plt
import pickle

from astropy.time import Time,  TimeJD
from datetime import datetime
from copy import deepcopy as cp
from scipy import stats
from photutils.aperture import aperture_photometry
from astropy.visualization import simple_norm
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.stats import sigma_clipped_stats
from photutils import find_peaks
from photutils.detection import DAOStarFinder
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from scipy import interpolate
from scipy.interpolate import interp1d
from reduct_funcs import funcs_utils
        
######## An Experiment

def EECep_light_curve_based_pol(target_data, MJD_L_cutoff , MJD_cutoff, point_MJD , lc, verb_t, mark_t ,txt_arg, true_pa, title_t):
    print("Reading:", 'ee_cep_2014.txt')
    file1 = open('./EE_Cep_light_curve/ee_cep_2014.txt', 'r')
    #file1 = open('./EE_Cep_light_curve/ee_cep_small_our_obs.txt', 'r')
    ee_lines = file1.readlines()
    
    if(verb_t):
        print("Plotting:",txt_arg)
        print("Len all:", len(target_data[0]), len(target_data[1]), len(target_data[2]))
        for h in range(0, len(target_data[0])):
            print(h, target_data[0][h] , target_data[1][h] , target_data[2][h].value)
    
    B_arr = []
    V_arr = []
    R_arr = []
    I_arr = []

    B_JD = []
    V_JD = []
    R_JD = []
    I_JD = []

    for things in ee_lines:
        if ',B,' in things:
            B_arr.append( float(things.split(",")[1]))
            B_JD.append(  float(things.split(",")[0]))
        elif ',V,' in things:
            V_arr.append( float(things.split(",")[1]))
            V_JD.append(  float(things.split(",")[0]))
        elif ',R,' in things:
            R_arr.append( float(things.split(",")[1]))
            R_JD.append(  float(things.split(",")[0]))
        elif ',I,' in things:
            I_arr.append( float(things.split(",")[1]))
            I_JD.append(  float(things.split(",")[0]))

    t_b = Time(B_JD , format='jd', scale='utc')
    t_v = Time(V_JD , format='jd', scale='utc')
    t_r = Time(R_JD , format='jd', scale='utc')
    t_i = Time(I_JD , format='jd', scale='utc')
        
    color_arr = [B_arr, V_arr, R_arr, I_arr]
    mjd_arr_k = [t_b, t_v, t_r, t_i, target_data[2]]
    
    #c_arr=['purple', 'red', 'green', 'blue']
    #c_labs=['I', 'R', 'V', 'B', 'pol_d']
    
    #originals
    c_arr = ['blue',  'green',  'red', 'purple']
    c_labs = ['B', 'V', 'R', 'I', 'pol_d']
    
    #c_labs = ['I', 'R', 'V', 'B', 'pol_d']    
    index_cutt = []
    index_L_cutt = []
    
    for jk in range(0, len(mjd_arr_k)):
        idx = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cutoff).argmin()
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_L_cutoff).argmin()
        index_cutt.append(int(idx))
        index_L_cutt.append(int(idxL))      

    if(verb_t):
        inspect_val = 500
        print("\n")
        print("The light curve is long:", len(mjd_arr_k[0]), len(mjd_arr_k[1]), len(mjd_arr_k[2]), len(mjd_arr_k[3]), "And of varying length")
        print("Starting at:", mjd_arr_k[0][0],  mjd_arr_k[1][0], mjd_arr_k[2][0], mjd_arr_k[3][0])
        print("Ending at:", mjd_arr_k[0][-1],  mjd_arr_k[1][-1], mjd_arr_k[2][-1], mjd_arr_k[3][-1])
        print("At those values:", (mjd_arr_k[0][inspect_val]), (mjd_arr_k[1][inspect_val]), (mjd_arr_k[2][inspect_val]), (mjd_arr_k[3][inspect_val]))
    
    fig, ax1 = plt.subplots(figsize=(36, 12))
    
    for ka in range(0, len(mjd_arr_k)-1): #this is a recursive call. Just recursive call the third plot
        ax1.scatter(mjd_arr_k[ka][index_L_cutt[ka]:index_cutt[ka]].mjd,  color_arr[ka][index_L_cutt[ka]:index_cutt[ka]], alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])
       
    fig.gca().invert_yaxis()
    
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)

    plt.xlabel('Time, (MJD)', fontsize=32)
    plt.ylabel('magnitude, (m)', fontsize=32)
    
    if(txt_arg == 'PD'):
        if(lc):
            if(title_t):
                ax1.set_title("EE Cep light curve and Polarization Degree ("+txt_arg+")", fontsize=32)
        else:
            if(title_t):
                ax1.set_title("EE Cep light curve ", fontsize=32)
    elif(txt_arg == 'PA'):
        if(lc):
            if(true_pa):
                if(title_t):
                    ax1.set_title("EE Cep light curve and true Position Angle ("+txt_arg+")", fontsize=32)
            else:
                if(title_t):
                    ax1.set_title("EE Cep light curve and Position Angle ("+txt_arg+")", fontsize=32)
        else:
            ax1.set_title("EE Cep light curve ", fontsize=32)
        
    ax1.grid()
    ax1.plot()
    
    if(lc):
        ax2 = ax1.twinx()

        if(txt_arg == 'PD'):              
            if(MJD_cutoff >= 59354.9925 ):#include
                markers, caps, bars = ax2.errorbar(target_data[2][index_L_cutt[4]:index_cutt[4]+1].mjd, target_data[0][index_L_cutt[4]:index_cutt[4]+1], yerr=target_data[1][index_L_cutt[4]:index_cutt[4]+1], xerr =[0]*len(target_data[0][index_L_cutt[4]:index_cutt[4]+1]), color='black' ,fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PD')
            else:
                markers, caps, bars = ax2.errorbar(target_data[2][index_L_cutt[4]:index_cutt[4]].mjd, target_data[0][index_L_cutt[4]:index_cutt[4]], yerr=target_data[1][index_L_cutt[4]:index_cutt[4]], xerr =[0]*len(target_data[0][index_L_cutt[4]:index_cutt[4]]), color='black', fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PD')
                
        elif(txt_arg == 'PA'):
            if( MJD_cutoff >= 59354.9925 ):#include
                markers, caps, bars = ax2.errorbar(target_data[2][index_L_cutt[4]:index_cutt[4]+1].mjd, target_data[0][index_L_cutt[4]:index_cutt[4]+1], yerr=target_data[1][index_L_cutt[4]:index_cutt[4]+1], xerr =[0]*len(target_data[0][index_L_cutt[4]:index_cutt[4]+1]), color='black', fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PA')
            else:
                markers, caps, bars = ax2.errorbar(target_data[2][index_L_cutt[4]:index_cutt[4]].mjd, target_data[0][index_L_cutt[4]:index_cutt[4]], yerr=target_data[1][index_L_cutt[4]:index_cutt[4]], xerr =[0]*len(target_data[0][index_L_cutt[4]:index_cutt[4]]), color='black', fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PA')

        if(txt_arg == 'PD'):
            ax2.set_ylabel('PD, (%)', fontsize=32)
            plt.yticks(fontsize=28)
        elif(txt_arg == 'PA'):
            if(true_pa):
                ax2.set_ylabel('true PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
            else:
                ax2.set_ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
            plt.yticks(fontsize=28)

        if(point_MJD):
            for hi in range(0, len(target_data[0])):
                if(hi%2 == 0):
                    plt.text(target_data[2][hi].mjd, target_data[0][hi],target_data[2][hi].mjd, fontsize=12, alpha=0.69, rotation=45    )

        [bar.set_alpha(0.2) for bar in bars]
        [cap.set_alpha(0.72) for cap in caps]

    lgnd = fig.legend(loc="upper right", fontsize=36, borderaxespad=3.2)#, borderaxespad=0.86)
    
    lgnd.legendHandles[0]._sizes = [999]
    lgnd.legendHandles[1]._sizes = [999]
    lgnd.legendHandles[2]._sizes = [999]
    lgnd.legendHandles[3]._sizes = [999]
    lgnd.legendHandles[4]._sizes = [999]
    
    fig.tight_layout()
    plt.show()

def EECep_stacked_based_pol(target_data_PD, target_data_PA, MJD_L_cutoff , MJD_cutoff, point_MJD , lc, verb_t, mark_td, mark_t, true_pa, title_t):
    """
    A function that plot the EE Cep light curve
    """
    file1 = open('./EE_Cep_light_curve/ee_cep_big_dump.txt', 'r')
    ee_lines = file1.readlines()
    
    B_arr = []
    V_arr = []
    R_arr = []
    I_arr = []

    B_JD = []
    V_JD = []
    R_JD = []
    I_JD = []

    for things in ee_lines:
        if ',B,' in things:
            B_arr.append( float(things.split(",")[1]))
            B_JD.append(  float(things.split(",")[0]))
        elif ',V,' in things:
            V_arr.append( float(things.split(",")[1]))
            V_JD.append(  float(things.split(",")[0]))
        elif ',R,' in things:
            R_arr.append( float(things.split(",")[1]))
            R_JD.append(  float(things.split(",")[0]))
        elif ',I,' in things:
            I_arr.append( float(things.split(",")[1]))
            I_JD.append(  float(things.split(",")[0]))

    t_b = Time(B_JD , format='jd', scale='utc')
    t_v = Time(V_JD , format='jd', scale='utc')
    t_r = Time(R_JD , format='jd', scale='utc')
    t_i = Time(I_JD , format='jd', scale='utc')
        
    color_arr = [B_arr, V_arr, R_arr, I_arr]
    mjd_arr_k = [t_b, t_v, t_r, t_i, target_data_PD[2]] #Which is just the time array
        
    #originals
    c_arr = ['blue',  'green',  'red', 'purple']
    c_labs = ['B', 'V', 'R', 'I', 'pol_d']
        
    index_cutt = []
    index_L_cutt = []
    
    for jk in range(0, len(mjd_arr_k)):
        idx = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cutoff).argmin()
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_L_cutoff).argmin()
        index_cutt.append(int(idx))
        index_L_cutt.append(int(idxL))
        
    #This makes things...

    x = np.linspace(0, 2 * np.pi, 400)
    y = np.sin(x ** 2)
    
    fig = plt.figure(figsize=(36, 12))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True, sharey=False)
    
    if(true_pa):
        if(title_t):
            fig.suptitle('Polarization degree (PD), true Position angle (PA) and EE CEP light curve', fontsize=32)    
    else:
        if(title_t):
            fig.suptitle('Polarization degree (PD), Position angle (PA) and EE CEP light curve', fontsize=32)
    
    #We should include last plot also 
    if(MJD_cutoff >= 59354.9925 ):#include
        markers, caps, bars = axs[0].errorbar(target_data_PD[2][index_L_cutt[4]:index_cutt[4]+1].mjd, target_data_PD[0][index_L_cutt[4]:index_cutt[4]+1], yerr=target_data_PD[1][index_L_cutt[4]:index_cutt[4]+1], xerr =[0]*len(target_data_PD[0][index_L_cutt[4]:index_cutt[4]+1]), color='black' ,fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PD')
        
        [bar.set_alpha(0.2) for bar in bars]
        [cap.set_alpha(0.72) for cap in caps]
        
            #axs[1].plot(target_data_PA[2][index_L_cutt[4]:index_cutt[4]+1].mjd, target_data_PA[0])
        markers, caps, bars = axs[1].errorbar(target_data_PA[2][index_L_cutt[4]:index_cutt[4]+1].mjd, target_data_PA[0][index_L_cutt[4]:index_cutt[4]+1], yerr=target_data_PA[1][index_L_cutt[4]:index_cutt[4]+1], xerr =[0]*len(target_data_PA[0][index_L_cutt[4]:index_cutt[4]+1]), color='black' ,fmt=mark_td, markersize=16, capsize=10, capthick=5, label='PA')
        axs[0].grid()
        axs[1].grid()
        
        [bar.set_alpha(0.2) for bar in bars]
        [cap.set_alpha(0.72) for cap in caps]
        
        axs[0].set_ylabel('PD, (%)', fontsize=32)
        axs[0].tick_params(axis="y", labelsize=28)
        
        if(true_pa):
            axs[1].set_ylabel('true PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
        else:
            axs[1].set_ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
            
        axs[1].tick_params(axis="y", labelsize=28)
        
    else:       
        markers, caps, bars = axs[0].errorbar(target_data_PD[2][index_L_cutt[4]:index_cutt[4]].mjd, target_data_PD[0][index_L_cutt[4]:index_cutt[4]], yerr=target_data_PD[1][index_L_cutt[4]:index_cutt[4]], xerr =[0]*len(target_data_PD[0][index_L_cutt[4]:index_cutt[4]]), color='black', fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PD')
            
        markers, caps, bars = axs[1].errorbar(target_data_PA[2][index_L_cutt[4]:index_cutt[4]].mjd, target_data_PA[0][index_L_cutt[4]:index_cutt[4]], yerr=target_data_PA[1][index_L_cutt[4]:index_cutt[4]], xerr =[0]*len(target_data_PA[0][index_L_cutt[4]:index_cutt[4]]), color='black', fmt=mark_td, markersize=16, capsize=10, capthick=5, label='PA')
        axs[0].grid()
        axs[1].grid()       
         
        axs[0].set_ylabel('PD, (%)', fontsize=32)
        axs[0].tick_params(axis="y", labelsize=28)
        axs[1].set_ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
        #axs[1].set_yticklabels(fontsize=28) #Disabled for debugging
        axs[1].tick_params(axis="y", labelsize=28)

    for ka in range(0, len(mjd_arr_k)-1): #this is a recursive call. Just recursive call the third plot
        axs[2].scatter(mjd_arr_k[ka][index_L_cutt[ka]:index_cutt[ka]].mjd,  color_arr[ka][index_L_cutt[ka]:index_cutt[ka]], alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])

    axs[2].grid()
    axs[2].set_ylabel('magnitude, (m)', fontsize=32)
    axs[2].invert_yaxis()
    axs[2].set_xlabel('Time, (MJD)', fontsize=32)

    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)

    # Hide x labels and tick labels for all but bottom plot.
    for ax in axs:
        ax.label_outer()

    lgnd = fig.legend(loc="upper right", fontsize=36, borderaxespad=1.75)#, borderaxespad=0.86)
    
    lgnd.legendHandles[0]._sizes = [999]
    lgnd.legendHandles[1]._sizes = [999]
    lgnd.legendHandles[2]._sizes = [999]
    lgnd.legendHandles[3]._sizes = [999]
    lgnd.legendHandles[4]._sizes = [999]
    lgnd.legendHandles[5]._sizes = [999]
        
    fig.tight_layout()
    plt.show()

def EECep_stacked_based(target_data_PD, target_data_PA, MJD_cuts, point_MJD , lc, verb_t, mark_td, mark_t, true_pa, title_t):
    """
    #Target_data PD = list of polarization degrees
    #Target_data PA = list of polarization angles
    #MJD_L_cutoffs
    #MJD_cutoffs
    """
    #print("Openning 2014")
    file1 = open('./EE_Cep_light_curve/ee_cep_2014.txt', 'r')
    ee_lines = file1.readlines()
    
    B_arr = []
    V_arr = []
    R_arr = []
    I_arr = []

    B_JD = []
    V_JD = []
    R_JD = []
    I_JD = []

    for things in ee_lines:
        if ',B,' in things:
            B_arr.append( float(things.split(",")[1]))
            B_JD.append(  float(things.split(",")[0]))
        elif ',V,' in things:
            V_arr.append( float(things.split(",")[1]))
            V_JD.append(  float(things.split(",")[0]))
        elif ',R,' in things:
            R_arr.append( float(things.split(",")[1]))
            R_JD.append(  float(things.split(",")[0]))
        elif ',I,' in things:
            I_arr.append( float(things.split(",")[1]))
            I_JD.append(  float(things.split(",")[0]))

    t_b = Time(B_JD , format='jd', scale='utc')
    t_v = Time(V_JD , format='jd', scale='utc')
    t_r = Time(R_JD , format='jd', scale='utc')
    t_i = Time(I_JD , format='jd', scale='utc')
        
    color_arr = [B_arr, V_arr, R_arr, I_arr]
    mjd_arr_k = [t_b, t_v, t_r, t_i, target_data_PD[2]]
        
    #originals
    c_arr = ['blue',  'green',  'red', 'purple']
    c_labs = ['B', 'V', 'R', 'I', 'pol_d']
        
    index_cutt = []
    index_L_cutt = []
    
    for jk in range(0, len(mjd_arr_k)):
        idx = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts[1]).argmin() #Left_side
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts[0]).argmin() #Right_Side
        index_cutt.append(int(idx))
        index_L_cutt.append(int(idxL))
        
    #This makes things...

    x = np.linspace(0, 2 * np.pi, 400)
    y = np.sin(x ** 2)
    
    fig = plt.figure(figsize=(36, 12))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True, sharey=False)
    
    if(true_pa):
        if(title_t):
            fig.suptitle('Polarization degree (PD), true Position angle (PA) and EE CEP light curve', fontsize=32)    
    else:
        if(title_t):
            fig.suptitle('Polarization degree (PD), Position angle (PA) and EE CEP light curve', fontsize=32)
    
    #We should include last plot also 
    if(MJD_cuts[1] >= 59354.9925 ):#include. Why is there
                                   #It is for the MJD_cuts_L
                                   #If the right cut is more than 59354.9925 whci is on the high end of tings
        markers, caps, bars = axs[0].errorbar(target_data_PD[2][index_L_cutt[4]:index_cutt[4]+1].mjd, target_data_PD[0][index_L_cutt[4]:index_cutt[4]+1], yerr=target_data_PD[1][index_L_cutt[4]:index_cutt[4]+1], xerr =[0]*len(target_data_PD[0][index_L_cutt[4]:index_cutt[4]+1]), color='black' ,fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PD')
        
        [bar.set_alpha(0.2) for bar in bars]
        [cap.set_alpha(0.72) for cap in caps]
        
            #axs[1].plot(target_data_PA[2][index_L_cutt[4]:index_cutt[4]+1].mjd, target_data_PA[0])
        markers, caps, bars = axs[1].errorbar(target_data_PA[2][index_L_cutt[4]:index_cutt[4]+1].mjd, target_data_PA[0][index_L_cutt[4]:index_cutt[4]+1], yerr=target_data_PA[1][index_L_cutt[4]:index_cutt[4]+1], xerr =[0]*len(target_data_PA[0][index_L_cutt[4]:index_cutt[4]+1]), color='black' ,fmt=mark_td, markersize=16, capsize=10, capthick=5, label='PA')
        axs[0].grid()
        axs[1].grid()
        
        [bar.set_alpha(0.2) for bar in bars]
        [cap.set_alpha(0.72) for cap in caps]
        
        axs[0].set_ylabel('PD, (%)', fontsize=32)
        axs[0].tick_params(axis="y", labelsize=28)
        
        if(true_pa):
            axs[1].set_ylabel('true PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
        else:
            axs[1].set_ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
            
        axs[1].tick_params(axis="y", labelsize=28)
        
    else:       
        markers, caps, bars = axs[0].errorbar(target_data_PD[2][index_L_cutt[4]:index_cutt[4]].mjd, target_data_PD[0][index_L_cutt[4]:index_cutt[4]], yerr=target_data_PD[1][index_L_cutt[4]:index_cutt[4]], xerr =[0]*len(target_data_PD[0][index_L_cutt[4]:index_cutt[4]]), color='black', fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PD')
            
        markers, caps, bars = axs[1].errorbar(target_data_PA[2][index_L_cutt[4]:index_cutt[4]].mjd, target_data_PA[0][index_L_cutt[4]:index_cutt[4]], yerr=target_data_PA[1][index_L_cutt[4]:index_cutt[4]], xerr =[0]*len(target_data_PA[0][index_L_cutt[4]:index_cutt[4]]), color='black', fmt=mark_td, markersize=16, capsize=10, capthick=5, label='PA')
        axs[0].grid()
        axs[1].grid()       
         
        axs[0].set_ylabel('PD, (%)', fontsize=32)
        axs[0].tick_params(axis="y", labelsize=28)
        axs[1].set_ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
        
        #axs[1].set_yticklabels(fontsize=28)        #Was this disabled for safety?
        #axs[1].tick_params(axis="y", labelsize=28)
        
        axs[1].tick_params(axis="y", labelsize=28)

    for ka in range(0, len(mjd_arr_k)-1): #this is a recursive call. Just recursive call the third plot
        axs[2].scatter(mjd_arr_k[ka][index_L_cutt[ka]:index_cutt[ka]].mjd,  color_arr[ka][index_L_cutt[ka]:index_cutt[ka]], alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])

    axs[2].grid()
    axs[2].set_ylabel('magnitude, (m)', fontsize=32)
    axs[2].invert_yaxis()
    axs[2].set_xlabel('Time, (MJD)', fontsize=32)

    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)

    # Hide x labels and tick labels for all but bottom plot.
    for ax in axs:
        ax.label_outer()

    lgnd = fig.legend(loc="upper right", fontsize=36, borderaxespad=1.75)#, borderaxespad=0.86)
    
    lgnd.legendHandles[0]._sizes = [999]
    lgnd.legendHandles[1]._sizes = [999]
    lgnd.legendHandles[2]._sizes = [999]
    lgnd.legendHandles[3]._sizes = [999]
    lgnd.legendHandles[4]._sizes = [999]
    lgnd.legendHandles[5]._sizes = [999]
        
    fig.tight_layout()
    plt.show()
    
def destroy_outliers(samp_MJD, cleaned_light_curve, line_a, line_b, line_c, line_d, line_e, line_f):
    destroy_list_a = [i for i in range(len(cleaned_light_curve)) if cleaned_light_curve[i] <= line_a] #
    destroy_list_b = [i for i in range(len(cleaned_light_curve)) if cleaned_light_curve[i] >= line_b and samp_MJD[i] <= line_c ]
    destroy_list_c = [i for i in range(len(cleaned_light_curve)) if cleaned_light_curve[i] >= line_b and samp_MJD[i] >= line_d and samp_MJD[i] <= line_e]
    destroy_list_d = [i for i in range(len(cleaned_light_curve)) if cleaned_light_curve[i] >= line_b and samp_MJD[i] >= line_f]
    
    f_destroy = destroy_list_a + destroy_list_b + destroy_list_c + destroy_list_d
    me_mod=0
    
    print("Length Original:", len(cleaned_light_curve),len(samp_MJD))
        
    for o in range(0, len(f_destroy)):
        del(cleaned_light_curve[f_destroy[o] - me_mod ])
        samp_MJD = np.delete(samp_MJD, f_destroy[o] - me_mod)
        me_mod+=1
    
    print("Removing:", me_mod, len(f_destroy))
    print("Length New:", len(cleaned_light_curve),len(samp_MJD) )
    
    return samp_MJD, cleaned_light_curve

def plot_things(samp_MJD, cleaned_light_curve, line_a, line_b, line_c, line_d, line_e, line_f, color_c):
    
    destroy_list_a = [i for i in range(len(cleaned_light_curve)) if cleaned_light_curve[i] <= line_a] #
    destroy_list_b = [i for i in range(len(cleaned_light_curve)) if cleaned_light_curve[i] >= line_b and samp_MJD[i] <= line_c ]
    destroy_list_c = [i for i in range(len(cleaned_light_curve)) if cleaned_light_curve[i] >= line_b and samp_MJD[i] >= line_d and samp_MJD[i] <= line_e]
    destroy_list_d = [i for i in range(len(cleaned_light_curve)) if cleaned_light_curve[i] >= line_b and samp_MJD[i] >= line_f]
   
    f_destroy = destroy_list_a + destroy_list_b + destroy_list_c + destroy_list_d
    
    scatter  =plt.scatter(samp_MJD,  cleaned_light_curve, alpha = 0.55, s=12, color=color_c)
    plt.gca().invert_yaxis()
    plt.grid()
    plt.legend( color_c)
    filter_l_alpha = 0.05
    filter_l_style = 'dashed'
    
    print("N filtered points:", len(f_destroy))
    for o in range(0, len(f_destroy)):
        #print(f_destroy[o])
        plt.axvline(x=samp_MJD[f_destroy[o]], c=color_c, alpha=filter_l_alpha,linestyle=filter_l_style)
        plt.axhline(y=cleaned_light_curve[f_destroy[o]], c=color_c, alpha=filter_l_alpha,linestyle=filter_l_style)
            
    plt.axhline(y=line_a, c=color_c, alpha=0.9,linestyle=filter_l_style) 
    plt.axhline(y=line_b, c=color_c, alpha=0.9,linestyle=filter_l_style) 
            
    plt.axvline(x=line_c, c=color_c, alpha=0.9,linestyle=filter_l_style) 
    plt.axvline(x=line_d, c=color_c, alpha=0.9,linestyle=filter_l_style) 
            
    plt.axvline(x=line_e, c=color_c, alpha=0.9,linestyle=filter_l_style) 
    plt.axvline(x=line_f, c=color_c, alpha=0.9,linestyle=filter_l_style) 
    plt.show()
    
    
def EECep_light_curve_loader_n_cleaner(target_data):
    filename = 'ee_cep_2014.txt'
    print("Reading:", filename )
    file1 = open('./EE_Cep_light_curve/'+ filename, 'r')
    ee_lines = file1.readlines()
    
    B_arr = []
    V_arr = []
    R_arr = []
    I_arr = []

    B_JD = []
    V_JD = []
    R_JD = []
    I_JD = []

    for things in ee_lines:
        if ',B,' in things:
            B_arr.append( float(things.split(",")[1]))
            B_JD.append(  float(things.split(",")[0]))
        elif ',V,' in things:
            V_arr.append( float(things.split(",")[1]))
            V_JD.append(  float(things.split(",")[0]))
        elif ',R,' in things:
            R_arr.append( float(things.split(",")[1]))
            R_JD.append(  float(things.split(",")[0]))
        elif ',I,' in things:
            I_arr.append( float(things.split(",")[1]))
            I_JD.append(  float(things.split(",")[0]))

    t_b = Time(B_JD , format='jd', scale='utc')
    t_v = Time(V_JD , format='jd', scale='utc')
    t_r = Time(R_JD , format='jd', scale='utc')
    t_i = Time(I_JD , format='jd', scale='utc')
        
    color_arr = [B_arr, V_arr, R_arr, I_arr]
    mjd_arr_k = [t_b, t_v, t_r, t_i, target_data[2]] #Time arg
    
    color_arr_cleaned= []
    mjd_arr_k_cleaned=[]
    
    #originals
    c_arr = ['blue',  'green',  'red', 'purple']
    c_labs = ['B', 'V', 'R', 'I', 'pol_data']
    
    for ka in range(0, len(mjd_arr_k)-1):
        scatter  =plt.scatter(mjd_arr_k[ka][:].mjd,  color_arr[ka][:], alpha = 0.55, s=12, color=c_arr[ka], label=c_labs[ka])
    plt.gca().invert_yaxis()
    plt.grid()
    plt.show()
    
    for ka in range(0, len(mjd_arr_k)-1):
        if(ka == 0): #B What was down will be up
            line_a = 11.0844
            line_b = 11.265            
            
            line_c = 56790
            line_d = 56940
            
            line_e = 58880
            line_f = 58994
            
            plot_things(mjd_arr_k[ka][:].mjd,  color_arr[ka][:], line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
            cleaned_mjd, cleaned_curve = destroy_outliers(mjd_arr_k[ka][:].mjd, color_arr[ka][:], line_a, line_b, line_c, line_d, line_e, line_f)
            #plot_things(cleaned_mjd,  cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
            
            cleaned_mjd, cleaned_curve = destroy_outliers(cleaned_mjd, cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f)
            plot_things(cleaned_mjd,  cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
            
            color_arr_cleaned.append(cleaned_curve)
            mjd_arr_k_cleaned.append(cleaned_mjd)
            
        if(ka == 1): #V but default Green            
            line_a = 10.745
            line_b = 10.855       
            
            line_c = 56800
            line_d = 56942
            
            line_e = 58884
            line_f = 58994
            
            plot_things(mjd_arr_k[ka][:].mjd,  color_arr[ka][:], line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
            cleaned_mjd, cleaned_curve = destroy_outliers(mjd_arr_k[ka][:].mjd, color_arr[ka][:], line_a, line_b, line_c, line_d, line_e, line_f)
            #plot_things(cleaned_mjd,  cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
            
            cleaned_mjd, cleaned_curve = destroy_outliers(cleaned_mjd, cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f)
            plot_things(cleaned_mjd,  cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
            
            color_arr_cleaned.append(cleaned_curve)
            mjd_arr_k_cleaned.append(cleaned_mjd)            
            
        elif(ka == 2): #R. The IR excess begins to emerge
            line_a = 10.474
            line_b = 10.620       
            
            line_c = 56800
            line_d = 56950
            
            line_e = 58890
            line_f = 58990
            
            plot_things(mjd_arr_k[ka][:].mjd,  color_arr[ka][:], line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
            cleaned_mjd, cleaned_curve = destroy_outliers(mjd_arr_k[ka][:].mjd, color_arr[ka][:], line_a, line_b, line_c, line_d, line_e, line_f)
            #plot_things(cleaned_mjd,  cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
            cleaned_mjd, cleaned_curve = destroy_outliers(cleaned_mjd, cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f)
            plot_things(cleaned_mjd,  cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
            
            color_arr_cleaned.append(cleaned_curve)
            mjd_arr_k_cleaned.append(cleaned_mjd)
                        
        elif(ka == 3): #I where the excess finna be
            line_a = 10.1092
            line_b = 10.230      
            
            line_c = 56820
            line_d = 56940
            
            line_e = 58885
            line_f = 58998
            
            plot_things(mjd_arr_k[ka][:].mjd,  color_arr[ka][:], line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
            cleaned_mjd, cleaned_curve = destroy_outliers(mjd_arr_k[ka][:].mjd, color_arr[ka][:], line_a, line_b, line_c, line_d, line_e, line_f)
            #plot_things(cleaned_mjd,  cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
            
            cleaned_mjd, cleaned_curve = destroy_outliers(cleaned_mjd, cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f)
            plot_things(cleaned_mjd,  cleaned_curve, line_a, line_b, line_c, line_d, line_e, line_f,c_arr[ka])
                       
            color_arr_cleaned.append(cleaned_curve)
            mjd_arr_k_cleaned.append(cleaned_mjd)         
        
    mjd_arr_k_cleaned.append(target_data[2]) #IT appends this item in the target
    
    for ka in range(0, len(mjd_arr_k_cleaned)-1):
        scatter  =plt.scatter(mjd_arr_k_cleaned[ka][:],  color_arr_cleaned[ka][:], alpha = 0.55, s=12, color=c_arr[ka], label=c_labs[ka])
    plt.gca().invert_yaxis()
    plt.grid()
    plt.show()
    
    return(color_arr_cleaned, mjd_arr_k_cleaned)   
    #return(color_arr[ka][:], mjd_arr_k[ka][:].mjd)
            
def EECep_light_curve_split_A(target_data, MJD_cuts_A, MJD_cuts_B, MJD_cuts_C, point_MJD , lc, verb_t, mark_td, mark_t, true_pa, title_t):
#def EECep_stacked_based(target_data_PD, target_data_PA, MJD_L_cutoff , MJD_cutoff, point_MJD , lc, verb_t, mark_td, mark_t, true_pa, title_t):
    #You can build off of this.
    #Takes in a tuple a tuple a tuple, and a tuple. You need another tuple for the bottom half of the data
    #And thechnically I guess you woul need another for 2020 as well. Like a four stack. It is the only way.
    filename = 'ee_cep_2014.txt'
    print("Reading:", filename )
    file1 = open('./EE_Cep_light_curve/'+ filename, 'r')
    ee_lines = file1.readlines()
    
    if(verb_t):
        print("Plotting:",txt_arg)
        print("Len all:", len(target_data[0]), len(target_data[1]), len(target_data[2]))
        for h in range(0, len(target_data[0])):
            print(h, target_data[0][h] , target_data[1][h] , target_data[2][h].value)
    
    B_arr = []
    V_arr = []
    R_arr = []
    I_arr = []

    B_JD = []
    V_JD = []
    R_JD = []
    I_JD = []

    for things in ee_lines:
        if ',B,' in things:
            B_arr.append( float(things.split(",")[1]))
            B_JD.append(  float(things.split(",")[0]))
        elif ',V,' in things:
            V_arr.append( float(things.split(",")[1]))
            V_JD.append(  float(things.split(",")[0]))
        elif ',R,' in things:
            R_arr.append( float(things.split(",")[1]))
            R_JD.append(  float(things.split(",")[0]))
        elif ',I,' in things:
            I_arr.append( float(things.split(",")[1]))
            I_JD.append(  float(things.split(",")[0]))

    t_b = Time(B_JD , format='jd', scale='utc')
    t_v = Time(V_JD , format='jd', scale='utc')
    t_r = Time(R_JD , format='jd', scale='utc')
    t_i = Time(I_JD , format='jd', scale='utc')
        
    color_arr = [B_arr, V_arr, R_arr, I_arr]
    mjd_arr_k = [t_b, t_v, t_r, t_i, target_data[2]] #Time arg
    
    #originals
    c_arr = ['blue',  'green',  'red', 'purple']
    c_labs = ['B', 'V', 'R', 'I', 'pol_d']
    
    index_cutt_A = []
    index_L_cutt_A = []
    
    for jk in range(0, len(mjd_arr_k)):
        idx = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_A[1]).argmin()
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_A[0]).argmin()
        index_cutt_A.append(int(idx))
        index_L_cutt_A.append(int(idxL))
        
    index_cutt_B = []
    index_L_cutt_B = []
        
    for jk in range(0, len(mjd_arr_k)):
        idx = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_B[1]).argmin()
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_B[0]).argmin()
        index_cutt_B.append(int(idx))
        index_L_cutt_B.append(int(idxL))
        
    index_cutt_C = []
    index_L_cutt_C = []
        
    for jk in range(0, len(mjd_arr_k)):
        idx = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_C[1]).argmin()
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_C[0]).argmin()
        index_cutt_C.append(int(idx))
        index_L_cutt_C.append(int(idxL)) 

    #print("Index cutts A (top):", index_L_cutt_A, index_cutt_A) #For all colours and the target I assume
    #print("Index cutts B (bottom):", index_L_cutt_B, index_cutt_B) #These are array indices
         
    #So this is where you do it.
    #Single Axis is time
    #Just plot me Light curves two slices...
    
    #Plot of Light curve A with index cutt A
    #Plot of Light curve B with index cutt B
    
    #Populate with data
    
    #Profit
    fig = plt.figure(figsize=(24, 18))
    fig, axs = plt.subplots(3, sharex=False, sharey=False)
    
    if(true_pa):
        if(title_t):
            fig.suptitle('Polarization degree (PD), true Position angle (PA) and EE CEP light curve', fontsize=32)    
    else:
        if(title_t):
            fig.suptitle('Polarization degree (PD), Position angle (PA) and EE CEP light curve', fontsize=32)
    
    for ka in range(0, len(mjd_arr_k)-1):
        axs[0].scatter(mjd_arr_k[ka][index_L_cutt_A[ka]:index_cutt_A[ka]].mjd,  color_arr[ka][index_L_cutt_A[ka]:index_cutt_A[ka]], alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])
        
    axs[0].grid()
    axs[0].tick_params(axis="y", labelsize=12)
    axs[0].tick_params(axis="x", labelsize=12)
    axs[0].invert_yaxis()
    axs[0].set_ylabel('Magnitude, (m)', fontsize=12)
            
    for ka in range(0, len(mjd_arr_k)-1):
        axs[1].scatter(mjd_arr_k[ka][index_L_cutt_B[ka]:index_cutt_B[ka]].mjd,  color_arr[ka][index_L_cutt_B[ka]:index_cutt_B[ka]], alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])
                    
    axs[1].grid()
    axs[1].tick_params(axis="y", labelsize=12)
    axs[1].tick_params(axis="x", labelsize=12) 
    axs[1].invert_yaxis() 
    axs[1].set_ylabel('Magnitude, (m)', fontsize=12)
    
    if(MJD_cuts_C[1] >= 59354.9925 ): #one is plus one. The other one is not
        markers, caps, bars = axs[2].errorbar(target_data[2][index_L_cutt_C[4]:index_cutt_A[4]+1].mjd, target_data[0][index_L_cutt_C[4]:index_cutt_C[4]+1], yerr=target_data[1][index_L_cutt_C[4]:index_cutt_C[4]+1], xerr =[0]*len(target_data[0][index_L_cutt_C[4]:index_cutt_C[4]+1]), color='black' ,fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PD')
        
    else:
        markers, caps, bars = axs[2].errorbar(target_data[2][index_L_cutt_C[4]:index_cutt_C[4]].mjd, target_data[0][index_L_cutt_C[4]:index_cutt_C[4]], yerr=target_data[1][index_L_cutt_C[4]:index_cutt_C[4]], xerr =[0]*len(target_data[0][index_L_cutt_C[4]:index_cutt_C[4]]), color='black', fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PD')
        
    [bar.set_alpha(0.2) for bar in bars]
    [cap.set_alpha(0.72) for cap in caps]
        
    axs[2].grid()
    axs[2].tick_params(axis="y", labelsize=12)
    axs[2].tick_params(axis="x", labelsize=12)
    if(target_data[3] == 'PD'):
        axs[2].set_ylabel('PD, (%)', fontsize=12)
    elif(target_data[3] == 'PA'):
        axs[2].set_ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=12)
    
    axs[2].set_xlabel('Time, (MJD)', fontsize=12)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    #for ax in axs:
    #    ax.label_outer()
        
    fig.tight_layout()
    plt.show()
    
def EECep_light_curve_split_B(target_data, MJD_cuts_A, MJD_cuts_B, MJD_cuts_C,cleaned_lc, new_arr ,  point_MJD , lc, verb_t, mark_td, mark_t, true_pa, title_t):         
    color_arr = cleaned_lc #This is cleaned light curve
    mjd_arr_k = new_arr #Time arg. This is new arr
    
    #originals
    c_arr = ['blue',  'green',  'red', 'purple']
    c_labs = ['B', 'V', 'R', 'I', 'pol_d']
    
    index_cutt_A = []
    index_L_cutt_A = []
    
    for jk in range(0, len(mjd_arr_k)):
        idx = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_A[1]).argmin()
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_A[0]).argmin()
        index_cutt_A.append(int(idx))
        index_L_cutt_A.append(int(idxL))
        
    index_cutt_B = []
    index_L_cutt_B = []
        
    for jk in range(0, len(mjd_arr_k)):
        idx = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_B[1]).argmin()
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_B[0]).argmin()
        index_cutt_B.append(int(idx))
        index_L_cutt_B.append(int(idxL))
        
    index_cutt_C = []
    index_L_cutt_C = []
        
    for jk in range(0, len(mjd_arr_k)):
        idx = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_C[1]).argmin()
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts_C[0]).argmin()
        index_cutt_C.append(int(idx))
        index_L_cutt_C.append(int(idxL)) 

    #print("Index cutts A (top):", index_L_cutt_A, index_cutt_A) #For all colours and the target I assume
    #print("Index cutts B (bottom):", index_L_cutt_B, index_cutt_B) #These are array indices
         
    #So this is where you do it.
    #Single Axis is time
    #Just plot me Light curves two slices...
    
    #Plot of Light curve A with index cutt A
    #Plot of Light curve B with index cutt B
    
    #Populate with data
    
    #Profit
    fig = plt.figure(figsize=(24, 18))
    fig, axs = plt.subplots(3, sharex=False, sharey=False)
    
    if(true_pa):
        if(title_t):
            fig.suptitle('Polarization degree (PD), true Position angle (PA) and EE CEP light curve', fontsize=32)    
    else:
        if(title_t):
            fig.suptitle('Polarization degree (PD), Position angle (PA) and EE CEP light curve', fontsize=32)
    
    for ka in range(0, len(mjd_arr_k)-1):
        axs[0].scatter(mjd_arr_k[ka][index_L_cutt_A[ka]:index_cutt_A[ka]].mjd,  color_arr[ka][index_L_cutt_A[ka]:index_cutt_A[ka]], alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])
        
    axs[0].grid()
    axs[0].tick_params(axis="y", labelsize=12)
    axs[0].tick_params(axis="x", labelsize=12)
    axs[0].invert_yaxis()
    axs[0].set_ylabel('Magnitude, (m)', fontsize=12)
            
    for ka in range(0, len(mjd_arr_k)-1):
        axs[1].scatter(mjd_arr_k[ka][index_L_cutt_B[ka]:index_cutt_B[ka]].mjd,  color_arr[ka][index_L_cutt_B[ka]:index_cutt_B[ka]], alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])
                    
    axs[1].grid()
    axs[1].tick_params(axis="y", labelsize=12)
    axs[1].tick_params(axis="x", labelsize=12) 
    axs[1].invert_yaxis() 
    axs[1].set_ylabel('Magnitude, (m)', fontsize=12)
    
    if(MJD_cuts_C[1] >= 59354.9925 ): #one is plus one. The other one is not
        markers, caps, bars = axs[2].errorbar(target_data[2][index_L_cutt_C[4]:index_cutt_A[4]+1].mjd, target_data[0][index_L_cutt_C[4]:index_cutt_C[4]+1], yerr=target_data[1][index_L_cutt_C[4]:index_cutt_C[4]+1], xerr =[0]*len(target_data[0][index_L_cutt_C[4]:index_cutt_C[4]+1]), color='black' ,fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PD')
        
    else:
        markers, caps, bars = axs[2].errorbar(target_data[2][index_L_cutt_C[4]:index_cutt_C[4]].mjd, target_data[0][index_L_cutt_C[4]:index_cutt_C[4]], yerr=target_data[1][index_L_cutt_C[4]:index_cutt_C[4]], xerr =[0]*len(target_data[0][index_L_cutt_C[4]:index_cutt_C[4]]), color='black', fmt=mark_t, markersize=16, capsize=10, capthick=5, label='PD')
        
    [bar.set_alpha(0.2) for bar in bars]
    [cap.set_alpha(0.72) for cap in caps]
        
    axs[2].grid()
    axs[2].tick_params(axis="y", labelsize=12)
    axs[2].tick_params(axis="x", labelsize=12)
    if(target_data[3] == 'PD'):
        axs[2].set_ylabel('PD, (%)', fontsize=12)
    elif(target_data[3] == 'PA'):
        axs[2].set_ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=12)
    
    axs[2].set_xlabel('Time, (MJD)', fontsize=12)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    #for ax in axs:
    #    ax.label_outer()
        
    fig.tight_layout()
    plt.show()