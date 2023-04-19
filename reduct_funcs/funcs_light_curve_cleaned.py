import astropy
import gzip
import glob
import numpy as np
import math
import xlrd
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
from astropy.visualization.mpl_normalize import ImageNormalize
#from astropy.visualization import SqrtStftch

from reduct_funcs import funcs_utils

#Functions that plot pol data with light curve?
#EECep_light_curve_based_pol    10 args
#EECep_stacked_based_pol        11 args
#EECep_stacked_based            10 args

#EECep_light_curve_loader_n_cleaner
#EECep_light_curve_split_A
#EECep_light_curve_split_B

#destroy outliers
#plot things

#Misc Funtions

def load_default_light_curve():
    filename = 'EECEPLightCurve.txt'
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

    color_arr = [B_arr, V_arr[57:], R_arr, I_arr]
    mjd_arr_k = [t_b, t_v[57:], t_r, t_i]
            
    return(color_arr, mjd_arr_k) 

def EECep_light_curve_stacked(PD_data, #PD, PDerror, MJD
                              PA_data,
                              light_curve_data,
                              MJD_cuts,
                              verb_t=False,
                              true_pa=False, 
                              title_t=False):
    """
    A function that takes in polarization master plot.
    
    Parameters
    ----------
    target_data : 
        Polarization Data
    MJD_L_cutoff : 
        Left side MJD cutoff. More backwards in time    
    MJD_cutoff  : Tuple  
        Right side MJD cutoff. More backwards in time 
    mark_t : 
        Marker type
    lc : bool
        Blah
    verb_t : bool
        Blah
    true_pa : bool
        Blah 
    title_t : bool
        Blah
    """
    if len(light_curve_data) == 0:
        print("Loading Default light curve")
        def_light_curve = load_default_light_curve()
        #TODO: Create function and routines that loads light curve data. Give me the same standard format
        color_arr = [def_light_curve[0][0], 
                     def_light_curve[0][1], 
                     def_light_curve[0][2], 
                     def_light_curve[0][3]]

        t_b = Time(def_light_curve[1][0] , format='mjd', scale='utc')
        t_v = Time(def_light_curve[1][1], format='mjd', scale='utc')
        t_r = Time(def_light_curve[1][2] , format='mjd', scale='utc')
        t_i = Time(def_light_curve[1][3] , format='mjd', scale='utc')
        t_data = [Time(x).mjd for x in PD_data[2]]         
    
    else:
        color_arr = [light_curve_data[0][0], 
                     light_curve_data[0][1], 
                     light_curve_data[0][2], 
                     light_curve_data[0][3]]

        t_b = Time(light_curve_data[1][0] , format='mjd', scale='utc')
        t_v = Time(light_curve_data[1][1], format='mjd', scale='utc')
        t_r = Time(light_curve_data[1][2] , format='mjd', scale='utc')
        t_i = Time(light_curve_data[1][3] , format='mjd', scale='utc')
        t_data = [Time(x).mjd for x in PD_data[2]]         

    mjd_arr_k = [t_b, t_v, t_r, t_i, Time(t_data , format='mjd', scale='utc') ]
    
    #originals
    c_arr = ['blue',  'green',  'red', 'purple']
    c_labs = ['B', 'V', 'R', 'I', 'pol_d']
        
    index_L_cutt = []
    index_R_cutt = []
    
    for jk in range(0, len(mjd_arr_k)):        
        #Find the index of the value that most corresponds to the cutoff
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float') - MJD_cuts[0]).argmin()
        idxR = np.abs(mjd_arr_k[jk].to_value('mjd', 'float') - MJD_cuts[1]).argmin() 
        
        index_L_cutt.append(int(idxL))  
        index_R_cutt.append(int(idxR))
        
    fig = plt.figure(figsize=(36, 12))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True, sharey=False)    
    
    markers, caps, bars = axs[0].errorbar(
            mjd_arr_k[4][index_L_cutt[4]:index_R_cutt[4]+1].mjd, 
            PD_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
            yerr=PD_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
            xerr =[0]*len(PD_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
            color='black', fmt='*', markersize=16, capsize=10, capthick=5, label='PD')
    
    [bar.set_alpha(0.15) for bar in bars]
    [cap.set_alpha(0.645) for cap in caps]
    
    markers, caps, bars = axs[1].errorbar(
            mjd_arr_k[4][index_L_cutt[4]:index_R_cutt[4]+1].mjd, 
            PA_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
            yerr=PA_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
            xerr =[0]*len(PA_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
            color='black', fmt='*', markersize=16, capsize=10, capthick=5, label='PA')
    
    axs[0].set_ylabel('PD, (%)', fontsize=32)
    axs[0].tick_params(axis="y", labelsize=28)

    if(true_pa):
        axs[1].set_ylabel('true PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
    else:
        axs[1].set_ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
    axs[1].tick_params(axis="y", labelsize=28)
    
    for ka in range(0, len(mjd_arr_k)-1): #this is a recursive call. Just recursive call the third plot
        axs[2].scatter(mjd_arr_k[ka][index_L_cutt[ka]:index_R_cutt[ka]].mjd, 
                    color_arr[ka][index_L_cutt[ka]:index_R_cutt[ka]], 
                    alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])
       
    if(true_pa):
        fig.suptitle('Polarization degree (PD), true Position angle (PA) and EE CEP light curve', 
                         fontsize=32)    
    else:
        fig.suptitle('Polarization degree (PD), Position angle (PA) and EE CEP light curve', fontsize=32)
        
    axs[0].grid()
    axs[1].grid()
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
            
    
def EECep_light_curve_superposed(pol_data, #PD, PDerror, MJD
                     light_curve_data, #colour = [0], MJD = [1]
                     MJD_cuts, #Tuple... Left and Right
                     txt_arg='',
                     lc=False,
                     verb_t=False,
                     true_pa=False, 
                     title_t=False):
    """
    A function that takes in polarization master plot.
    
    Parameters
    ----------
    target_data : 
        Polarization Data
    MJD_L_cutoff : 
        Left side MJD cutoff. More backwards in time    
    MJD_cutoff  : Tuple  
        Right side MJD cutoff. More backwards in time 
    mark_t : 
        Marker type
    lc : bool
        Blah
    verb_t : bool
        Blah
    true_pa : bool
        Blah 
    title_t : bool
        Blah
    """
    if len(light_curve_data) == 0:
        print("Loading Default light curve")
        def_light_curve = load_default_light_curve()
        #TODO: Create function and routines that loads light curve data. Give me the same standard format
        color_arr = [def_light_curve[0][0], 
                     def_light_curve[0][1], 
                     def_light_curve[0][2], 
                     def_light_curve[0][3]]

        t_b = Time(def_light_curve[1][0] , format='mjd', scale='utc')
        t_v = Time(def_light_curve[1][1], format='mjd', scale='utc')
        t_r = Time(def_light_curve[1][2] , format='mjd', scale='utc')
        t_i = Time(def_light_curve[1][3] , format='mjd', scale='utc')
        t_data = [Time(x).mjd for x in pol_data[2]]         
    
    else:
        color_arr = [light_curve_data[0][0], 
                     light_curve_data[0][1], 
                     light_curve_data[0][2], 
                     light_curve_data[0][3]]

        t_b = Time(light_curve_data[1][0] , format='mjd', scale='utc')
        t_v = Time(light_curve_data[1][1], format='mjd', scale='utc')
        t_r = Time(light_curve_data[1][2] , format='mjd', scale='utc')
        t_i = Time(light_curve_data[1][3] , format='mjd', scale='utc')
        t_data = [Time(x).mjd for x in pol_data[2]]         

    mjd_arr_k = [t_b, t_v, t_r, t_i, Time(t_data , format='mjd', scale='utc') ]
    
    c_arr = ['blue',  'green',  'red', 'purple']
    c_labs = ['B', 'V', 'R', 'I', 'pol_d']
        
    index_L_cutt = []
    index_R_cutt = []
    
    for jk in range(0, len(mjd_arr_k)):
        #For each data string
        
        #Find the index of the value that most corresponds to the cutoff
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float') - MJD_cuts[0]).argmin()
        idxR = np.abs(mjd_arr_k[jk].to_value('mjd', 'float') - MJD_cuts[1]).argmin() 
        
        index_L_cutt.append(int(idxL))  
        index_R_cutt.append(int(idxR))
    
    #fig = plt.figure(figsize=(36, 12))    
    fig, ax1 = plt.subplots(figsize=(36, 12))
    
    for ka in range(0, len(mjd_arr_k)-1): #this is a recursive call. Just recursive call the third plot
        ax1.scatter(mjd_arr_k[ka][index_L_cutt[ka]:index_R_cutt[ka]].mjd, 
                    color_arr[ka][index_L_cutt[ka]:index_R_cutt[ka]], 
                    alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])

    fig.gca().invert_yaxis()
    #"""
    ax2 = ax1.twinx()
    markers, caps, bars = ax2.errorbar(
            mjd_arr_k[4][index_L_cutt[4]:index_R_cutt[4]+1].mjd, 
            pol_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
            yerr=pol_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
            xerr =[0]*len(pol_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
            color='black', fmt='*', markersize=16, capsize=10, capthick=5, label='PA')
    
    #plt.set_ylabel('PD, (%)', fontsize=32)
    #plt.tick_params(axis="y", labelsize=28)
    #"""
    [bar.set_alpha(0.2) for bar in bars]
    [cap.set_alpha(0.72) for cap in caps]  

    
    """
    if(txt_arg == 'PD'):
        ax2.set_title("EE Cep light curve and Polarization Degree ("+txt_arg+")", fontsize=32)
        ax2.set_ylabel('PD, (%)', fontsize=32)
        ax2.tick_params(axis="y", labelsize=28)
        
    elif(txt_arg == 'PA'):
        if(true_pa):
            ax1.set_title("EE Cep light curve and true Position Angle ("+txt_arg+")", fontsize=32)

        else
            ax1.set_title("EE Cep light curve and Position Angle ("+txt_arg+")", fontsize=32)
    else:
        print("Incorrect Arguement Supplied")
    """
    ax2.set_ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
    ax2.tick_params(axis="y", labelsize=28)

    
    ax1.set_ylabel('magnitude, (m)', fontsize=32)
    ax1.set_xlabel('Time, (MJD)', fontsize=32)
    ax1.tick_params(axis="y", labelsize=28)
    ax1.tick_params(axis="x", labelsize=28)
    #plt.xticks(fontsize=28)
    #plt.yticks(fontsize=28)
    
    plt.grid()
    
    
        
    """
    fig = plt.figure(figsize=(36, 12))
    #gs = fig.add_gridspec(3, hspace=0)
    #axs = gs.subplots(sharex=True, sharey=False)    
    
    markers, caps, bars = axs[0].errorbar(
            mjd_arr_k[4][index_L_cutt[4]:index_R_cutt[4]+1].mjd, 
            PD_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
            yerr=PD_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
            xerr =[0]*len(PD_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
            color='black', fmt='*', markersize=16, capsize=10, capthick=5, label='PD')
    
    [bar.set_alpha(0.15) for bar in bars]
    [cap.set_alpha(0.645) for cap in caps]
    
    markers, caps, bars = axs[1].errorbar(
            mjd_arr_k[4][index_L_cutt[4]:index_R_cutt[4]+1].mjd, 
            PA_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
            yerr=PA_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
            xerr =[0]*len(PA_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
            color='black', fmt='*', markersize=16, capsize=10, capthick=5, label='PA')
    
    axs[0].set_ylabel('PD, (%)', fontsize=32)
    axs[0].tick_params(axis="y", labelsize=28)

    if(true_pa):
        axs[1].set_ylabel('true PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
    else:
        axs[1].set_ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
    axs[1].tick_params(axis="y", labelsize=28)
    
    for ka in range(0, len(mjd_arr_k)-1): #this is a recursive call. Just recursive call the third plot
        axs[2].scatter(mjd_arr_k[ka][index_L_cutt[ka]:index_R_cutt[ka]].mjd, 
                    color_arr[ka][index_L_cutt[ka]:index_R_cutt[ka]], 
                    alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])
       
    if(true_pa):
        fig.suptitle('Polarization degree (PD), true Position angle (PA) and EE CEP light curve', 
                         fontsize=32)    
    else:
        fig.suptitle('Polarization degree (PD), Position angle (PA) and EE CEP light curve', fontsize=32)
        
    axs[0].grid()
    axs[1].grid()
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
    """
    """    
    plt.show()
    """ 
    