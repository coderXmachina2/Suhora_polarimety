import astropy
import gzip
import glob
import math
import pickle
import re
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astropy.time import Time, TimeJD
from astropy.visualization import simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from copy import deepcopy as cp
from datetime import datetime, timedelta

from mpl_toolkits.axes_grid1 import make_axes_locatable
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.aperture import EllipticalAperture, EllipticalAnnulus
from photutils.aperture import aperture_photometry
from reduct_funcs import funcs_utils
from scipy import stats
from scipy import interpolate
from scipy.interpolate import interp1d

   
def check_pol_std():
    """
    A function that prints out high polarization and zero polarization standards. Takes in no input arguments. Please double check this list because some of them are targets. Does not return anything.

    Parameters
    ----------
    None.
    """
    high_pol_std = ['BDp64106', 
                'BDp59389', 
                'HD19820', 
                'HD25443',
                'HD215806',
                'BDp25727', 
                'HD251204', #spotted
                'HD160529', 
                'HD161056',
                'Hiltner960',
                'VICyg12',
                'HD204827' #spotted
               ]
    
    low_pol_std = ['HD12021', 
                   'HD14069', 
                   'HD21447', 
                   'G191B2B', #spotted
                   'HD64299', 
                   'HD94851', 
                   'GD319', 
                   'BDp33_2642', 
                   'HD154892', 
                   'BDp32_3739',
                   'BDp28_4211',
                   'HD212311'] #spotted
    
    print("High Pol standard:", high_pol_std, "\n")
    print("Low Pol standard:", low_pol_std, "\n")
    
def calc_pd2(input_data, 
             plot_title, 
             plot_c, 
             perc_arg= False, 
             verbose_calc_pd=False, 
             verbose_data=False,
             sv_pold_img=False):
    """
    A function that takes in data and calculates polarization degree (PD). Does not return anything. Can be used in preliminary analysis of polarization data.

    Parameters
    ----------
    input_data : tuple
        Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
    plot_title : str,
        String that defines the title of the plot. Used in the filename if the pol data is saved to file.
    plot_c : str,
        String that defines the color of the plot.
    perc_arg : bool, optional
        Multiply calculated pd values by 100 to get percentage polarization. 
    verbose_calc_pd : bool, optional
        Prints calculations for extra verbosity.  False by default
    verbose_data : bool, optional
        Prints calculations for extra verbosity.  False by default 
    sv_arg : bool, optional
        Saves image to file.  False by default 
    """
    dates = sorted(input_data[1])
    
    print("Calculate and plot polarization degree (PD) for duration", dates[0], "to", dates[-1],"\nwithout returning data")
    pol_d_array = []
    pol_d_err_array = []
    
    for dats in input_data[0]:        
        result = re.search('(.*)_', list(dats.keys())[0][:12])
        targ_name_str = re.search('_(.*)', list(dats.keys())[0])
                
        ###Math
        mean_q_squared = np.mean(dats[list(dats.keys())[0]][0][1:])**2 #mean q squared
        mean_u_squared = np.mean(dats[list(dats.keys())[0]][2][1:])**2 #mean u squared
        
        mean_q_err_squared = np.std(dats[list(dats.keys())[0]][1][1:])**2 #std q squared
        mean_u_err_squared = np.std(dats[list(dats.keys())[0]][3][1:])**2 #std u squared
        
        sum_o_squares = mean_q_squared + mean_u_squared #mean q squared + mean u squared
        
        pol_d = math.sqrt(sum_o_squares)

        pol_d_err  = math.sqrt(((mean_q_squared)/(sum_o_squares))*(mean_q_err_squared)+ ((mean_u_squared)/(sum_o_squares))*(mean_u_err_squared))
        ###
        
        if(perc_arg):
            if(verbose_calc_pd):
                print(targ_name_str.group(1), "MJD:", result.group(1), pol_d*100, u"\u00B1",  pol_d_err*100)
            pol_d_array.append(pol_d*100)
            pol_d_err_array.append(pol_d_err*100)
        else:
            if(verbose_calc_pd):
                print("MJD:", result.group(1), pol_d, u"\u00B1",  pol_d_err)
            pol_d_array.append(pol_d)
            pol_d_err_array.append(pol_d_err)

    t = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in input_data[1]], format='isot', scale='utc')
 
    plt.errorbar(t.mjd , 
                 pol_d_array, 
                 xerr=[0]*len(pol_d_array), 
                 yerr=pol_d_err_array, 
                 lw=0.75, fmt="^", capsize=10,color=plot_c, alpha=0.9)
        
    plt.title(plot_title)
    
    if(verbose_data):
        for i in range(0, len(t)):
            plt.text(t[i].mjd, 
                     pol_d_array[i], 
                     str(np.round(t[i].mjd))+"\nPD:"+str(np.round(pol_d_array[i],2)),
                     fontsize=12, rotation=45)
            
    if(perc_arg):
        plt.ylabel("PD %")
    else:
        plt.ylabel("PD")        
    plt.xlabel("Data list index")
    plt.grid()
    if(sv_pold_img):
        plt.savefig(plot_title+'.png', bbox_inches='tight',pad_inches=0.1)
    plt.show()
    
def calc_pa2(input_data, 
             plot_title, 
             plot_c,
             deg_arg=False, 
             verbose_calc_pa=False, 
             verbose_data=False,
             sv_polpa_img=False):
    """
    A function that takes in data and calculates position angle (PA). Does not return anything. Can be used in preliminary analysis of polarization data.

    Parameters
    ----------
    input_data : tuple
        Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
    plot_title : str,
        String that defines the title of the plot. Used in the filename if the pol data is saved to file.
    plot_c : str,
        String that defines the color of the plot.
    deg_arg : bool, optional
        Multiply calculated pd values by 180/pi to get pol angle degree. 
    verbose_calc_pa : bool, optional
        Prints calculations for extra verbosity.  False by default
    verbose_data : bool, optional
        Prints calculations for extra verbosity.  False by default 
    sv_polpa_img : bool, optional
        Saves image to file.  False by default 
    """
    dates = sorted(input_data[1])
    
    print("Calculate and plot position angle (PA) for duration", dates[0], "to", dates[-1],"\nwithout returning data")
    pol_pa_array = []
    pol_pa_err_array = []

    for dats in input_data[0]:        
        result = re.search('(.*)_', list(dats.keys())[0][:12])
        targ_name_str = re.search('_(.*)', list(dats.keys())[0] )
        
        ###Please check
        mean_q_err_squared = np.std(dats[list(dats.keys())[0]][1][1:])**2
        mean_u_err_squared = np.std(dats[list(dats.keys())[0]][3][1:])**2
        
        mean_q_squared = np.mean(dats[list(dats.keys())[0]][0][1:])**2
        mean_u_squared = np.mean(dats[list(dats.keys())[0]][2][1:])**2
        
        mean_q = np.mean(dats[list(dats.keys())[0]][0][1:])
        mean_u = np.mean(dats[list(dats.keys())[0]][2][1:])
        
        sum_o_squares = mean_q_squared + mean_u_squared
        
        pol_pa = 0.5*math.atan2(mean_u , mean_q)
        pol_pa_err = math.sqrt(((1/(2*mean_q*(1 + (mean_u_squared/mean_q_squared))))**2 )*(mean_u_err_squared) + ((-1*((mean_u)/(2*sum_o_squares)))**2 )*(mean_q_err_squared))
        #Please check
        
        if(deg_arg):
            if(verbose_calc_pa): #
                print("MJD:", result.group(1), pol_pa*(180/3.142), u"\u00B1",  pol_pa_err*(180/3.142))
            pol_pa_array.append(pol_pa*(180/3.142)) 
            pol_pa_err_array.append(pol_pa_err*(180/3.142)) 
        else:
            if(verbose_calc_pa):
                print("MJD:", result.group(1), pol_pa, u"\u00B1",  pol_pa_err)
            pol_pa_array.append(pol_pa)
            pol_pa_err_array.append(pol_pa_err)

    t = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in input_data[1]], format='isot', scale='utc')

    plt.errorbar(t.mjd, 
                 pol_pa_array, 
                 xerr=[0]*len(pol_pa_array), 
                 yerr=pol_pa_err_array, 
                 lw=0.75, fmt="^", capsize=10, color=plot_c, alpha=0.9) 
    
    plt.title(plot_title)
    plt.grid()
    
    if(verbose_data):
        for i in range(0, len(t)):
            plt.text(t[i].mjd, 
                     pol_pa_array[i], 
                     str(np.round(t[i].mjd))+"\nPA:"+str(np.round(pol_pa_array[i])),
                     fontsize=14, rotation=45)
            
    plt.ylabel("PA")
    plt.xlabel("MJD")
    if(sv_polpa_img):
        plt.savefig(sv_im_str, bbox_inches='tight',pad_inches=0.1)
    plt.show()  

def calc_PD_stability(input_data,              
                      targ_corr_MJD='',
                      verbose_calc_pd=False, 
                      verbose_mjd_align_check=False, 
                      perc_arg=False,
                      to_excel=False, 
                      corr_MJD=False):
    
    """
    A function that takes in data and calculates polarization degree (PD). Returns data to be used in other plot tools

    Parameters
    ----------
    input_data : tuple
        Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
    targ_corr_MJD : str
        String to identify data file for correcting MJD
    verbose_calc_pd : bool, optional
        Multiply calculated pd values by 100 to get percentage polarization degree. 
    verbose_mjd_align_check : bool, optional
        Prints calculations for extra verbosity.  False by default 
    perc_arg : bool, optional
        Saves image to file.  False by default
    to_excel : bool, optional
        Saves image to file.  False by default 
    corr_MJD : bool, optional
        Saves image to file.  False by default 
    """
    dates = sorted(input_data[1])
    objn = list(input_data[0][0].keys())[0][11:]
    
    print("Calculate and return polarization degree (PD) for duration", dates[0], "to", dates[-1], "\nfor", objn)
    
    if(corr_MJD):
        if(targ_corr_MJD=='EECEP'):
            with open('./data_pickles/MJD_corr_pix.pickle', 'rb') as fid: #This is saved in pickle
                corr_da = pickle.load(fid)
        elif(targ_corr_MJD=='bd64'):
            with open('./data_pickles/MJD_corr_64106.pickle', 'rb') as fid: #This is saved in pickle
                corr_da = pickle.load(fid)
        elif(targ_corr_MJD=='g191'):
            with open('./data_pickles/MJD_corr_G191.pickle', 'rb') as fid: #This is saved in pickle
                corr_da = pickle.load(fid)            
        elif(targ_corr_MJD=='215806'):
            with open('./data_pickles/MJD_corr_hd215806.pickle', 'rb') as fid: #This is saved in pickle
                corr_da = pickle.load(fid)
    
    means_arr = []
    means_err_arr = []
    
    mjd_strs = []
    
    for dats in input_data[0]:        
        result = re.search('(.*)_', list(dats.keys())[0][:12])
        targ_name_str = re.search('_(.*)', list(dats.keys())[0])

        if(corr_MJD):
            for keyes in corr_da.keys():
                if(result.group(1) == keyes):
                    mjd_strs.append(corr_da[keyes][0]) #Its straight up in mjd format
        
        m_q_sq = np.mean(dats[list(dats.keys())[0]][0][1:])**2 #mean q squared
        m_u_sq = np.mean(dats[list(dats.keys())[0]][2][1:])**2 #mean u squared 

        std_q_err_sq = np.std(dats[list(dats.keys())[0]][0][1:])**2 #std q squared
        std_u_err_sq = np.std(dats[list(dats.keys())[0]][2][1:])**2 #std u squared
        
        sum_sq = m_q_sq + m_u_sq #mean q squared plus mean u squared
                
        pol_d = math.sqrt(sum_sq)        
        pol_d_err = math.sqrt((std_q_err_sq)*((m_q_sq)/(sum_sq)) + ((std_u_err_sq)*((m_u_sq)/(sum_sq))))
        
        if(verbose_mjd_align_check):
            print(targ_name_str.group(1), Time([result.group(1)+'T00:00:00.000000000'])[0].mjd  )
            
        if(perc_arg):
            if(verbose_calc_pd):
                print("Pol D:", pol_d*100,u"\u00B1",pol_d_err*100 )
            means_arr.append(pol_d*100)
            means_err_arr.append(pol_d_err*100)
        else:
            if(verbose_calc_pd):
                print("Pol D:", pol_d,u"\u00B1",pol_d_err)
            means_arr.append(pol_d)
            means_err_arr.append(pol_d_err)
            
    if(corr_MJD):
        t = Time(mjd_strs, scale='utc',format='mjd')
    else:
        t = input_data[1]

    if(verbose_mjd_align_check):
        for a in range(0, len(input_data[0])):
            print(list(input_data[0][a].keys())[0], t[a], t.mjd[a])

    if(to_excel):
        funcs_utils.data_to_excel((t.mjd, means_arr, means_err_arr),'cal_eecep', 'pol_deg')

    return(means_arr,means_err_arr, t)

def calc_PA_stability(input_data,
                      targ_corr_MJD='', 
                      verbose_calc_pa=False,
                      verbose_mjd_align_check=False, 
                      deg_arg=False,
                      to_excel=False, 
                      corr_MJD=False,
                      PA_shift=False):
    """
    A function that takes in data and calculates position angle (PA). Returns data to be used in other plot tools

    Parameters
    ----------
    input_data : tuple
        Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
    targ_corr_MJD : str
        String to identify data file for correcting MJD
    verbose_calc_pa : bool, optional
        Multiply calculated pd values by 180/pi to get position angle in degree degree. 
    verbose_mjd_align_check : bool, optional
        Prints calculations for extra verbosity.  False by default 
    perc_arg : bool, optional
        Saves image to file.  False by default
    to_excel : bool, optional
        Saves image to file.  False by default 
    corr_MJD : bool, optional
        Saves image to file.  False by default 
    """
    dates = sorted(input_data[1])
    objn = list(input_data[0][0].keys())[0][11:]
    
    print("Calculate and return position angle (PA) for duration", dates[0], "to", dates[-1], "\nfor", objn)
    
    if( corr_MJD ):
        if(targ_corr_MJD=='EECEP'):
            with open('./data_pickles/MJD_corr_pix.pickle', 'rb') as fid: #This is saved in pickle
                corr_da = pickle.load(fid)
        elif(targ_corr_MJD=='bd64'):
            with open('./data_pickles/MJD_corr_64106.pickle', 'rb') as fid: #This is saved in pickle
                corr_da = pickle.load(fid)
        elif(targ_corr_MJD=='g191'):
            with open('./data_pickles/MJD_corr_G191.pickle', 'rb') as fid: #This is saved in pickle
                corr_da = pickle.load(fid)            
        elif(targ_corr_MJD=='215806'):
            with open('./data_pickles/MJD_corr_hd215806.pickle', 'rb') as fid: #This is saved in pickle
                corr_da = pickle.load(fid)
            
    means_arr = []
    means_err_arr = []
    
    mjd_strs = []
    
    for dats in input_data[0]:        
        result = re.search('(.*)_', list(dats.keys())[0][:12])
        targ_name_str = re.search('_(.*)', list(dats.keys())[0])
        
        if(corr_MJD):
            for keyes in corr_da.keys():
                if(result.group(1) == keyes):
                    mjd_strs.append(corr_da[keyes][0]) #Its straight up in mjd format
        
        mean_q_err_squared = np.std(dats[list(dats.keys())[0]][0][1:])**2 #taking 
        mean_u_err_squared = np.std(dats[list(dats.keys())[0]][2][1:])**2
        
        mean_q_squared = np.mean(dats[list(dats.keys())[0]][0][1:])**2
        mean_u_squared = np.mean(dats[list(dats.keys())[0]][2][1:])**2
        
        mean_q = np.mean(dats[list(dats.keys())[0]][0][1:])
        mean_u = np.mean(dats[list(dats.keys())[0]][2][1:])
        
        sum_o_squares = mean_q_squared + mean_u_squared
        
        pol_pa = 0.5*math.atan2(mean_u , mean_q )
        
        pol_pa_err = math.sqrt(((1/(2*mean_q*(1 + (mean_u_squared/mean_q_squared))))**2 )*(mean_u_err_squared) + ((-1*( (mean_u)/(2*sum_o_squares)))**2 )*(mean_q_err_squared))
        
        if(verbose_mjd_align_check):
            print(targ_name_str.group(1), Time([result.group(1)+'T00:00:00.000000000'])[0].mjd  )
            
        if(deg_arg):
            if(verbose_calc_pa):
                print("Pol PA:", pol_pa*(180/3.142),u"\u00B1",pol_pa*(180/3.142) )
            means_arr.append(  pol_pa*(180/3.142)) #mean of qs. Calculate PA. PA SJift.
            means_err_arr.append(pol_pa_err*(180/3.142)) #std of list. 
        else:
            if(verbose_calc_pa):
                print("Pol PA:", pol_pa,u"\u00B1",pol_pa_err)
            means_arr.append(pol_pa)
            means_err_arr.append(pol_pa_err)
                
    if(corr_MJD):
        t = Time(mjd_strs, scale='utc',format='mjd') #watch out
    else:
        t = input_data[1]       
    
    if(PA_shift):
        print("Shifting PA...")
        means_arr= np.array(means_arr) + 40
    
    if(verbose_mjd_align_check):
        for a in range(0, len(input_data[0])):
            print(list(input_data[0][a].keys())[0], t[a], t.mjd[a])
            
    if(to_excel):
        funcs_utils.data_to_excel((t.mjd, means_arr, means_err_arr), 'cal_eecep', 'pol_PA')
            
    return(means_arr, means_err_arr, t) #This is where it happens.

def plot_pol_stab(MJD_track,
                  obj_pol, 
                  obj_pol_err, 
                  plot_data,  
                  toggle=False):
    """
    A function that takes in data and calculates position angle (PA). Returns data to be used in other plot tools. More presentable publication quality style.

    Parameters
    ----------
    input_data : tuple
        Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
    targ_corr_MJD : str
        String to identify data file for correcting MJD
    verbose_calc_pa : bool, optional
        Multiply calculated pd values by 180/pi to get position angle in degree degree. 
    verbose_mjd_align_check : bool, optional
        Prints calculations for extra verbosity.  False by default 
    perc_arg : bool, optional
        Saves image to file.  False by default
    to_excel : bool, optional
        Saves image to file.  False by default 
    corr_MJD : bool, optional
        Saves image to file.  False by default 
    """
    
    fig, ax = plt.subplots(figsize=(36, 12))
    
    MJD_track = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in MJD_track], format='isot', scale='utc')
    
    if(plot_data=='PD'):
        markers, caps, bars = ax.errorbar(MJD_track.mjd, obj_pol, yerr=obj_pol_err, xerr =[0]*len(obj_pol),
                        fmt='o',markersize=16, ecolor='blue',capsize=10, capthick=5,  label='PD')
    elif(plot_data=='PA'):
        markers, caps, bars = ax.errorbar(MJD_track.mjd, obj_pol, yerr=obj_pol_err, xerr =[0]*len(obj_pol),
                        fmt='o',markersize=16, ecolor='blue',capsize=10, capthick=5,  label='PA')

    plt.grid()
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    
    if(plot_data=='PD'):
        plt.title('Polarization Degree (PD) versus Time (MJD)', fontsize=32)
        plt.ylabel('PD, (%)', fontsize=28)
        plt.xlabel('Time, (MJD)', fontsize=28)
        if(toggle):
            bounce = 1
            mjds_arr = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in MJD_track], scale='utc',format='isot')
            for k in range(0, len(MJD_track.value)):
                if (bounce == 0):
                    plt.text(float(mjds_arr[k].mjd), 2.8+0.025, k) #The 
                    bounce = 1
                elif (bounce == 1):
                    plt.text(float(mjds_arr[k].mjd), 2.8-0.025, k)
                    bounce = 2
                elif (bounce == 2):
                    plt.text(float(mjds_arr[k].mjd), 2.8+0.06, k)
                    bounce = 3
                elif (bounce == 3):
                    plt.text(float(mjds_arr[k].mjd), 2.8-0.06, k)
                    bounce = 0
                        
    elif(plot_data=='PA'):
        plt.title('Positiona Angle (PA) versus time (MJD)', fontsize=32)
        plt.ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=28)
        plt.xlabel('Time, (MJD)', fontsize=28)
        if(toggle):
            bounce = 1
            mjds_arr = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in MJD_track], scale='utc',format='isot')
            for k in range(0, len(MJD_track.value)):
                if (bounce == 0):
                    plt.text(mjds_arr[k].mjd, -10+0.525, k)
                    bounce = 1
                elif (bounce == 1):
                    plt.text(mjds_arr[k].mjd, -10-0.525, k)
                    bounce = 2
                elif (bounce == 2):
                    plt.text(mjds_arr[k].mjd, -10+1.06, k)
                    bounce = 3
                elif (bounce == 3):
                    plt.text(mjds_arr[k].mjd, -10-1.06, k)
                    bounce = 0
    else:
        print("error. Please supply correct pol type.")
        
    [bar.set_alpha(0.2) for bar in bars]
    [cap.set_alpha(0.6) for cap in caps]
    
    fig.legend(loc="upper right", fontsize=64, borderaxespad=0.86)
    fig.tight_layout()
    plt.show()
    
def plot_pol_stab_double(input_data, 
                         fmt_icon, 
                         color_icon):
    
    """
    A function that takes in data and calculates position angle (PA). Returns data to be used in other plot tools. More presentable publication quality style.

    Parameters
    ----------
    input_data : list
        Lis of data where expectation is:
        input_data[0] = MJD_track_PD
        input_data[1] = MJD_track_PA
        input_data[2] = obj_pol_PD
        input_data[3] = obj_pol_PA
        input_data[4] = obj_pol_PD_err
    icon : list
        List of strings where expectation is:
        fmt_icon[0] = '*'
        fmt_icon[1] = '-'
        Or whatever format icon is preferred
    color_icon : list
        List of strings where expectation is:
        color_icon[0] = 'red'
        color_icon[1] = 'blue'
        Or whatever format color is preferred
    """
    fig, ax1 = plt.subplots(figsize=(36, 12))
    mjds_PD_arr = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in input_data[0]], 
                       scale='utc',
                       format='isot')
    
    markers, caps, bars = ax1.errorbar(mjds_PD_arr.mjd, 
                                       input_data[2], 
                                       yerr=input_data[4], 
                                       xerr =[0]*len(input_data[4]), 
                                       fmt=fmt_icon[0], 
                                       markersize=32, 
                                       capsize=10, 
                                       capthick=5, 
                                       color = color_icon[0], 
                                       label='PD')

    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.xlabel('Time (MJD)', fontsize=28)
    plt.ylabel('PD, (%)', fontsize=24)
    plt.title('Polarization Degree (PD) and Polarization Angle (PA) versus Time (MJD)', fontsize=28)
        
    ax1.grid()
    ax1.plot()
    
    [bar.set_alpha(0.2) for bar in bars]
    [cap.set_alpha(0.6) for cap in caps]
    
    ax2 = ax1.twinx()
    
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.xlabel('Time, (MJD)', fontsize=28)
    plt.ylabel('PA, ('+u'\N{DEGREE SIGN}'+')' , fontsize=24)

    mjds_PA_arr = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in input_data[1]],
                       scale='utc',
                       format='isot')
    
    markers, caps, bars = ax2.errorbar(mjds_PA_arr.mjd, 
                                       input_data[3], 
                                       yerr=input_data[5], 
                                       xerr =[0]*len(input_data[4]),
                                       fmt=fmt_icon[1], 
                                       markersize=24, 
                                       capsize=10, 
                                       capthick=5, 
                                       color = color_icon[1],
                                       label='PA')
    
    [bar.set_alpha(0.2) for bar in bars]
    [cap.set_alpha(0.6) for cap in caps]
    
    ax2.plot()

    fig.legend(loc="upper right", fontsize=36, borderaxespad=3.0)

    fig.tight_layout()
    plt.show()

def plot_q_u_stability(input_data, 
                       q_u_check, 
                       sv_im='', 
                       plot_verbose=False, 
                       verbose=False, 
                       m_plot=False):
    """
    A function that takes in data and calculates position angle (PA). Returns data to be used in other plot tools. More presentable publication quality style.

    Parameters
    ----------
    input_data : tuple
        Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
    q_u_check : str
        String to identify data file for correcting MJD
    sv_im : bool, optional
        Multiply calculated pd values by 180/pi to get position angle in degree degree. 
    plot_verbose : bool, optional
        Prints calculations for extra verbosity.  False by default 
    m_plot : bool, optional
        Saves image to file.  False by default 
    """
    
    means_arr = []
    means_err_arr = []
    
    for dats in input_data[0]:        
        result = re.search('(.*)_', list(dats.keys())[0][:12])
        targ_name_str = re.search('_(.*)', list(dats.keys())[0])
                
        if(q_u_check=='q'):
            if(verbose):
                print(np.mean(dats[list(dats.keys())[0]][0][1:]),u"\u00B1",np.std(dats[list(dats.keys())[0]][0][1:]))
            means_arr.append(np.mean(dats[list(dats.keys())[0]][0][1:]))
            means_err_arr.append(np.std(dats[list(dats.keys())[0]][0][1:]))
        elif(q_u_check=='u'):
            if(verbose):
                print(np.mean(dats[list(dats.keys())[0]][2][1:]),u"\u00B1",np.std(dats[list(dats.keys())[0]][2][1:]))
            means_arr.append(np.mean(dats[list(dats.keys())[0]][2][1:]))
            means_err_arr.append(np.std(dats[list(dats.keys())[0]][2][1:]))
            
    t = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in input_data[1]], format='isot', scale='utc')
       
    if(m_plot):
        print("Plot", q_u_check,"stability")
        fig, ax = plt.subplots()
        markers, caps, bars = ax.errorbar(t.mjd, means_arr, yerr=means_err_arr, xerr =[0]*len(means_arr),
                fmt='o', ecolor='blue',capsize=2, capthick=2)

        plt.title(q_u_check + " stability over time", fontsize=24)
        if(plot_verbose):
            for l in range(0, len(t.mjd)):
                plt.text(t.mjd[l], means_arr[l], 
                         str(round(means_arr[l],6))+u"\u00B1"+str(np.round(means_err_arr[l],4)), 
                         fontsize=12, rotation=45)

        plt.yticks(fontsize = 22)
        plt.xticks(fontsize = 22)
        if(q_u_check=='q'):
            plt.ylabel('q', fontsize=24)
            plt.axhline(np.mean(means_arr), linestyle='dashed', alpha = 0.45 )
            locs, labels = plt.yticks() 

            plt.text(t.mjd[0],
                     np.mean(means_arr) - (locs[1]-locs[0])/6 ,
                     'Mean Stokes Q:'+ str(np.round(np.mean(means_arr),4))+u"\u00B1"+str(np.round(np.mean(means_err_arr), 4  ) ),
                     fontsize=24)
        elif(q_u_check=='u'):
            plt.ylabel('u', fontsize=24)
            plt.axhline(np.mean(means_arr), linestyle='dashed', alpha = 0.45)
            locs, labels = plt.yticks() 

            plt.text(t.mjd[0],
                     np.mean(means_arr) - (locs[1]-locs[0])/6,
                     'Mean Stokes U:'+ str(np.round(np.mean(means_arr),4))+u"\u00B1"+str(np.round(np.mean(means_err_arr), 4)), 
                     fontsize=24)

        plt.axhline(y=0, color = 'black')
        plt.xlabel('MJD',fontsize=24)
        plt.grid()

        [bar.set_alpha(0.22) for bar in bars]
        [cap.set_alpha(0.5) for cap in caps]

        if(sv_im != ''):
            plt.savefig(sv_im,bbox_inches='tight',pad_inches=0.1 )
        plt.show()   
        
    return([means_arr, means_err_arr]) #returns the plot here
    
    """
    print("Plot versus time!")
    plt.errorbar(t.mjd, means_arr, yerr=means_err_arr, xerr =[0,0,0,0,0])
    plt.title(q_u_check + " stability")
    if(plot_verbose):
        for l in range(0, len(t.mjd)):
            plt.text(t.mjd[l], means_arr[l], str(round(means_arr[l], 6))+u"\u00B1"+str(np.round(means_err_arr[l],4)), fontsize=24)
    if(q_u_check=='q'):
        plt.ylabel('q', fontsize=24)
    elif(q_u_check=='u'):
        plt.ylabel('u', fontsize=24)
    plt.xlabel('MJD')
    plt.grid()
    plt.show()
    """

'''
input_data[0] = pol_data
input_data[0] = tstamps
'''
    
def q_n_u_single_plot_v1(input_data,
                         plot_c,   
                         sv_im='', 
                         verbose_MJD_arg=False,
                         pol_deg=False, 
                         only_means=False, 
                         to_excel=False, 
                         retrun_plotelems=False,
                         key_verb=False):
    """
    A function that takes in data and calculates position angle (PA). Returns data to be used in other plot tools. More presentable publication quality style.

    Parameters
    ----------
    input_data : tuple
        Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
    q_u_check : str
        String to identify data file for correcting MJD
    sv_im : bool, optional
        Multiply calculated pd values by 180/pi to get position angle in degree degree. 
    plot_verbose : bool, optional
        Prints calculations for extra verbosity.  False by default 
    m_plot : bool, optional
        Saves image to file.  False by default 
    """    
    if(key_verb):
        print("N data:", len(input_data[0]))
    
    targ_qmeans = []
    targ_umeans = []
    
    targ_qmeans_err = []
    targ_umeans_err = []
    
    targ_qs = []
    targ_us = []
    targ_qstds = []
    targ_ustds = []
        
    for things in input_data[0]:
        for k in range(0, len(things[list(things.keys())[0]][0])):
            targ_qs.append(things[list(things.keys())[0]][0][k])
            targ_us.append(things[list(things.keys())[0]][2][k])
            targ_qstds.append(things[list(things.keys())[0]][1][k])
            targ_ustds.append(things[list(things.keys())[0]][3][k])

        targ_qmeans.append(np.mean(things[list(things.keys())[0]][0][:])) #qmeans
        targ_umeans.append(np.mean(things[list(things.keys())[0]][2][:])) #umeans
        
        targ_qmeans_err.append(np.std(things[list(things.keys())[0]][0][:])) #qstd
        targ_umeans_err.append(np.std(things[list(things.keys())[0]][2][:])) #ustd

    t = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in input_data[1]], format='isot', scale='utc')

    if(only_means):
        plt.scatter(targ_qmeans, targ_umeans,
                    alpha=0.9, color = plot_c,)
        plt.errorbar(targ_qmeans, targ_umeans, xerr=targ_qmeans_err, yerr=targ_umeans_err,
                     color = plot_c, lw=0.75, fmt="o", alpha=0.9)
    else:
        plt.scatter(targ_qs, targ_us, 
                    color = plot_c, alpha=0.22)
        plt.errorbar(targ_qs, targ_us, 
                     xerr=targ_qstds, yerr=targ_ustds, 
                     lw=0.75, fmt="o", color=plot_c, alpha=0.22)
    if(pol_deg):
        for z in range(0, len(targ_qmeans)):
            plt.plot([0,targ_qmeans[z]], [0, targ_umeans[z]], 
                     'k-', lw=1.75, alpha=0.4, linestyle = '--')
    if(verbose_MJD_arg):
        for z in range(0, len(targ_qmeans)):
            plt.text(targ_qmeans[z], targ_umeans[z], 
                     int(t[z].mjd), rotation=45, fontsize=16)       
    else:
        for z in range(0, len(targ_qmeans)):
            plt.text(targ_qmeans[z], targ_umeans[z], 
                     str(t[z]), rotation=45, fontsize=16)
    
    plt.grid()
    plt.title("Pol Scatter")        
        
    plt.yticks(fontsize = 22)
    plt.xticks(fontsize = 22)        
    plt.ylabel('u', fontsize = 24)
    plt.xlabel('q', fontsize = 24)
    plt.axhline(y=0, color = 'black')
    plt.axvline(x=0, color = 'black')
    
    #It should take in data and error.
    if(to_excel):
        funcs_utils.data_to_excel((t.mjd, targ_qmeans, targ_qmeans_err), 'cal_eecep', 'Stokes_q')
        funcs_utils.data_to_excel((t.mjd, targ_umeans, targ_umeans_err), 'cal_eecep', 'Stokes_u')
            
    if(sv_im != ''):
        plt.title(" Pol Scatter")
        plt.savefig(sv_im,bbox_inches='tight',pad_inches=0.1 )
    
    plt.show()
    
    #Returns what what you plot... It strips.
    if(retrun_plotelems):
        return( targ_qmeans, targ_qmeans_err, targ_umeans, targ_umeans_err)
    
def q_n_u_stack_plot_v2( pol_data, sv_im_str ,pol_deg, launch_verb, key_verb):
    """
    #This is the version that takes in the generalised data format... the dictionary.
    """
    
    if(key_verb):
        print("Length of pol data:", len(pol_data))
    
    targ_date_strs = []
    zero_pol_date_strs = []
    high_pol_date_strs = []
    
    target_qs = []    #This is what is plotted
    target_us = []
    
    z_pol_qs = []    
    z_pol_us = []  
    
    h_pol_qs = []
    h_pol_us = []
    
    targ_qmeans = []
    targ_qstds = []
    targ_umeans = []
    targ_ustds = []
    
    z_pol_qmeans = []
    z_pol_qstds = []
    z_pol_umeans = []
    z_pol_ustds = []
    
    h_pol_qmeans = []
    h_pol_qstds = []
    h_pol_umeans = []
    h_pol_ustds = []
    
    simp_counter = 1
    for things in pol_data:
        print(simp_counter, list(things.keys())[0])#, ,things[list(things.keys())[0]], type(things[list(things.keys())[0]]))
        simp_counter+=1
        if('ee' in list(things.keys())[0] or 'EE'  in list(things.keys())[0]):
            targ_qmeans.append(np.mean(things[list(things.keys())[0]][0][1:])) #means
            targ_umeans.append(np.mean(things[list(things.keys())[0]][2][1:])) #ummeans
                                                                               #mean std
            targ_date_strs.append(list(things.keys())[0])
            
            target_qs = target_qs + things[list(things.keys())[0]][0][1:]
            targ_qstds = targ_qstds + things[list(things.keys())[0]][1][1:]
            target_us = target_us + things[list(things.keys())[0]][2][1:] 
            targ_ustds = targ_ustds + things[list(things.keys())[0]][3][1:]
            
        elif('g191' in list(things.keys())[0] or 'G191' in list(things.keys())[0] or 'hd212311' in list(things.keys())[0]):
            z_pol_qmeans.append(np.mean(things[list(things.keys())[0]][0][1:]))
            z_pol_umeans.append(np.mean(things[list(things.keys())[0]][2][1:]))
            zero_pol_date_strs.append(list(things.keys())[0])
            
            z_pol_qs = z_pol_qs + things[list(things.keys())[0]][0][1:]
            z_pol_qstds = z_pol_qstds + things[list(things.keys())[0]][1][1:]
            z_pol_us = z_pol_us + things[list(things.keys())[0]][2][1:]
            z_pol_ustds = z_pol_ustds + things[list(things.keys())[0]][3][1:]
        
        elif('215806' in list(things.keys())[0] or '287' in list(things.keys())[0] or '204827' in list(things.keys())[0] or '251204' in list(things.keys())[0] or '64106' in list(things.keys())[0]):
            h_pol_qmeans.append(np.mean(things[list(things.keys())[0]][0][1:]))
            h_pol_umeans.append(np.mean(things[list(things.keys())[0]][2][1:]))
            high_pol_date_strs.append(list(things.keys())[0])
            
            h_pol_qs = h_pol_qs + things[list(things.keys())[0]][0][1:]
            h_pol_qstds = h_pol_qstds + things[list(things.keys())[0]][1][1:]
            h_pol_us = h_pol_us + things[list(things.keys())[0]][2][1:]
            h_pol_ustds = h_pol_ustds + things[list(things.keys())[0]][3][1:]           
    
    if(launch_verb):
        fig, axs = plt.subplots(2, 2)
        fig.suptitle(" Pol Scatter ")

        axs[0, 0].scatter(target_qs, target_us, color = 'red', alpha=0.11)
        axs[0, 0].errorbar(target_qs, target_us, xerr=targ_qstds, yerr=targ_ustds, lw=0.75, fmt="o", color="r", alpha=0.1)

        axs[0, 0].grid()
        axs[0, 0].set_title('Target')
        axs[0, 0].tick_params(axis='both', which='major', labelsize=20)
        axs[0, 0].set_ylabel('u', fontsize = 24)
        axs[0, 0].set_xlabel('q', fontsize = 24)
        axs[0, 0].axhline(y=0, color = 'black')
        axs[0, 0].axvline(x=0, color = 'black')
        
        
        if(pol_deg):
            for z in range(0, len(targ_qmeans)):
                axs[0, 0].plot([0,targ_qmeans[z]], [0, targ_umeans[z]], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
                axs[0, 0].text(targ_qmeans[z], targ_umeans[z], targ_date_strs[z].replace("_","\n"), rotation=-45, fontsize=10)

        axs[1, 0].scatter(z_pol_qs,  z_pol_us, color = 'blue', alpha=0.11)
        axs[1, 0].errorbar(z_pol_qs,  z_pol_us, xerr=z_pol_qstds, yerr=z_pol_ustds, lw=0.75, fmt="o", color="blue", alpha=0.1)
        axs[1, 0].grid()
        axs[1, 0].set_title('Zero Polarization Standard')
        axs[1, 0].tick_params(axis='both', which='major', labelsize=20)
        axs[1, 0].set_ylabel('u', fontsize = 24)
        axs[1, 0].set_xlabel('q', fontsize = 24)
        axs[1, 0].axhline(y=0, color = 'black')
        axs[1, 0].axvline(x=0, color = 'black')
        if(pol_deg):
            for z in range(0, len(z_pol_qmeans)):
                axs[1, 0].plot([0,z_pol_qmeans[z]], [0, z_pol_umeans[z]], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
                axs[1, 0].text(z_pol_qmeans[z], z_pol_umeans[z], zero_pol_date_strs[z].replace("_","\n"), rotation=-45  , fontsize=10)

        axs[0, 1].scatter(h_pol_qs,  h_pol_us, color = 'green', alpha=0.11)
        axs[0, 1].errorbar(h_pol_qs,  h_pol_us, xerr=h_pol_qstds, yerr=h_pol_ustds, lw=0.75, fmt="o", color="green", alpha=0.1)        
        
        axs[0, 1].grid()
        axs[0, 1].set_title('High Polarization Standard')
        axs[0, 1].tick_params(axis='both', which='major', labelsize=20)
        axs[0, 1].set_ylabel('u', fontsize = 24)
        axs[0, 1].set_xlabel('q', fontsize = 24)
        axs[0, 1].axhline(y=0, color = 'black')
        axs[0, 1].axvline(x=0, color = 'black')
        if(pol_deg):
            for z in range(0, len(h_pol_qmeans)):
                axs[0, 1].plot([0,h_pol_qmeans[z]], [0, h_pol_umeans[z]], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
                axs[0, 1].text(h_pol_qmeans[z], h_pol_umeans[z], high_pol_date_strs[z].replace("_","\n"), rotation=-45  , fontsize=10)

        axs[1, 1].scatter(target_qs, target_us, color = 'red', alpha=0.11)
        axs[1, 1].errorbar(target_qs, target_us, xerr=targ_qstds, yerr=targ_ustds, lw=0.75, fmt="o", color="red", alpha=0.1)
        axs[1, 1].scatter(h_pol_qs,  h_pol_us, color = 'green', alpha=0.11)
        axs[1, 1].errorbar(h_pol_qs,  h_pol_us, xerr=h_pol_qstds, yerr=h_pol_ustds, lw=0.75, fmt="o", color="green", alpha=0.1)
        axs[1, 1].scatter(z_pol_qs,  z_pol_us, color = 'blue', alpha=0.11)
        axs[1, 1].errorbar(z_pol_qs,  z_pol_us, xerr=z_pol_qstds, yerr=z_pol_ustds, lw=0.75, fmt="o", color="blue", alpha=0.1)
        if(pol_deg):
            for z in range(0, len(targ_qmeans)):
                axs[1, 1].plot([0,targ_qmeans[z]], [0, targ_umeans[z]], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
            for z in range(0, len(z_pol_qmeans)):
                axs[1, 1].plot([0,z_pol_qmeans[z]], [0, z_pol_umeans[z]], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
            for z in range(0, len(h_pol_qmeans)):
                axs[1, 1].plot([0,h_pol_qmeans[z]], [0, h_pol_umeans[z]], 'k-', lw=1.75, alpha=0.4, linestyle = '--')

        axs[1, 1].grid()
        axs[1, 1].set_title('All objects combined')
        axs[1, 1].tick_params(axis='both', which='major', labelsize=20)
        axs[1, 1].set_ylabel('u', fontsize = 24)
        axs[1, 1].set_xlabel('q', fontsize = 24)
        axs[1, 1].axhline(y=0, color = 'black')
        axs[1, 1].axvline(x=0, color = 'black')

        fig.tight_layout()
        
        if(sv_im_str):
            plt.savefig(sv_im_str,bbox_inches='tight',
                            pad_inches=0.1)

def q_n_u_stack_plot(target_data, zero_pol_std, high_pol_std, MJD_obs , name_array ,title , pol_deg, sv_im ):
    """
    #Function that scatter plots q nad u
    #This thing just plots. It does not do anything fance such as compute mean blah
    #Error should be taken in quadrature
    
    #with all the shit going into this function I have a feeling that it will easily break
    
    #Plot all in 4
    """
    #print("Plot all 4!")

    combined_q = target_data[0][1:] + zero_pol_std[0][1:] +  high_pol_std[0][1:]
    combined_u = target_data[2][1:] + zero_pol_std[2][1:] +  high_pol_std[2][1:]
    
    #I already have it haha! Simple!
    #This usesless in the context of how it is used.
    #So go ahead and use std again
    targ_mean_q, targ_mean_u, median_q, median_u, targ_q_std, targ_u_std = q_u_stats(target_data[0][1:] ,target_data[1][1:] ,target_data[2][1:] ,target_data[3][1:])
    zero_pol_mean_q, zero_pol_mean_u, median_q, median_u, zero_pol_q_std, zero_pol_u_std = q_u_stats(zero_pol_std[0][1:] ,zero_pol_std[1][1:] ,zero_pol_std[2][1:] ,zero_pol_std[3][1:])
    high_pol_mean_q, high_pol_mean_u, median_q, median_u, high_pol_q_std, high_pol_u_std = q_u_stats(high_pol_std[0][1:] ,high_pol_std[1][1:] ,high_pol_std[2][1:] ,high_pol_std[3][1:])

    fig, axs = plt.subplots(2, 2)
    fig.suptitle(MJD_obs+" Pol Scatter "+title)
    
    axs[0, 0].scatter(target_data[0][1:], target_data[2][1:], color = 'red', alpha=0.1)
    axs[0, 0].errorbar(target_data[0][1:],  target_data[2][1:], xerr=target_data[1][1:], yerr=target_data[3][1:], lw=0.75, fmt="o", color="r", alpha=0.1)
    
    axs[0, 0].scatter([targ_mean_q], [targ_mean_u], color = 'red', alpha=0.98)
    axs[0, 0].errorbar([targ_mean_q], [targ_mean_u], xerr=[targ_q_std], yerr=[targ_u_std], lw=0.75, fmt="o", color="r", alpha=0.98)
    
    #Here
    axs[0, 0].text(targ_mean_q, targ_mean_u, name_array[0]+ "\nq:"+str(np.round(targ_mean_q, 4))+u"\u00B1"+ str(np.round(targ_q_std, 3))+"\nu:"+str(np.round(targ_mean_u,4))+u"\u00B1"+ str(np.round(targ_u_std, 3)), fontsize=20    )
    
    axs[0, 0].set_title("Target")
    axs[0, 0].grid()
    axs[0, 0].axvline(x=0, lw = 1, color = 'black', alpha=0.6)
    axs[0, 0].axhline(y=0, lw = 1, color = 'black', alpha=0.6)
    if(pol_deg):
        axs[0, 0].plot([0,targ_mean_q], [0, targ_mean_u], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
        
    axs[1, 0].scatter(zero_pol_std[0][1:], zero_pol_std[2][1:], color = 'blue', alpha=0.1)
    axs[1, 0].errorbar(zero_pol_std[0][1:],  zero_pol_std[2][1:], xerr=zero_pol_std[1][1:], yerr=zero_pol_std[3][1:], lw=0.75, fmt="o", color="blue", alpha=0.1)
    
    axs[1, 0].scatter([zero_pol_mean_q], [zero_pol_mean_u], color = 'blue', alpha=0.98)
    axs[1, 0].errorbar([zero_pol_mean_q], [zero_pol_mean_u], xerr=[zero_pol_q_std], yerr=[zero_pol_u_std], lw=0.75, fmt="o", color="b", alpha=0.98)
    
    #Here
    axs[1, 0].text(zero_pol_mean_q, zero_pol_mean_u, name_array[1]+ "\nq:"+str(np.round(zero_pol_mean_q, 4))+u"\u00B1"+ str(np.round(zero_pol_q_std, 3))+"\nu:"+str(np.round(zero_pol_mean_u,4))+u"\u00B1"+ str(np.round(zero_pol_u_std, 3)), fontsize=20)
    
    axs[1, 0].set_title("Zero Polarization Standard Star")
    if(pol_deg):
        axs[1, 0].plot([0,zero_pol_mean_q], [0, zero_pol_mean_u], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
    axs[1, 0].axvline(x=0, lw = 1, color = 'black', alpha=0.6)
    axs[1, 0].axhline(y=0, lw = 1, color = 'black', alpha=0.6)
    axs[1, 0].grid()
    
    axs[0, 1].scatter(high_pol_std[0][1:], high_pol_std[2][1:], color = 'green', alpha=0.11)
    axs[0, 1].errorbar(high_pol_std[0][1:],  high_pol_std[2][1:], xerr=high_pol_std[1][1:], yerr=high_pol_std[3][1:], lw=0.75, fmt="o", color="green", alpha=0.11)
    
    axs[0, 1].scatter([high_pol_mean_q], [high_pol_mean_u], color = 'green', alpha=0.98)
    axs[0, 1].errorbar([high_pol_mean_q], [high_pol_mean_u], xerr=[high_pol_q_std], yerr=[high_pol_u_std], lw=0.75, fmt="o", color="g", alpha=0.98)
    
    #Here
    axs[0, 1].text(high_pol_mean_q, high_pol_mean_u, name_array[2]+ "\nq:"+str(np.round(high_pol_mean_q, 4))+u"\u00B1"+ str(np.round(high_pol_q_std, 3))+"\nu:"+str(np.round(high_pol_mean_u,4))+u"\u00B1"+ str(np.round(high_pol_u_std, 3)), fontsize=20)
    
    axs[0, 1].set_title("High Polarization Standard Star")
    
    #axs[1, 1].text(zero_pol_mean_q, zero_pol_mean_u, str(np.round(zero_pol_mean_q, 4))+u"\u00B1"+ str(np.round(zero_pol_q_std, 3))+","+str(np.round(zero_pol_mean_u,4))+u"\u00B1"+ str(np.round(zero_pol_u_std, 3)), fontsize=20)
    
    if(pol_deg):
        axs[0, 1].plot([0,high_pol_mean_q], [0, high_pol_mean_u], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
    axs[0, 1].axvline(x=0, lw = 1, color = 'black', alpha=0.6)
    axs[0, 1].axhline(y=0, lw = 1, color = 'black', alpha=0.6)
    axs[0, 1].grid()
    
    axs[1, 1].scatter(target_data[0][1:], target_data[2][1:], color = 'red', alpha=0.11)
    axs[1, 1].errorbar(target_data[0][1:],  target_data[2][1:], xerr=target_data[1][1:], yerr=target_data[3][1:], lw=0.75, fmt="o", color="red", alpha=0.11)
    axs[1, 1].scatter([targ_mean_q], [targ_mean_u], color = 'red', alpha=0.98)
    axs[1, 1].errorbar([targ_mean_q], [targ_mean_u], xerr=[targ_q_std], yerr=[targ_u_std], lw=0.75, fmt="o", color="r", alpha=0.98)
    
    #Here
    axs[1, 1].text(targ_mean_q, targ_mean_u, name_array[0]+ "\nq:"+str(np.round(targ_mean_q, 4))+u"\u00B1"+ str(np.round(targ_q_std, 3))+"\nu:"+str(np.round(targ_mean_u,4))+u"\u00B1"+ str(np.round(targ_u_std, 3)), fontsize=20    )
    
    if(pol_deg):
        axs[1, 1].plot([0,targ_mean_q], [0, targ_mean_u], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
        
    axs[1, 1].scatter(zero_pol_std[0][1:], zero_pol_std[2][1:], color = 'blue', alpha=0.11)
    axs[1, 1].errorbar(zero_pol_std[0][1:],  zero_pol_std[2][1:], xerr=zero_pol_std[1][1:], yerr=zero_pol_std[3][1:], lw=0.75, fmt="o", color="blue", alpha=0.11)
    axs[1, 1].scatter([zero_pol_mean_q], [zero_pol_mean_u], color = 'blue', alpha=0.98)
    axs[1, 1].errorbar([zero_pol_mean_q], [zero_pol_mean_u], xerr=[zero_pol_q_std], yerr=[zero_pol_u_std], lw=0.75, fmt="o", color="b", alpha=0.98)
    
    axs[1, 1].text(zero_pol_mean_q, zero_pol_mean_u, name_array[1]+ "\nq:"+str(np.round(zero_pol_mean_q, 4))+u"\u00B1"+ str(np.round(zero_pol_q_std, 3))+"\nu:"+str(np.round(zero_pol_mean_u,4))+u"\u00B1"+ str(np.round(zero_pol_u_std, 3)), fontsize=20)
    
    if(pol_deg):
        axs[1, 1].plot([0,zero_pol_mean_q], [0, zero_pol_mean_u], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
    axs[1, 1].scatter(high_pol_std[0][1:], high_pol_std[2][1:], color = 'green', alpha=0.11)
    axs[1, 1].errorbar(high_pol_std[0][1:],  high_pol_std[2][1:], xerr=high_pol_std[1][1:], yerr=high_pol_std[3][1:], lw=0.75, fmt="o", color="green", alpha=0.11)
    
    axs[1, 1].scatter([high_pol_mean_q], [high_pol_mean_u], color = 'green', alpha=0.98)
    axs[1, 1].errorbar([high_pol_mean_q], [high_pol_mean_u], xerr=[high_pol_q_std], yerr=[high_pol_u_std], lw=0.75, fmt="o", color="g", alpha=0.98)
    
    axs[1, 1].text(high_pol_mean_q, high_pol_mean_u, name_array[2]+ "\nq:"+str(np.round(high_pol_mean_q, 4))+u"\u00B1"+ str(np.round(high_pol_q_std, 3))+"\nu:"+str(np.round(high_pol_mean_u,4))+u"\u00B1"+ str(np.round(high_pol_u_std, 3)), fontsize=20)
    
    if(pol_deg):
        axs[1, 1].plot([0,high_pol_mean_q], [0, high_pol_mean_u], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
    axs[1, 1].axvline(x=0, lw = 1, color = 'black', alpha=0.6)
    axs[1, 1].axhline(y=0, lw = 1, color = 'black', alpha=0.6)
        
    axs[1, 1].set_title("Combined")
    axs[1, 1].grid()
    fig.tight_layout()
    
    if(sv_im):
        plt.savefig(MJD_obs+" pol_scatter "+title,bbox_inches='tight',
                        pad_inches=0.1)        
        
def calib_data(inp_data, 
               instrumental_pol, 
               plt_show = False,
               verbose=False):
    """
    A function that takes in data and calculates position angle (PA). Returns data to be used in other plot tools. More presentable publication quality style.

    Parameters
    ----------
    input_data : tuple
        Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
    instrumental_pol : tuple
        Tuple of list
    plot_show : bool, optional
        Prints calculations for extra verbosity.  False by default 
    verbose : bool, optional
        Saves image to file.  False by default 
    """  
    cal_prod = cp(inp_data) #I have copied the data
    
    print("Calibrating data...") #calibration point an
    
    if(verbose):
        print("Data (pre cal):", inp_data)
        print("Instrumental Polarization:", instrumental_pol) 
        
    for k in range(0, len(cal_prod[0])):    #runs through the calib product    
        for l in range(0, len(cal_prod[0][k][ list(cal_prod[0][k].keys())[0]][0][1:])):
            if(verbose): #goes through each q value and shifts them
                print("Old val:", cal_prod[0][ list(cal_prod[0][k].keys())[0]][0][l])

            cal_prod[0][k][list(cal_prod[0][k].keys())[0]][0][l] = cal_prod[0][k][ list(cal_prod[0][k].keys())[0]][0][l] - instrumental_pol[0][0]

            if(verbose):            
                print("New val:", cal_prod[0][k][list(cal_prod[0][k].keys())[0]][0][l])

        for l in range(0, len(cal_prod[0][k][list(cal_prod[0][k].keys())[0]][2][1:])):
            if(verbose):
                print("Old val:", cal_prod[0][k][ list(cal_prod[0][k].keys())[0]][0][l])
                
            cal_prod[0][k][list(cal_prod[0][k].keys())[0]][2][l] = cal_prod[0][k][ list(cal_prod[0][k].keys())[0]][2][l] - instrumental_pol[1][0]
            if(verbose):
                print("New val:", cal_prod[k][ list(cal_prod[k].keys())[0]][0][l])
    
    return(cal_prod)
    
    """
    for dats in inp_data:
        for h in range(0, len(   dats[list(dats.keys())[0]][0][1:]   )):
            if(verbose):
                print("qs", (dats[list(dats.keys())[0]][0][1:][h]) )
                print("qs cal:",  (dats[list(dats.keys())[0]][0][1:][h])  - instrumental_pol[0])
                print("us", (dats[list(dats.keys())[0]][2][1:][h]))
                print("us cal", (dats[list(dats.keys())[0]][2][1:][h]) - instrumental_pol[2])
            
            qs_cal.append((dats[list(dats.keys())[0]][0][1:][h])  - instrumental_pol[0])
            us_cal.append((dats[list(dats.keys())[0]][2][1:][h])  - instrumental_pol[2])
            
            #overwrite the points
            
        qs_uncal += dats[list(dats.keys())[0]][0][1:]
        us_uncal += dats[list(dats.keys())[0]][2][1:]
        
        qs_err += dats[list(dats.keys())[0]][1][1:]
        us_err += dats[list(dats.keys())[0]][3][1:]

     
    if(plt_show):
        plt.scatter( qs_uncal , us_uncal) #orange 
        plt.scatter( qs_cal, us_cal)   #blue

        plt.scatter(qs_uncal , us_uncal)#, lw=0.75, fmt="o", alpha=0.9)
        plt.scatter(qs_cal , us_cal)
        #plt.errorbar(qs_cal, us_cal, xerr=qs_err, yerr=us_err, lw=0.75, fmt="o", alpha=0.9)
        plt.title("Calibrated Pol Scatter")        
        plt.yticks(fontsize = 22)
        plt.xticks(fontsize = 22)        
        plt.ylabel('u', fontsize = 24)
        plt.xlabel('q', fontsize = 24)
        plt.axhline(y=0, color = 'black')
        plt.axvline(x=0, color = 'black')
        plt.grid()
        plt.show()

        plt.errorbar(qs_uncal , us_uncal, xerr=qs_err, yerr=us_err, lw=0.75, fmt="o", alpha=0.9)
        plt.errorbar(qs_cal, us_cal, xerr=qs_err, yerr=us_err, lw=0.75, fmt="o", alpha=0.9)
        plt.title("Calibrated Pol Scatter")        
        plt.yticks(fontsize = 22)
        plt.xticks(fontsize = 22)        
        plt.ylabel('u', fontsize = 24)
        plt.xlabel('q', fontsize = 24)
        plt.axhline(y=0, color = 'black')
        plt.axvline(x=0, color = 'black')
        plt.grid()
        plt.show()

        fig, ax = plt.subplots()
        markers, caps, bars = ax.errorbar(qs_uncal , us_uncal, xerr=qs_err, yerr=us_err,
                fmt='o', ecolor='blue',capsize=1, capthick=1)
        markers, caps, bars = ax.errorbar(qs_cal, us_cal, xerr=qs_err, yerr=us_err,
                fmt='o', ecolor='orange',capsize=1, capthick=1)
        plt.title("Calibrated Pol Scatter", fontsize=24)
        plt.yticks(fontsize = 22)
        plt.xticks(fontsize = 22)
        plt.axhline(y=0, color = 'black')
        plt.axvline(x=0, color = 'black')

        plt.grid()
        [bar.set_alpha(0.1) for bar in bars]
        [cap.set_alpha(0.95) for cap in caps]
        plt.show()
        #if(sv_im != ''):
        #    plt.savefig(sv_im,bbox_inches='tight',pad_inches=0.1 )
    """            
def mean_q_u_check(inp_data, n, q_u_ret, verb_arg):
    q_top = []
    q_bot = []

    u_top = []
    u_bot = []

    q_mean_alt = np.mean(inp_data[n][list(inp_data[n].keys())[0]][0][1:])
    u_mean_alt = np.mean(inp_data[n][list(inp_data[n].keys())[0]][2][1:]) 

    qs = inp_data[n][list(inp_data[n].keys())[0]][0][1:]
    q_errs = inp_data[n][list(inp_data[n].keys())[0]][1][1:]
    us = inp_data[n][list(inp_data[n].keys())[0]][2][1:]
    u_errs = inp_data[n][list(inp_data[n].keys())[0]][3][1:]

    for l in range(0, len(qs)):
        top = qs[l]/(q_errs[l]**2)
        q_top.append(top)

        top = us[l]/(u_errs[l]**2)
        u_top.append(top)

        bot = (1/q_errs[l]**2)
        q_bot.append(bot)

        bot = (1/u_errs[l]**2)
        u_bot.append(bot)

    q_mean = sum(q_top)/sum(q_bot)
    u_mean = sum(u_top)/sum(u_bot)

    if(verb_arg):
        print("q_mean, u_mean:",q_mean , u_mean)
        print("q_alt, u_alt:", q_mean_alt, u_mean_alt)
        print("absolute difference:", abs(q_mean_alt-q_mean), abs(u_mean_alt-u_mean))
    
    if(q_u_ret == 'mean'):
        return (q_mean_alt, u_mean_alt)
    elif(q_u_ret == 'RINGO'):
        return (q_mean , u_mean)
    else:
        print("Invalid")
        
