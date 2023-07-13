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

plt.rcParams['savefig.facecolor'] = 'black'
def check_pol_std():
    """
    A function that prints out high polarization and zero polarization standards. Takes in no input arguments. Please double check this list because some of them are targets. Does not return anything.

    Parameters
    ----------
        None.
    """
    #perhaps we could scrape this!
    #https://www.not.iac.es/instruments/turpol/std/hpstd.html
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
    
    #perhaps we could scrape this!
    #https://www.not.iac.es/instruments/turpol/std/zpstd.html
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
             plot_legend=[],
             plot_title='',
             sv_pold_img='',
             perc_arg= False, 
             verbose_calc_pd=False, 
             verbose_data=False
             ):
    """
    A function that takes in data and calculates polarization degree (PD). Does not return anything. Can be used in preliminary analysis of polarization data.

    Parameters
    ----------
        input_data : List
            List of standardised polarimetric data containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
        plot_legend : List,
            List of strings that define the legends on the plot. Manually input the target name.
        plot_title : str,
            String that defines the title of the plot. Used in the filename if the pol data is saved to file.
        perc_arg : bool, optional
            Multiply calculated pd values by 100 to get percentage polarization.
        verbose_calc_pd : bool, optional
            Prints calculations for extra verbosity.  False by default.
        verbose_data : bool, optional
            Prints calculations for extra verbosity.  False by default.
        sv_arg : bool, optional
            Saves image to file.  False by default.
    """   
    print("Calculate and plot polarization degree for", len(input_data), "inputs without returning data.")
    for k in range(0, len(input_data)): 
        print("Data", k+1, "has length", len(input_data[k][0]) )
        
    pol_d_array = []
    pol_d_err_array = []
    t_array = []
    
    for data2d in input_data:
        pol_d = []
        pol_d_err = []
        
        for dats in data2d[0]:
            result = re.search('(.*)_', list(dats.keys())[0][:12])
            targ_name_str = re.search('_(.*)', list(dats.keys())[0])

            ###Math
            mean_q_squared = np.mean(dats[list(dats.keys())[0]][0][1:])**2 #mean q squared
            mean_u_squared = np.mean(dats[list(dats.keys())[0]][2][1:])**2 #mean u squared

            mean_q_err_squared = np.std(dats[list(dats.keys())[0]][1][1:])**2 #std q squared
            mean_u_err_squared = np.std(dats[list(dats.keys())[0]][3][1:])**2 #std u squared

            sum_o_squares = mean_q_squared + mean_u_squared #mean q squared + mean u squared

            pol_dval = math.sqrt(sum_o_squares)
            pol_d_errval  = math.sqrt(((mean_q_squared)/(sum_o_squares))*(mean_q_err_squared)+ ((mean_u_squared)/(sum_o_squares))*(mean_u_err_squared))
            ###
            
            if(perc_arg):
                if(verbose_calc_pd):
                    print(targ_name_str.group(1), "MJD:", result.group(1), 
                          pol_dval*100, u"\u00B1",  pol_d_errval*100)
                pol_d.append(pol_dval*100)
                pol_d_err.append(pol_d_errval*100)
            else:
                if(verbose_calc_pd):
                    print("MJD:", result.group(1), 
                          pol_dval, u"\u00B1",  pol_d_errval)
                pol_d.append(pol_dval)
                pol_d_err.append(pol_d_errval)

        pol_d_array.append(pol_d)
        pol_d_err_array.append(pol_d_err)
                
        t = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in data2d[1]], format='isot', scale='utc')
        t_array.append(t)    
    
    #Start of plot part
    color_arr = ['red', 'blue', 'green', 'purple', 'magenta' , 'orange', 'aquamarine']
    for l in range(0, len(input_data)):
        plt.errorbar(t_array[l].mjd , 
                 pol_d_array[l], 
                 xerr=[0]*len(pol_d_array[l]), 
                 yerr=pol_d_err_array[l], 
                 lw=0.75, fmt="^", capsize=10,
                 color=color_arr[l], alpha=0.9)

    plt.title(plot_title)
    plt.legend(plot_legend)
    plt.grid()
    
    if(verbose_data):
        for k in range(0, len(t_array)):
            for i in range(0, len(t_array[k])):
                plt.text(t_array[k][i].mjd, 
                         pol_d_array[k][i], 
                         str(np.round(t_array[k][i].mjd))+"\nPD:"+str(np.round(pol_d_array[k][i],2)),
                         fontsize=12, 
                         rotation=45)
                
    plt.xlabel("time (MJD)")
    if(perc_arg):
        plt.ylabel("PD %")
    else:
        plt.ylabel("PD")      
            
    if(sv_pold_img != ''):
        plt.savefig(sv_pold_img,bbox_inches='tight',pad_inches=0.1 , facecolor='w' )        
    plt.show()
    
def calc_pa2(input_data,
             plot_legend=[],
             plot_title='',
             deg_arg=False, 
             verbose_calc_pa=False, 
             verbose_data=False,
             sv_polpa_img=False):
    """
    A function that takes in data and calculates position angle (PA). Does not return anything. Can be used in preliminary analysis of polarization data.

    Parameters
    ----------
        input_data : List
            List of standardised polarimetric data containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
        plot_legend : List,
            List of strings that define the legends on the plot. Manually input the target name.
        plot_title : str,
            String that defines the title of the plot. Used in the filename if the pol data is saved to file.
        deg_arg : bool, optional
            Multiply calculated pd values by 180/pi to get pol angle degree. 
        verbose_calc_pa : bool, optional
            Prints calculations for extra verbosity.  False by default
        verbose_data : bool, optional
            Prints calculations for extra verbosity.  False by default 
        sv_polpa_img : bool, optional
            Saves image to file.  False by default 
    """
    
    print("Calculate and plot polarization degree for", len(input_data), "input types without returning data")
    for k in range(0, len(input_data)): 
        print("Data", k+1, "has length", len(input_data[k][0]) )
    
    pol_pa_array = []
    pol_pa_err_array = []
    t_array = []
    
    for data2d in input_data:
        pol_pa = []
        pol_pa_err = []
    
        for dats in data2d[0]:        
            result = re.search('(.*)_', list(dats.keys())[0][:12])
            targ_name_str = re.search('_(.*)', list(dats.keys())[0] )

            ###Math
            mean_q_err_squared = np.std(dats[list(dats.keys())[0]][1][1:])**2
            mean_u_err_squared = np.std(dats[list(dats.keys())[0]][3][1:])**2

            mean_q_squared = np.mean(dats[list(dats.keys())[0]][0][1:])**2
            mean_u_squared = np.mean(dats[list(dats.keys())[0]][2][1:])**2

            mean_q = np.mean(dats[list(dats.keys())[0]][0][1:])
            mean_u = np.mean(dats[list(dats.keys())[0]][2][1:])

            sum_o_squares = mean_q_squared + mean_u_squared

            pol_paval = 0.5*math.atan2(mean_u , mean_q)
            pol_pa_errval = math.sqrt(((1/(2*mean_q*(1 + (mean_u_squared/mean_q_squared))))**2 )*(mean_u_err_squared) + ((-1*((mean_u)/(2*sum_o_squares)))**2 )*(mean_q_err_squared))
            ###

            if(deg_arg):
                if(verbose_calc_pa): #
                    print("MJD:", result.group(1), 
                          pol_paval*(180/3.142), u"\u00B1",  pol_pa_errval*(180/3.142))
                pol_pa.append(pol_paval*(180/3.142)) 
                pol_pa_err.append(pol_pa_errval*(180/3.142)) 
            else:
                if(verbose_calc_pa):
                    print("MJD:", result.group(1), pol_paval, u"\u00B1",  pol_pa_errval)
                pol_pa.append(pol_paval)
                pol_pa_err.append(pol_pa_errval)

        pol_pa_array.append(pol_pa)
        pol_pa_err_array.append(pol_pa_err)
        
        t = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in data2d[1]], format='isot', scale='utc')
        t_array.append(t)
        
    color_arr = ['red', 'blue', 'green', 'purple', 'magenta' , 'orange', 'aquamarine']
    for l in range(0, len(input_data)):
        plt.errorbar(t_array[l].mjd , 
                 pol_pa_array[l], 
                 xerr=[0]*len(pol_pa_array[l]), 
                 yerr=pol_pa_err_array[l], 
                 lw=0.75, fmt="^", capsize=10,
                 color=color_arr[l], alpha=0.9)
            
    plt.title(plot_title)
    plt.legend(plot_legend)
    plt.grid()
    
    if(verbose_data):
        for k in range(0, len(t_array)):
            for i in range(0, len(t_array[k])):
                plt.text(t_array[k][i].mjd, 
                         pol_pa_array[k][i], 
                         str(np.round(t_array[k][i].mjd))+"\nPA:"+str(np.round(pol_pa_array[k][i],2)),
                         fontsize=12, rotation=45)
                
    plt.xlabel("time (MJD)")
    if(deg_arg):
        plt.ylabel("PA " + '(\u00b0)')
    else:
        plt.ylabel("PA rad")  
        
    if(sv_polpa_img):
        plt.savefig(sv_im_str, bbox_inches='tight',pad_inches=0.1)
        
    plt.show()  

def calc_PD_stability(input_data,              
                      targ_corr_MJD='',
                      perc_arg=False,
                      to_excel=False, 
                      corr_MJD=False,
                      verbose_calc_pd=False, 
                      verbose_mjd_align_check=False,):
    
    """
    A function that takes in data and calculates polarization degree (PD). Returns data as list to be used in other plot tools.

    Parameters
    ----------
        input_data : tuple
            Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
        targ_corr_MJD : str
            String to identify filename for correcting MJD. Used inconjunctuiib if corr_MJD=True
        corr_MJD : bool, optional
            Corrects MJD if MJD needs correcting
        perc_arg : bool, optional
            Multiply calculated pd values by 100 to get percentage polarization. 
        to_excel : bool, optional
            Saves data to file.  False by default 
        verbose_calc_pd : bool, optional
            Multiply calculated pd values by 100 to get percentage polarization degree. 
        verbose_mjd_align_check : bool, optional
            Prints calculations for extra verbosity.  False by default 
    """
    dates = sorted(input_data[1])
    objn = list(input_data[0][0].keys())[0][11:]
    
    print("Calculate and return polarization degree (PD) for duration", dates[0], "to", dates[-1], "\nfor", objn)
    
    if(corr_MJD):
        if(targ_corr_MJD=='EECEP'):
            with open('./data_pickles/MJD_corr_pix.pickle', 'rb') as fid: #
                corr_da = pickle.load(fid)
        elif(targ_corr_MJD=='bd64'):
            with open('./data_pickles/MJD_corr_64106.pickle', 'rb') as fid: #
                corr_da = pickle.load(fid)
        elif(targ_corr_MJD=='g191'):
            with open('./data_pickles/MJD_corr_G191.pickle', 'rb') as fid: #
                corr_da = pickle.load(fid)            
        elif(targ_corr_MJD=='215806'):
            with open('./data_pickles/MJD_corr_hd215806.pickle', 'rb') as fid: #
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
        ###
        m_q_sq = np.mean(dats[list(dats.keys())[0]][0][1:])**2 #mean q squared
        m_u_sq = np.mean(dats[list(dats.keys())[0]][2][1:])**2 #mean u squared 

        std_q_err_sq = np.std(dats[list(dats.keys())[0]][0][1:])**2 #std q squared
        std_u_err_sq = np.std(dats[list(dats.keys())[0]][2][1:])**2 #std u squared
        
        sum_sq = m_q_sq + m_u_sq #mean q squared plus mean u squared
                
        pol_d = math.sqrt(sum_sq)        
        pol_d_err = math.sqrt((std_q_err_sq)*((m_q_sq)/(sum_sq)) + ((std_u_err_sq)*((m_u_sq)/(sum_sq))))
        ###
        
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
                      corr_MJD=False,
                      deg_arg=False,
                      to_excel=False, 
                      PA_shift=False,
                      verbose_calc_pa=False,
                      verbose_mjd_align_check=False,):
    """
    A function that takes in data and calculates position angle (PA). Returns data to be used in other plot tools

    Parameters
    ----------
        input_data : tuple
            Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
        targ_corr_MJD : str
            String to identify filename for correcting MJD. Used inconjunctuiib if corr_MJD=True
        corr_MJD : bool, optional
            Corrects MJD if MJD needs correcting
        deg_arg : bool, optional
            Multiply calculated pd values by 180/pi to get pol angle degree. 
        to_excel : bool, optional
            Saves image to file.  False by default 
        verbose_calc_pa : bool, optional
            Multiply calculated pd values by 180/pi to get position angle in degree degree. 
        verbose_mjd_align_check : bool, optional
            Prints calculations and time for extra verbosity. False by default 
    """
    dates = sorted(input_data[1])
    objn = list(input_data[0][0].keys())[0][11:]
    
    print("Calculate and return position angle (PA) for duration", dates[0], "to", dates[-1], "\nfor", objn)
    
    if( corr_MJD ):
        if(targ_corr_MJD=='EECEP'):
            with open('./data_pickles/MJD_corr_pix.pickle', 'rb') as fid: 
                corr_da = pickle.load(fid)
        elif(targ_corr_MJD=='bd64'):
            with open('./data_pickles/MJD_corr_64106.pickle', 'rb') as fid: 
                corr_da = pickle.load(fid)
        elif(targ_corr_MJD=='g191'):
            with open('./data_pickles/MJD_corr_G191.pickle', 'rb') as fid: 
                corr_da = pickle.load(fid)            
        elif(targ_corr_MJD=='215806'):
            with open('./data_pickles/MJD_corr_hd215806.pickle', 'rb') as fid: 
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
        
        ###
        mean_q_err_squared = np.std(dats[list(dats.keys())[0]][0][1:])**2 #taking 
        mean_u_err_squared = np.std(dats[list(dats.keys())[0]][2][1:])**2
        
        mean_q_squared = np.mean(dats[list(dats.keys())[0]][0][1:])**2
        mean_u_squared = np.mean(dats[list(dats.keys())[0]][2][1:])**2
        
        mean_q = np.mean(dats[list(dats.keys())[0]][0][1:])
        mean_u = np.mean(dats[list(dats.keys())[0]][2][1:])
        
        sum_o_squares = mean_q_squared + mean_u_squared
        
        pol_pa = 0.5*math.atan2(mean_u , mean_q )
        pol_pa_err = math.sqrt(((1/(2*mean_q*(1 + (mean_u_squared/mean_q_squared))))**2 )*(mean_u_err_squared) + ((-1*( (mean_u)/(2*sum_o_squares)))**2 )*(mean_q_err_squared))
        ###
        
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
    A function that takes in computed data and calculates position angle (PA). Returns data to be used in other plot tools. More presentable publication quality style.

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
    
def plot_pol_stab_PDPA(input_data, 
                         fmt_icon, 
                         color_icon):
    
    """
    A function that takes in data and calculates position angle (PA) with double axis left and right. Returns data to be used in other plot tools. More presentable publication quality style.

    Parameters
    ----------
        input_data : list
            List of data where expectation is:
                input_data[0] = MJD_track_PD
                input_data[1] = MJD_track_PA
                input_data[2] = obj_pol_PD
                input_data[3] = obj_pol_PA
                input_data[4] = obj_pol_PD_err
        icon : list
            List of strings where expectation is:
                fmt_icon[0] = '*'
                fmt_icon[1] = '-'
            Or whatever format icon is preferred/ accepted
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
                       in_title = '',
                       sv_im='', 
                       plot_verbose=False, 
                       verbose=False, 
                       m_plot=False):
    """
    A function that takes in data and calculates position angle (PA). Returns data to be used in other plot tools. More presentable publication quality style. Returns a list  containing the mean and the mean error.

    Parameters
    ----------
        input_data : tuple
            Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps).
        q_u_check : str
            String to identify data file for correcting MJD
        sv_im : bool, optional
            Saves image to file.  False by default  
        plot_verbose : bool, optional
            Prints calculations for extra verbosity.  False by default 
        m_plot : bool, optional
            Makes plot
    """
    
    means_arr = []
    means_err_arr = []
    c = 0
    
    for dats in input_data[0]:        
        result = re.search('(.*)_', list(dats.keys())[0][:12])
        targ_name_str = re.search('_(.*)', list(dats.keys())[0])
                
        if(q_u_check=='q'):
            if(verbose):
                print(c, np.mean(dats[list(dats.keys())[0]][0][1:]),u"\u00B1",np.std(dats[list(dats.keys())[0]][0][1:]))
            means_arr.append(np.mean(dats[list(dats.keys())[0]][0][1:]))
            means_err_arr.append(np.std(dats[list(dats.keys())[0]][0][1:]))
        elif(q_u_check=='u'):
            if(verbose):
                print(c, np.mean(dats[list(dats.keys())[0]][2][1:]),u"\u00B1",np.std(dats[list(dats.keys())[0]][2][1:]))
            means_arr.append(np.mean(dats[list(dats.keys())[0]][2][1:]))
            means_err_arr.append(np.std(dats[list(dats.keys())[0]][2][1:]))
        c+=1
            
    t = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in input_data[1]], format='isot', scale='utc')
    
    #use this to check for weird things
    #print(means_arr  )
    #print(means_err_arr  )
    #means_arr = [x for x in means_arr if str(x) != 'nan']
    #means_err_arr= [x for x in  means_err_arr if str(x) != 'nan']
    #for k in range(0, len(means_arr)):
    #    print(k, means_arr[k], means_err_arr[k], t[k])
    
    if(m_plot):
        fig, ax = plt.subplots()
        markers, caps, bars = ax.errorbar(t.mjd, means_arr, yerr=means_err_arr, xerr =[0]*len(means_arr),
                fmt='o', ecolor='blue',capsize=2, capthick=2)

        plt.title(in_title+ ' ' + q_u_check + " stability over time", fontsize=24)
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
            
        #plt.axhline(y=np.mean(means_arr  ), color = 'black')
        plt.axhline(y=0, color = 'black')
        plt.xlabel('MJD',fontsize=24)
        plt.grid()

        [bar.set_alpha(0.22) for bar in bars]
        [cap.set_alpha(0.5) for cap in caps]

        if(sv_im != ''):
            plt.savefig(sv_im,bbox_inches='tight',pad_inches=0.1 , facecolor='w')
        plt.show()   
        
    return ([means_arr, means_err_arr]) #returns the plot here
    """
    """        
def q_n_u_single_plot_v1(input_data,
                         plot_c='blue',  
                         in_title = '',
                         sv_im='', 
                         reflection_axis =  '',
                         time_ord = False,
                         verbose_MJD_arg=False,
                         pol_deg=False, 
                         only_means=False, 
                         to_excel=False, 
                         retrun_plotelems=False,
                         key_verb=False):
    """
    A function that takes in data and plots q and u scatter. Returns an array containing mean q, mean q error, mean u, mean u error

    Parameters
    ----------
        input_data : tuple
            Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
        plot_c : str
            String to specify the color of the plot. Blue by default.
        sv_im : bool, optional
            Multiply calculated pd values by 180/pi to get position angle in degree degree. 
        verbose_MJD_arg : bool, optional
            Prints MJDs 
        pol_deg : bool, optional
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
        #Draw a line from first date to next date
        if time_ord:
            
            
            #plt.arrow(targ_qmeans[0], targ_umeans[0],
            #          targ_qmeans[1]-targ_qmeans[0], targ_umeans[1]-targ_umeans[0],
            #          width=0.00005, head_width=0.0001,
            #          head_length=0.0001, color='black',
            #          alpha=0.2, linestyle='--')
                      
                #          targ_qmeans[i+1], targ_umeans[i+1], width=0.00005,
                #          head_width=0.0005, head_length=0.0000005, color='black')
            for i in range(0, len(targ_qmeans)-1):
            # (starting_x, starting_y, dx, dy, ...)
                
                #plt.plot( [targ_qmeans[i], targ_qmeans[i+1]], [targ_umeans[i], targ_umeans[i+1]], linestyle='--', alpha=0.2 ) 
                
                #obviously some scaling computations would be needed. Compute arrow size as some fraction of what?
                xlen = targ_qmeans[i+1]-targ_qmeans[i]
                ylen = targ_umeans[i+1]-targ_umeans[i]
                
                plt.arrow(targ_qmeans[i], targ_umeans[i],
                          xlen, ylen,
                          width=0.0001, head_width=0.0005,
                          length_includes_head=True,
                          head_length=(1/5)*np.sqrt(xlen**2 + ylen**2), color='black',
                          alpha=0.2, linestyle='--')
                
                #plt.arrow(targ_qmeans[i], targ_umeans[i], 
                #          targ_qmeans[i+1], targ_umeans[i+1] , width=0.00005,
                #          head_width=0.0005, head_length=0.0000005, color='black')
                #print( targ_qmeans[i], targ_umeans[i], "--------------->" ,  targ_qmeans[i+1], targ_umeans[i+1]      )
    else:
        plt.scatter(targ_qs, targ_us, 
                    color = plot_c, alpha=0.22)
        plt.errorbar(targ_qs, targ_us, 
                     xerr=targ_qstds, yerr=targ_ustds, 
                     lw=0.75, fmt="o", color=plot_c, alpha=0.22)
       
    if(reflection_axis=='pos'):
        plt.plot( [0, 0.01], [0, 0.01] , linestyle='--', alpha=0.7 ) #left side x y    
    elif(reflection_axis=='neg'):
        plt.plot( [0,-0.01], [0, -0.01] , linestyle='--', alpha=0.7 ) #left side x y
    elif(reflection_axis=='xy'):
        plt.plot( [-0.01,0.01], [-0.01,0.01] , linestyle='--', alpha=0.7 ) #left side x y        
           
    if(pol_deg):
        for z in range(0, len(targ_qmeans)):
            plt.plot([0,targ_qmeans[z]], [0, targ_umeans[z]], 
                     lw=1.75, alpha=0.4, linestyle = '--')
    if(verbose_MJD_arg):
        for z in range(0, len(targ_qmeans)):
            plt.text(targ_qmeans[z], 
                     targ_umeans[z], 
                     str(int(t[z].mjd)), rotation=45, fontsize=16)   #'int' object is not subscriptable    
    
    plt.grid()
    plt.title(in_title + " Q and U Scatter")        
        
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
        #plt.title(" Pol Scatter")
        plt.savefig(sv_im,bbox_inches='tight',pad_inches=0.1 , facecolor='w' )
    
    plt.show()
    
    #Returns what what you plot... It strips.
    if(retrun_plotelems):
        return(targ_qmeans, targ_qmeans_err, targ_umeans, targ_umeans_err)
    
def q_n_u_stack_plot_v2(pol_data, 
                        cal_data=[],
                        sv_im_str='',
                        pol_deg=False, 
                        key_verb=False):
    """
    A function that takes in data and plots q and u scatter and plots it in all separate sections. Does not return anything

    Parameters
    ----------
        input_data : tuple
            Data structure that contains all the data. Must inclue zero pol, high pol, and target
        sv_im_str  : str
            String to specify the filename of the image output
        pol_deg : bool, optional 
            Plots a line from 0 to the mean scatter. Not used? 
        key_verb : bool, optional
            Counts all data 
    """
    dtypes = ['zero pol', 'high pol', 'target']
    if(key_verb):
        for k in range(0,3):
            print( dtypes[k] ,  len(pol_data[k][0]))
    
    #####
    targ_date_strs = []
    zero_pol_date_strs = []
    high_pol_date_strs = []
    
    #####
    target_qs = []
    target_us = []
    
    #####
    z_pol_qs = []    
    z_pol_us = []  
    
    #####
    h_pol_qs = []
    h_pol_us = []
    
    #####targ_qmeans_err
    targ_qmeans = []
    targ_qmeans_err = []
    targ_umeans = []
    targ_umeans_err = []
    
    #####
    cal_targ_qmeans = []
    cal_targ_qmeans_err = []
    cal_targ_umeans = []
    cal_targ_umeans_err = []
    
    #####
    z_pol_qmeans = []
    z_pol_qmeans_err = []
    z_pol_umeans = []
    z_pol_umeans_err = []
    
    #####
    h_pol_qmeans = []
    h_pol_qmeans_err = []
    h_pol_umeans = []
    h_pol_umeans_err = []
    
    #####
    cal_h_pol_qmeans = []
    cal_h_pol_qmeans_err = []
    cal_h_pol_umeans = []
    cal_h_pol_umeans_err = []
    

    for i in range(0,3):
        if(i==0):
            for l in range(0, len(pol_data[i][0])):
                z_pol_qmeans.append(np.mean(pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][0]))
                z_pol_qmeans_err.append(np.mean(pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][1]))
                z_pol_umeans.append(np.mean(  pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][2]))
                z_pol_umeans_err.append(np.mean(  pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][3]))         
        elif(i==1):
            for l in range(0, len(pol_data[i][0])):
                h_pol_qmeans.append(np.mean(  pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][0]))
                h_pol_qmeans_err.append(np.mean(  pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][1]))
                h_pol_umeans.append(np.mean(  pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][2]))
                h_pol_umeans_err.append(np.mean(  pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][3]))
                
                if(len(cal_data) !=0):
                    cal_h_pol_qmeans.append(np.mean( cal_data[1][0][l][list(cal_data[1][0][l].keys())[0]][0] ))
                    cal_h_pol_qmeans_err.append(np.mean( cal_data[1][0][l][list(cal_data[1][0][l].keys())[0]][1] ))
                    cal_h_pol_umeans.append(np.mean( cal_data[1][0][l][list(cal_data[1][0][l].keys())[0]][2] ))
                    cal_h_pol_umeans_err.append(np.mean( cal_data[1][0][l][list(cal_data[1][0][l].keys())[0]][3] ))
        elif(i==2):
            for l in range(0, len(pol_data[i][0])):
                targ_qmeans.append(np.mean(  pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][0]))
                targ_qmeans_err.append(   np.mean(  pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][1]   )     )
                targ_umeans.append(np.mean(  pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][2]))
                targ_umeans_err.append(np.mean(  pol_data[i][0][l][list(pol_data[i][0][l].keys())[0]][3])       )
                
                if(len(cal_data) !=0):
                    cal_targ_qmeans.append(np.mean( cal_data[0][0][l][list(cal_data[0][0][l].keys())[0]][0] ))
                    cal_targ_qmeans_err.append(np.mean( cal_data[0][0][l][list(cal_data[0][0][l].keys())[0]][1] ))
                    cal_targ_umeans.append(np.mean( cal_data[0][0][l][list(cal_data[0][0][l].keys())[0]][2] ))
                    cal_targ_umeans_err.append(np.mean( cal_data[0][0][l][list(cal_data[0][0][l].keys())[0]][3] ))               
    
    simp_counter = 1
    
    #print(len(targ_qmeans),len(targ_umeans),len(targ_qmeans_err),len(targ_qmeans_err) )
    c_wheel = ['red','green','blue']
    cal_c_wheel = ['cyan', 'magenta', 'yellow'] 
    
    fig, axs = plt.subplots(2, 2)
    fig.suptitle("Pol Scatter Master Plot")
    alpha_c = 0.8
    black_lw = 0.8
    for r in range(0,2):
        for o in range(0,2):
            if(r==0 and o==0):            
                axs[r, o].scatter(targ_qmeans, targ_umeans, color = c_wheel[0], alpha=alpha_c)
                axs[r, o].errorbar(targ_qmeans, targ_umeans, xerr=targ_qmeans_err, yerr=targ_umeans_err, 
                                   lw=0.75, fmt="o", color=c_wheel[0], alpha=0.1)
                
                if(len(cal_data) !=0):
                    axs[r, o].scatter(cal_targ_qmeans, cal_targ_umeans, color = cal_c_wheel[0], alpha=alpha_c)
                    axs[r, o].errorbar(cal_targ_qmeans, cal_targ_umeans, xerr=cal_targ_qmeans_err, yerr=cal_targ_umeans_err, 
                                   lw=0.75, fmt="*", color=cal_c_wheel [0], alpha=0.1)
                    axs[r, o].legend(['Uncalibrated Target', 'Calibrated Target'])
                    
                axs[r, o].legend(['Uncalibrated Target'] )
                axs[r, o].axhline(y=0, color = 'black', lw=black_lw)
                axs[r, o].axvline(x=0, color = 'black', lw=black_lw)
                axs[r, o].grid()
                axs[r, o].set_title('Target')
                
            elif(r==1 and o==0):
                axs[r, o].scatter(z_pol_qmeans, z_pol_umeans, color = c_wheel[1], alpha=alpha_c)
                axs[r, o].errorbar(z_pol_qmeans, z_pol_umeans, xerr=z_pol_qmeans_err, yerr=z_pol_umeans_err, 
                                   lw=0.75, fmt="o", color=c_wheel[1], alpha=0.1)
                
                axs[r, o].axhline(y=0, color = 'black', lw=black_lw)
                axs[r, o].axvline(x=0, color = 'black', lw=black_lw)
                axs[r, o].grid()
                axs[r, o].set_title('Zero Polarization')
                
            elif(r==0 and o==1):
                axs[r, o].scatter(h_pol_qmeans, h_pol_umeans, color = c_wheel[2], alpha=alpha_c )
                axs[r, o].errorbar(h_pol_qmeans, h_pol_umeans, xerr=h_pol_qmeans_err, yerr=h_pol_umeans_err, 
                                   lw=0.75, fmt="o", color=c_wheel[2], alpha=0.1)
                
                if(len(cal_data) !=0):
                    axs[r, o].scatter(cal_h_pol_qmeans, cal_h_pol_umeans, color = cal_c_wheel[2], alpha=alpha_c)
                    axs[r, o].errorbar(cal_h_pol_qmeans, cal_h_pol_umeans, xerr=cal_h_pol_qmeans_err, yerr=cal_h_pol_umeans_err, 
                                   lw=0.75, fmt="*", color=cal_c_wheel [2], alpha=0.1)
                    axs[r, o].legend(['Uncalibrated High Pol', 'Calibrated  High Pol'] )
                    
                axs[r, o].legend(['Uncalibrated High Pol'])
                axs[r, o].axhline(y=0, color = 'black', lw=black_lw)
                axs[r, o].axvline(x=0, color = 'black', lw=black_lw)
                axs[r, o].grid()
                axs[r, o].set_title('High Polarization')
            elif(r==1 and o==1):
                axs[r, o].scatter(targ_qmeans, targ_umeans, color = c_wheel[0], alpha=alpha_c )
                axs[r, o].errorbar(targ_qmeans, targ_umeans, xerr=targ_qmeans_err, yerr=targ_umeans_err, 
                                   lw=0.75, fmt="o", color=c_wheel[0], alpha=0.1)
                
                axs[r, o].scatter(z_pol_qmeans, z_pol_umeans, color = c_wheel[1], alpha=alpha_c )
                axs[r, o].errorbar(z_pol_qmeans, z_pol_umeans, xerr=z_pol_qmeans_err, yerr=z_pol_umeans_err, 
                                   lw=0.75, fmt="o", color=c_wheel[1], alpha=0.1)
                
                axs[r, o].scatter(h_pol_qmeans, h_pol_umeans, color = c_wheel[2], alpha=alpha_c )
                axs[r, o].errorbar(h_pol_qmeans, h_pol_umeans, xerr=h_pol_qmeans_err, yerr=h_pol_umeans_err, 
                                   lw=0.75, fmt="o", color=c_wheel[2], alpha=0.1)
                
                if(len(cal_data) !=0):
                    axs[r, o].scatter(cal_targ_qmeans, cal_targ_umeans, color = cal_c_wheel[0], alpha=alpha_c)
                    axs[r, o].errorbar(cal_targ_qmeans, cal_targ_umeans, xerr=cal_targ_qmeans_err, yerr=cal_targ_umeans_err, 
                                   lw=0.75, fmt="*", color=cal_c_wheel [0], alpha=0.1)
                    
                    axs[r, o].scatter(cal_h_pol_qmeans, cal_h_pol_umeans, color = cal_c_wheel[2], alpha=alpha_c)
                    axs[r, o].errorbar(cal_h_pol_qmeans, cal_h_pol_umeans, xerr=cal_h_pol_qmeans_err, yerr=cal_h_pol_umeans_err, 
                                   lw=0.75, fmt="*", color=cal_c_wheel [2], alpha=0.1)
                   
                    axs[r, o].legend(['Uncalibrated Target', 
                                      'Zero pol', 
                                      'Uncalibrated High Pol',
                                      'Calibrated High Pol',
                                      'Calibrated Target'])     
                    
                axs[r, o].legend(['Uncalibrated Target', 'Zero pol', 'Uncalibrated High Pol'])
                axs[r, o].axhline(y=0, color = 'black', lw=black_lw)
                axs[r, o].axvline(x=0, color = 'black', lw=black_lw)

                axs[r, o].grid()
                axs[r, o].set_title('Data Combined')  
    if(sv_im_str != ''):
        plt.savefig(sv_im_str, bbox_inches='tight', pad_inches=0.1 , facecolor='w')
        plt.show()   

def mean_q_u_check(inp_data, 
                   n, 
                   q_u_ret, 
                   verb_arg=False):
    """
    A function that checks the computation of the mean 

    Parameters
    ----------
        input_data : tuple
            Data structure that contains all the data. Must inclue zero pol, high pol, and target
        n  : int
            integer
        q_u_ret : str
            Plots a line from 0 to the mean scatter 
        verb_arg : bool, optional
            Counts all data 
    """
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