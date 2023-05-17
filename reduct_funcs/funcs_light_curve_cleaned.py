import astropy
import gzip
import glob
import numpy as np
import math
import xlrd
import datetime
import matplotlib.pyplot as plt
import pickle

from astropy.time import Time,  TimeJD
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
    """
    #
    """
    filename = 'EECEPLightCurve.txt'
    #print("Reading:", filename )
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

    color_arr = [B_arr, V_arr[57:], R_arr, I_arr] #This have already been taken out is it?
    mjd_arr_k = [t_b, t_v[57:], t_r, t_i]
            
    return(color_arr, mjd_arr_k) 

def compute_phase(PD_data,
              PA_data,
              light_curve_data,
              relative_phase=False,
              verbose_text=False):
    """
    #Inputs are double polarization data
    """
    if len(light_curve_data) == 0:
        #print("Loading Default light curve")
        def_light_curve = load_default_light_curve()
        
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
    pcolor = ['B', 'V', 'R', 'I']
    
    if(verbose_text):
        print("Time to phase translation:")
        for k in range(0, len(mjd_arr_k[:-1])):#things in mjd_arr_k[:-1]:
            print("Length", pcolor[k], ":" , len(mjd_arr_k[k]))    

    phase_axis = []
    x_axis = []
    
    tz = Time(2456894, format='jd') #this is your zero
    mjdz = tz.mjd
        
    # Calculate the number of seconds in 5.6 years
    seconds_in_year = 31536000
    seconds_in_5_6_years = int(5.6 * seconds_in_year)

    # Create a timedelta object representing 5.6 years
    delta_years = datetime.timedelta(seconds=seconds_in_5_6_years)
    
    tz2020ecplipse = tz.to_datetime() + delta_years #2020
    tz2025ecplipse = tz2020ecplipse + delta_years   #2025
    twentynineclipse = tz .to_datetime() - delta_years #2009
    twentyofoureclipse = twentynineclipse - delta_years #2014
    
    #The computation is taking ages. Instead of appending to a new data structure
    #Just overwrite the values.
    for k in mjd_arr_k[4]: #This will be a problem if there are other eclipses.
        #TODO: Our data does not span beyond 2020... well technically there is that little scrap of the 2014 eclipse that exists
        #This needs to be made more robust. To include more
        if(k.mjd<Time(tz2020ecplipse).mjd and k.mjd >  mjdz):
            if(verbose_text):
                print("Left!",k, k.mjd, (k.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz) )
            if(relative_phase):
                x_axis.append((k.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz) -1  )
            else:
                x_axis.append((k.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz))   #This is correct. It is left
        elif(k.mjd>Time(tz2020ecplipse).mjd ):
            if(verbose_text):
                print("Right!",k, k.mjd, (k.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz)-1 )
            if(relative_phase):
                x_axis.append( (k.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz)-1)
            else:
                x_axis.append( (k.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz)-1)
            
        
    #This is the light curve
    for timelists in mjd_arr_k[:-1]: #this is a recursive call. Just recursive call the third plot   
        lctrans = []
        for things in timelists:
            if(things.mjd < Time(twentynineclipse).mjd): #phase --3 - 2
                if(verbose_text):
                    print(things, 1+(things.mjd-Time(twentynineclipse).mjd)/(Time(tz2020ecplipse).mjd - mjdz), 'before 2009' )
                if(relative_phase):
                    lctrans.append(  (1+ (things.mjd-Time(twentynineclipse).mjd)/(Time(tz2020ecplipse).mjd - mjdz)) - 3  )                
                else:
                    lctrans.append(  1+ (things.mjd-Time(twentynineclipse).mjd)/(Time(tz2020ecplipse).mjd - mjdz)  )
        
            elif(things.mjd < mjdz and things.mjd > Time(twentynineclipse).mjd): #phase -2 - 1
                if(verbose_text):
                    print(things,  (things.mjd-Time(twentynineclipse).mjd)/(Time(tz2020ecplipse).mjd - mjdz), 'after 2009 before 2014' )
                if(relative_phase):
                    lctrans.append(   ((things.mjd-Time(twentynineclipse).mjd)/(Time(tz2020ecplipse).mjd - mjdz)) -2  )                    
                else:
                    lctrans.append((things.mjd-Time(twentynineclipse).mjd)/(Time(tz2020ecplipse).mjd - mjdz))
                    
            elif(things.mjd<Time(tz2020ecplipse).mjd and things.mjd >  mjdz): #phase -1 - 0
                if(verbose_text):
                    print(things, (things.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz), 'after 2014 before 2020'  )
                if(relative_phase):
                    lctrans.append( ( (things.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz) ) -1                       )              
                else:
                    lctrans.append((things.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz)) 
                    
            elif(things.mjd>Time(tz2020ecplipse).mjd and things.mjd < Time(tz2025ecplipse).mjd): #phase 0 - 1
                if(verbose_text):
                    print(things,  (things.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz)-1, 'after 2020 before 2025') 
                if(relative_phase):
                    lctrans.append((things.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz)-1)
                else:
                    lctrans.append((things.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz)-1)
            elif(things.mjd>Time(tz2025ecplipse).mjd ): #phase 1 - 2
                if(verbose_text):
                    print(things,  (things.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz)-1, 'after 2025')
                if(relative_phase):
                    lctrans.append(    ((things.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz)-1)+1    )
                else:
                    lctrans.append((things.mjd-mjdz)/(Time(tz2020ecplipse).mjd - mjdz)-1)
            else:
                if(verbose_text):
                    print("No Date computed. Error here:", things.mjd) 
                                
        phase_axis.append(lctrans)
        
    return(phase_axis, x_axis)# , color_arr)

def EECep_light_curve_stacked(PD_data, #PD, PDerror, MJD
                              PA_data,
                              light_curve_data,
                              MJD_cuts,
                              phase_plot=False,
                              verb_t=False,
                              true_pa=False, 
                              title_t=False):
    """
    A function that takes in polarization master plot.
    
    Parameters
    ----------
    PD_data : 
        Polarization Degree data
    PA_data : 
        Polarization Angle data
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
        
    if(phase_plot):
        print("Phase Shift Translation:")
        
        phase_shift_data = process_A(PD_data, #PD, PDerror, MJD
                              PA_data,
                              [],
                              verbose_text=False,
                              relative_phase=True)
        
        #print("Expectation is 3",len(phase_shift_data))
        
        #phase_shift_axis=phase_shift[0]
        """
        phase_axis = []
        t = Time(2456894, format='jd') #this is your zero
        mjdz = t.mjd
        
        # Calculate the number of seconds in 5.6 years
        seconds_in_year = 31536000
        seconds_in_5_6_years = int(5.6 * seconds_in_year)

        # Create a timedelta object representing 5.6 years
        delta_years = datetime.timedelta(seconds=seconds_in_5_6_years)
        twtyeclipse = t.to_datetime() + delta_years
        twentynineclipse = t.to_datetime() - delta_years
        
        print("Spans 0 - 2009:", mjdz - Time(twentynineclipse).mjd )
        print("Spans 2020 - 2014:", Time(twtyeclipse).mjd - mjdz )
        print(twtyeclipse, twentynineclipse)
        
        x_axis = []
        for k in mjd_arr_k[4]:
            if(k.mjd<Time(twtyeclipse).mjd and k.mjd >  mjdz):
                x_axis.append( (k.mjd-mjdz)/(Time(twtyeclipse).mjd - mjdz))
                
            elif(k.mjd>Time(twtyeclipse).mjd ):
                x_axis.append( (k.mjd-mjdz)/(Time(twtyeclipse).mjd - mjdz)-1)
        
        #This is the light curve
        for timelists in mjd_arr_k[:-1]: #in range(0, len(mjd_arr_k)-1): #this is a recursive call. Just recursive call the third plot   
            lctrans = []
            for things in timelists:
                #print(things)
                if(things.mjd< mjdz and things.mjd > Time(twentynineclipse).mjd ):
                    print(things,  (things.mjd-Time(twentynineclipse).mjd)/(Time(twtyeclipse).mjd - mjdz), 'before 2014'    )
                    lctrans.append((things.mjd-Time(twentynineclipse).mjd)/(Time(twtyeclipse).mjd - mjdz))
                    
                elif(things.mjd<Time(twtyeclipse).mjd and things.mjd >  mjdz):
                    print(things, (things.mjd-mjdz)/(Time(twtyeclipse).mjd - mjdz), 'after 2014 before 2020'  )
                    lctrans.append((things.mjd-mjdz)/(Time(twtyeclipse).mjd - mjdz))
                    
                elif(things.mjd>Time(twtyeclipse).mjd ):
                    print(things,  (things.mjd-mjdz)/(Time(twtyeclipse).mjd - mjdz)-1, 'after 2020')
                    lctrans.append((things.mjd-mjdz)/(Time(twtyeclipse).mjd - mjdz)-1)
            phase_axis.append(lctrans)
        """
        
    #print("xaxis:",  x_axis) 
    #x=range(len(x_axis))
    #return_val[0] #BVRI phase shifts
    #return_val[1] #PD data phase shifts
    #return_val[2]#BVRI magnitudes
    #print("This is your new x axis:",phase_shift_data[0])
    #print("This should be the same size:",  
    #      len(mjd_arr_k[4][index_L_cutt[4]:index_R_cutt[4]+1].mjd),
    #      len(phase_shift_data[1])  )
    #print( mjd_arr_k[4][index_L_cutt[4]:index_R_cutt[4]+1].mjd )
    #print(phase_shift_data[1]  )
    #print(len(phase_shift_data[0][0]), len(color_arr[0]) ) #should be the same.
    #print(len(phase_shift_data[0][1]), len(color_arr[1]) ) 
    #print(len(phase_shift_data[0][2]), len(color_arr[2]) ) 
    #print(len(phase_shift_data[0][3]), len(color_arr[3]) ) 
    
    #You are still calculating the old x axis. 
        
    fig = plt.figure(figsize=(36, 12))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True, sharey=False)    

    if(phase_plot):
        markers, caps, bars = axs[0].errorbar(
                phase_shift_data[1] , 
                PD_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
                yerr=PD_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
                xerr =[0]*len(PD_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
                color='black', fmt='*', markersize=16, capsize=10, capthick=5, label='PD')
        
        #You had attempted it earlier. Bravo
        #markers, caps, bars = axs[0].errorbar(
        #    x, 
        #    PD_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
        #    yerr=PD_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
        #    xerr =[0]*len(PD_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
        #    color='black', fmt='*', markersize=16, capsize=10, capthick=5, label='PD')
        
    else:
        markers, caps, bars = axs[0].errorbar(
                mjd_arr_k[4][index_L_cutt[4]:index_R_cutt[4]+1].mjd, 
                PD_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
                yerr=PD_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
                xerr =[0]*len(PD_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
                color='black', fmt='*', markersize=16, capsize=10, capthick=5, label='PD')
    
    [bar.set_alpha(0.15) for bar in bars]
    [cap.set_alpha(0.645) for cap in caps]
    
    if(phase_plot):
        markers, caps, bars = axs[1].errorbar(
                phase_shift_data[1] , 
                PA_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
                yerr=PA_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
                xerr =[0]*len(PA_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
                color='black', fmt='*', markersize=16, capsize=10, capthick=5, label='PA')
        #You had attempted it earlier. Bravo
        #markers, caps, bars = axs[1].errorbar(
        #        x, 
        #        PA_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
        #        yerr=PA_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
        #        xerr =[0]*len(PA_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
        #        color='black', fmt='*', markersize=16, capsize=10, capthick=5, label='PA')
        
    else:
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
    
    if(phase_plot):
        for ka in range(0, len(mjd_arr_k)-1): #this is a recursive call. Just recursive call the third plot
            axs[2].scatter(phase_shift_data[0][ka][index_L_cutt[ka]:index_R_cutt[ka]], 
                        color_arr[ka][index_L_cutt[ka]:index_R_cutt[ka]], 
                        alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])    
    
    #This is the light curve
    #    for ka in range(0, len(phase_axis)): #this is a recursive call. Just recursive call the third plot
    #        axs[2].scatter(phase_axis[ka][:], 
    #                    color_arr[ka][:], 
    #                    alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])
    else:   
        for ka in range(0, len(mjd_arr_k)-1): #this is a recursive call. Just recursive call the third plot
            axs[2].scatter(mjd_arr_k[ka][index_L_cutt[ka]:index_R_cutt[ka]].mjd, 
                        color_arr[ka][index_L_cutt[ka]:index_R_cutt[ka]], 
                        alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])
    
    if(true_pa):
        fig.suptitle('Polarization degree (PD), true Position angle (PA) and EE CEP light curve', 
                         fontsize=32)    
    else:
        fig.suptitle('Polarization degree (PD), Position angle (PA) and EE CEP light curve', fontsize=32)
        
    if(phase_plot): 
        k=0
        #axs[2].axvline(Time(twtyeclipse).mjd) #this is definitely one year later
        #axs[2].axvline(x=mjdz )
        #axs[2].set_xticks(x)
    else:
        k=0
               
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
                     phase_plot = False,
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
    
    if(phase_plot):
        print("Phase Shift:")
        
        phase_shift_data = process_A(pol_data, #PD, PDerror, MJD
                              pol_data,
                              [],
                              verbose_text=False,
                              relative_phase=True)
    
    #fig = plt.figure(figsize=(36, 12))    
    fig, ax1 = plt.subplots(figsize=(36, 12))
    
    #Plot colour magnitudes
    if(phase_plot):    
        #I do not need to return the magnitudes.
        for ka in range(0, len(mjd_arr_k)-1): #this is a recursive call. Just recursive call the third plot
            ax1.scatter(phase_shift_data[0][ka][index_L_cutt[ka]:index_R_cutt[ka]], 
                        color_arr[ka][index_L_cutt[ka]:index_R_cutt[ka]], 
                        alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])
    else:
        for ka in range(0, len(mjd_arr_k)-1): #this is a recursive call. Just recursive call the third plot
            ax1.scatter(mjd_arr_k[ka][index_L_cutt[ka]:index_R_cutt[ka]].mjd, 
                        color_arr[ka][index_L_cutt[ka]:index_R_cutt[ka]], 
                        alpha = 0.72, s=16, color=c_arr[ka], label=c_labs[ka])

    fig.gca().invert_yaxis()
    #"""
    ax2 = ax1.twinx()    
    
    #Plot pol data
    if(phase_plot):
        markers, caps, bars = ax2.errorbar(
                phase_shift_data[1], 
                pol_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
                yerr=pol_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
                xerr =[0]*len(pol_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
                color='black', fmt='*', markersize=16, capsize=10, capthick=5, label=txt_arg)        
    else:
        markers, caps, bars = ax2.errorbar(
                mjd_arr_k[4][index_L_cutt[4]:index_R_cutt[4]+1].mjd, 
                pol_data[0][index_L_cutt[4]:index_R_cutt[4]+1], 
                yerr=pol_data[1][index_L_cutt[4]:index_R_cutt[4]+1], 
                xerr =[0]*len(pol_data[0][index_L_cutt[4]:index_R_cutt[4]+1]), 
                color='black', fmt='*', markersize=16, capsize=10, capthick=5, label=txt_arg)
    
    #plt.set_ylabel('PD, (%)', fontsize=32)
    #plt.tick_params(axis="y", labelsize=28)
    #"""
    [bar.set_alpha(0.2) for bar in bars]
    [cap.set_alpha(0.72) for cap in caps]  

    if(txt_arg == 'PD'):
        ax2.set_title("EE Cep light curve and Polarization Degree ("+txt_arg+")", fontsize=32)
        ax2.set_ylabel('PD, (%)', fontsize=32)
        ax2.tick_params(axis="y", labelsize=28)
        
    elif(txt_arg == 'PA'):
        if(true_pa):
            ax1.set_title("EE Cep light curve and true Position Angle ("+txt_arg+")", fontsize=32)

        else:
            ax1.set_title("EE Cep light curve and Position Angle ("+txt_arg+")", fontsize=32)    
        ax2.set_ylabel('PA, ('+u'\N{DEGREE SIGN}'+')', fontsize=32)
    else:
        print("Insufficient polarization input")

    """
    """

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
    