import gzip
import shutil
import glob
import astropy
import matplotlib.pyplot as plt
import numpy as np
import astropy
import math
import xlrd
import re

from astropy.time import Time
from datetime import datetime
from copy import deepcopy as cp
from scipy import stats
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.aperture import EllipticalAperture, EllipticalAnnulus
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

def source_peak_finder_pol_std(fits_data_1, siegma, search_array, trim, plot_peaks, verbose):
    """
    Function that finds coordinates of stars.
    Star Finder used for the really faint stars of the pol std group. Haha you were so off.
    So I guess this faint source detection thing is kind of pointless now. All of the targets are observeable
    
    #How did you parse that sum bitch?
    """
    #Some things that hold data
    x_peak = []
    y_peak = []
    peak_val =[]

    x_interest = []
    y_interest = []
    peak_interest = []
    index_x = []

    #First layer thresholding
    mean, median, std = sigma_clipped_stats(fits_data_1[0].data, sigma=siegma) #searches the whole image
    first_threshold = median + (10. * std)
    tbl = find_peaks(fits_data_1[0].data, first_threshold, box_size=40)
    tbl['peak_value'].info.format = '%.3g'  #format threshold
    
    #stripper
    for things in tbl:
        x_peak.append(things['x_peak'])
        y_peak.append(things['y_peak'])
        peak_val.append(things['peak_value'])

    for i in range(0, len(peak_val)): #looks inside zone
        """
        #TODO replace with a circle
        #Calculate a distance r
        #if distnace to point (d = sqrt( (x-x)^2 + (y-y)^2) is less than some input d#search offset is the center
        """
        if(y_peak[i] > 512-search_array[0] and 
           y_peak[i] < 512+search_array[1] and 
           x_peak[i] > 512-search_array[2] and 
           x_peak[i] < 512+search_array[3]):
            x_interest.append(x_peak[i])
            y_interest.append(y_peak[i])#theres this really bright pixel off to the side
            peak_interest.append(peak_val[i])
            index_x.append(i)
    
    x_targ = []
    y_targ = []
    peak_targ = []
    sub_sample = []
    
    for peaks in peak_val:
        if peaks < trim: #we only append faint sources
            sub_sample.append(peaks) #But it cannot identify it out of this sample...
            
    second_threshold = np.mean(sub_sample)+5*np.std(sub_sample) #This thing becomes nan and shit is fucked
    if math.isnan(second_threshold):
        second_threshold = 50 #if it is nan then we set to some threshold. Sometimes that shit is so faint that not even the 
                               #manual threshold can get a fix on it
        
    #second_threshold = some hard coded value. I just want coordinates!
    
    for z in range(0, len(peak_interest)):
        if(peak_interest[z] > second_threshold): #then this becomes nan and shit is truly fucked
            x_targ.append(x_interest[z])
            y_targ.append(y_interest[z])
            peak_targ.append(peak_interest[z])
    #Output
    if(plot_peaks):
        plt.plot(peak_val[:])
        plt.title("Flux peaks")
        plt.grid()
        
        plt.text(0, second_threshold,
                 'Region of interest threshold ', 
                 fontsize=16, alpha=10, color = 'red')
        plt.axhline(y=second_threshold, 
                    color = 'red', linestyle='--', lw=2.5) #horizontal line y is constan
        plt.text(0, trim,
                         'Trim', 
                         fontsize=16, alpha=10, color = 'red')
        plt.axhline(y=trim, 
                            color = 'red', linestyle='--', lw=2.5) #horizontal line y is constan
        for i in range(0, len(index_x)):
            plt.plot([index_x[i],index_x[i]], [peak_interest[i]-300, peak_interest[i]+300], 'k-', lw=2.5) #vertical line x is constant
            plt.plot([index_x[i]-2.5,index_x[i]+2.5], [peak_interest[i], peak_interest[i]] , 'k-', lw=2.5) #horizontal line y is constant
        plt.show()

    if(verbose):
        print(len(tbl), 
              "peaks detected from image of size", 
              fits_data_1[0].header['NAXIS1'], "x", fits_data_1[0].header['NAXIS2'],
              "with sigma:", siegma, "and second threshold:", second_threshold )

        print("Targets within region of interest: ", len(x_targ))

    return(x_targ, y_targ, peak_targ)
    
def check_pol_std():
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
    
#def load_pol_data(path, f_name, verbose):
def load_pol_data(in_str, verbose):
    """
    #This should return a dictioanary with the key being a date.
    #
    """
    
    #Parse it top layer first...
    path_s = re.search('(.*)master_', in_str) #squared away
    f_name_s = re.search('master_(.*)', in_str)
    
    f_name = 'master_' + f_name_s.group(1)
    
    result_MJD = re.search('master_(.*)_', f_name[:19])
    result_obj = re.search('_(.*)_P', f_name[15:])
    if(verbose):
        print("File path:", path)
        print("File name:", f_name)
        print("MJD:", result_MJD.group(1))
        print("obj:", result_obj.group(1))
    
    file_dir = path_s.group(1)+f_name
    loc = (file_dir)

    wb = xlrd.open_workbook(loc)
    sheet = wb.sheet_by_index(0)

    q_arr = []
    q_err_arr = []
    u_arr = []
    u_err_arr = []

    for j in range(0, sheet.nrows ):
        q_arr.append(sheet.cell_value(j, 26))
        q_err_arr.append(sheet.cell_value(j, 27))
        u_arr.append(sheet.cell_value(j, 28))
        u_err_arr.append(sheet.cell_value(j, 29))
        
    ret_dict = {result_MJD.group(1)+"_"+ result_obj.group(1): (q_arr,  q_err_arr, u_arr, u_err_arr)} #inside the plotting script
    #pare the left side of the _ as MJD and the right side as the target.
    #the you just plot for each
    return ret_dict
    
def q_u_stats(q, q_err, u, u_err):
    mean_q = np.mean(q)
    mean_u = np.mean(u)
    
    median_q = np.median(q)
    median_u = np.median(u)
    
    std_q = np.std(q)
    std_u = np.std(u)
    
    return(mean_q, mean_u, median_q, median_u, std_q, std_u)

def correct_q_u(target_data, zero_pol_std, high_pol_std, zero_pol_offset):   
    target_data[0][1:] = target_data[0][1:] - zero_pol_offset[0]
    target_data[2][1:] = target_data[2][1:] - zero_pol_offset[1]
    
    zero_pol_std[0][1:] = zero_pol_std[0][1:] - zero_pol_offset[0]
    zero_pol_std[2][1:] = zero_pol_std[2][1:] - zero_pol_offset[1]
    
    high_pol_std[0][1:] = high_pol_std[0][1:] - zero_pol_offset[0]
    high_pol_std[2][1:] = high_pol_std[2][1:] - zero_pol_offset[1]
    
    return(target_data, zero_pol_std, high_pol_std)

def calc_pd2(input_data, plot_title, plot_c,sv_im_str ,perc_arg ,calc_pd_verbose, sv_arg):
    """
    #Input 
    
    #Does not plot as a function of MJD
    """
    print("Calculate Polarization Degree:")
    pol_d_array = []
    pol_d_err_array = []
    array_des_dates = []
    array_des_nom = []
    mjd_strs = []
    
    for dats in input_data:        
        result = re.search('(.*)_', list(dats.keys())[0][:12])
        targ_name_str = re.search('_(.*)', list(dats.keys())[0])
        
        array_des_nom.append(targ_name_str.group(1))
        array_des_dates.append(result.group(1))
        
        mjd_strs.append(result.group(1)+'T00:00:00.000000000')
        
        mean_q_squared = np.mean(dats[list(dats.keys())[0]][0][1:])**2
        mean_u_squared = np.mean(dats[list(dats.keys())[0]][2][1:])**2
        
        mean_q_err_squared = np.mean(dats[list(dats.keys())[0]][1][1:])**2
        mean_u_err_squared = np.mean(dats[list(dats.keys())[0]][3][1:])**2
        
        sum_o_squares = mean_q_squared + mean_u_squared
        
        pol_d = math.sqrt(sum_o_squares)
        #From Slowikowska et al 2016
        pol_d_err  = math.sqrt( ((mean_q_squared)/(sum_o_squares))*(mean_q_err_squared)+ ((mean_u_squared)/(sum_o_squares))*(mean_u_err_squared))
        
        if(perc_arg):
            if(calc_pd_verbose):
                print(targ_name_str.group(1), "MJD:", result.group(1), pol_d*100, u"\u00B1",  pol_d_err*100)
            pol_d_array.append(pol_d*100)
            pol_d_err_array.append(pol_d_err*100)
        else:
            if(calc_pd_verbose):
                print("MJD:", result.group(1), pol_d, u"\u00B1",  pol_d_err)#, "type date:", type(result.group(1) 
            pol_d_array.append(pol_d)
            pol_d_err_array.append(pol_d_err)

    t = Time(mjd_strs, format='isot', scale='utc')
    
    #If i'm plotting pd then  when it moves to another one it will change shapes
    plt.scatter(np.linspace(0, len(pol_d_array)-1, num=len(pol_d_array)), pol_d_array, marker ="^")
    plt.errorbar(np.linspace(0, len(pol_d_array)-1, num=len(pol_d_array)), pol_d_array, xerr=[0]*len(pol_d_array), yerr=pol_d_err_array, lw=0.75, fmt="^", color=plot_c, alpha=0.9)
    plt.title(plot_title)
    
    for h in range(0, len(t)):
        plt.text(h, pol_d_array[h], array_des_nom[h]+"\nMJD: "+str(int(t[h].mjd))+"\nPD:"+str(round(pol_d_array[h], 2))+u"\u00B1"+str(round(pol_d_err_array[h], 1)), fontsize=14 ,rotation = 45)

    plt.grid()
    plt.savefig(sv_im_str, bbox_inches='tight',pad_inches=0.1)
    plt.show()
    
def calc_pa2(input_data, plot_title, plot_c,sv_im_str  ,deg_arg , calc_pa_verbose, sv_arg):
    """
    #Input 
    """
    print("Calculate Polarization Angle:")
    pol_pa_array = []
    pol_pa_err_array = []
    array_des_dates = []
    array_des_noms = []
    
    #targ_date_strs = [] #array des dats is this
    mjd_strs = []
    for dats in input_data:        
        result = re.search('(.*)_', list(dats.keys())[0][:12])
        targ_name_str = re.search('_(.*)', list(dats.keys())[0] )
        
        array_des_noms.append(targ_name_str.group(1))
        mjd_strs.append(result.group(1)+'T00:00:00.000000000')
        array_des_dates.append(result.group(1)) #This is 
        
        mean_q_err_squared = np.mean(dats[list(dats.keys())[0]][1][1:])**2
        mean_u_err_squared = np.mean(dats[list(dats.keys())[0]][3][1:])**2
        
        mean_q_squared = np.mean(dats[list(dats.keys())[0]][0][1:])**2
        mean_u_squared = np.mean(dats[list(dats.keys())[0]][2][1:])**2
        
        mean_q = np.mean(dats[list(dats.keys())[0]][0][1:])
        mean_u = np.mean(dats[list(dats.keys())[0]][2][1:])
        
        sum_o_squares = mean_q_squared + mean_u_squared
        
        #So what do we need?
        #For RINGO Slowikoska et al
        pol_pa = 0.5*math.atan2(mean_u , mean_q )
        pol_pa_err = math.sqrt(((1/(2*mean_q*(1 + (mean_u_squared/mean_q_squared))))**2 )*(mean_u_err_squared) + ((-1*( (mean_u)/(2*sum_o_squares)))**2 )*(mean_q_err_squared))
                
        if(deg_arg):
            if(calc_pa_verbose): #
                print("MJD:", result.group(1), pol_pa*(180/3.142), u"\u00B1",  pol_pa_err*(180/3.142))
            pol_pa_array.append(pol_pa*(180/3.142)) 
            pol_pa_err_array.append(pol_pa_err*(180/3.142)) 
        else:
            if(calc_pa_verbose):
                print("MJD:", result.group(1), pol_pa, u"\u00B1",  pol_pa_err)#, "type date:", type(result.group(1) 
            pol_pa_array.append(pol_pa)
            pol_pa_err_array.append(pol_pa_err)

    t = Time(mjd_strs, format='isot', scale='utc')
    
    #This should be universal.
    #Plot_MJD
    #The input for PD and PA should be the same as q and u stability
    plt.scatter( np.linspace(0, len(pol_pa_array)-1, num=len(pol_pa_array)), pol_pa_array, marker ="^")
    plt.errorbar( np.linspace(0, len(pol_pa_array)-1, num=len(pol_pa_array)), pol_pa_array, xerr=[0]*len(pol_pa_array), yerr=pol_pa_err_array, lw=0.75, fmt="^", color=plot_c, alpha=0.9)
    plt.title(plot_title)
    plt.grid()
    
    for h in range(0, len( t )):
        plt.text(h, pol_pa_array[h] ,  array_des_noms[h]+"\nMJD: "+str(int(t[h].mjd))+"\nPA:"+str( round(pol_pa_array[h], 2))+u"\u00B1"+str(round(pol_pa_err_array[h], 1)), fontsize=14 ,rotation = 45)
    plt.savefig(sv_im_str, bbox_inches='tight',pad_inches=0.1)
    plt.show()  

    
    

    
#These are rather deprecated
#
#
#
#
#
#
def calc_pd(target_data, zero_pol_std, high_pol_std, name_array):
    """
    #Just do all?
    """
    print("Target Instances:", len(target_data[0]))
    print("Zero Pol Instances:", len(zero_pol_std[0]))
    print("High Pol Instances:", len( high_pol_std[0]), "\n")
    
    target_pd = math.sqrt( np.mean(target_data[0][1:])**2 + np.mean(target_data[2][1:])**2) #This works using scatter bro
    zero_pol_pd = math.sqrt(np.mean(zero_pol_std[0][1:])**2 + np.mean(zero_pol_std[2][1:])**2)
    high_pol_pd = math.sqrt(np.mean(high_pol_std[0][1:])**2 + np.mean(high_pol_std[2][1:])**2)
    
    targ_mean_q, targ_mean_u, median_q, median_u, targ_q_std, targ_u_std = q_u_stats(target_data[0][1:] ,target_data[1][1:] ,target_data[2][1:] ,target_data[3][1:])
    zero_pol_mean_q, zero_pol_mean_u, median_q, median_u, zero_pol_q_std, zero_pol_u_std = q_u_stats(zero_pol_std[0][1:] ,zero_pol_std[1][1:] ,zero_pol_std[2][1:] ,zero_pol_std[3][1:])
    high_pol_mean_q, high_pol_mean_u, median_q, median_u, high_pol_q_std, high_pol_u_std = q_u_stats(high_pol_std[0][1:] ,high_pol_std[1][1:] ,high_pol_std[2][1:] ,high_pol_std[3][1:])
    
    target_PD_err = math.sqrt( targ_q_std**2 +    targ_u_std**2           )
    zero_pol_PD_err = math.sqrt( zero_pol_q_std**2 +    zero_pol_u_std**2           )
    high_pol_PD_err = math.sqrt(  high_pol_q_std**2 +     high_pol_u_std **2           )
            
    print(name_array[0], "Target PD:", target_pd, u"\u00B1", target_PD_err, "% Pol:", target_pd*100, u"\u00B1", target_PD_err*100)
    print(name_array[1], "Zero pol PD:", zero_pol_pd, u"\u00B1", zero_pol_PD_err, "% Pol:",zero_pol_pd*100, u"\u00B1", zero_pol_PD_err*100)
    print(name_array[2], "High pol PD:", high_pol_pd, u"\u00B1", high_pol_PD_err, "% Pol:", high_pol_pd*100, u"\u00B1", high_pol_PD_err*100)
    
    return (target_pd, zero_pol_pd, high_pol_pd )

def calc_pa(target_data, zero_pol_std, high_pol_std, name_array):
    
    target_pa = 0.5*math.atan2(np.mean(target_data[2][1:]),np.mean(target_data[0][1:])) #2 then 0. u then q, y then x
    zero_pol_pa = 0.5*math.atan2(np.mean(zero_pol_std[2][1:]),np.mean(zero_pol_std[0][1:]))
    high_pol_pa = 0.5*math.atan2(np.mean(high_pol_std[2][1:]),np.mean(high_pol_std[0][1:]))
         
    print(name_array[0], "Target PA (rad):", target_pa, "Target PA (%):" ,target_pa*(180/3.142) )
    print(name_array[1], "Zero pol std PA (rad):", zero_pol_pa, "Zero pol std PA (%):", zero_pol_pa*(180/3.142))
    print(name_array[2], "High pol std PA (rad):", high_pol_pa , "High pol PA (%):", high_pol_pa*(180/3.142) )
    
    return (target_pa, zero_pol_pa, high_pol_pa)
#
#
#
#
#
#These are rather deprecated

def plot_PA_stability(input_data, q_u_check, plot_verbose , mjd_align_check, verbose, sv_im):
    #MJD PLOT SHOULD BE ITS OWN FUNCTION AYE
    print("Plot Position Angle stability")

    means_arr = []
    means_err_arr = []
    array_des_dates = []
    
    array_des_noms = []
    mjd_strs = []
    
    for dats in input_data:        
        result = re.search('(.*)_', list(dats.keys())[0][:12])
        targ_name_str = re.search('_(.*)', list(dats.keys())[0])

        array_des_noms.append(targ_name_str.group(1))
        mjd_strs.append(result.group(1)+'T00:00:00.000000000')
        array_des_dates.append(result.group(1)) #This is 
        
        
        mean_q_err_squared = np.mean(dats[list(dats.keys())[0]][1][1:])**2
        mean_u_err_squared = np.mean(dats[list(dats.keys())[0]][3][1:])**2
        
        mean_q_squared = np.mean(dats[list(dats.keys())[0]][0][1:])**2
        mean_u_squared = np.mean(dats[list(dats.keys())[0]][2][1:])**2
        
        mean_q = np.mean(dats[list(dats.keys())[0]][0][1:])
        mean_u = np.mean(dats[list(dats.keys())[0]][2][1:])
        
        sum_o_squares = mean_q_squared + mean_u_squared
        
        #So what do we need?
        #For RINGO Slowikoska et al
        pol_pa = 0.5*math.atan2(mean_u , mean_q )
        pol_pa_err = math.sqrt(((1/(2*mean_q*(1 + (mean_u_squared/mean_q_squared))))**2 )*(mean_u_err_squared) + ((-1*( (mean_u)/(2*sum_o_squares)))**2 )*(mean_q_err_squared))
        
        if(verbose and mjd_align_check):
            print(targ_name_str.group(1), Time([result.group(1)+'T00:00:00.000000000'])[0].mjd  )
        elif(verbose):
            print(targ_name_str.group(1), result.group(1))
            
        if(q_u_check=='degree'):
            if(verbose):
                print("Pol PA!!!", pol_pa*(180/3.142),u"\u00B1",pol_pa*(180/3.142) )

                means_arr.append(  pol_pa*(180/3.142)   ) #mean of qs. Calculate PD
                means_err_arr.append(pol_pa_err*(180/3.142)) #std of list. 
        elif(q_u_check=='rad'):
            if(verbose):
                print("Pol PA!!!", pol_pa,u"\u00B1",pol_pa_err)

                means_arr.append(pol_pa)
                means_err_arr.append(pol_pa_err)
            
    t = Time(mjd_strs, format='isot', scale='utc')
    
    if(mjd_align_check):
        for a in range(0, len(input_data)):
            print(list(input_data[a].keys())[0], t[a], t.mjd[a])

    fig, ax = plt.subplots()
    markers, caps, bars = ax.errorbar(t.mjd, means_arr, yerr=means_err_arr, xerr =[0]*len(means_arr),
            fmt='o', ecolor='blue',capsize=2, capthick=2)
    plt.title(array_des_noms[0]+" Position Angle stability", fontsize=24)
    if(plot_verbose):
        for l in range(0, len(t.mjd)):
            plt.text(t.mjd[l], means_arr[l], str(round(means_arr[l], 2))+u"\u00B1"+str(np.round(means_err_arr[l],4)), rotation = 45, fontsize=24)
    plt.yticks(fontsize = 22)
    plt.xticks(fontsize = 22)
    if(q_u_check=='q'):
        plt.ylabel('q', fontsize=24)
    elif(q_u_check=='u'):
        plt.ylabel('u', fontsize=24)
    plt.axhline(y=0, color = 'black')
    plt.xlabel('MJD',fontsize=24)
    plt.grid()

    [bar.set_alpha(0.22) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]
    
    if(sv_im != ''):
        plt.savefig(sv_im,bbox_inches='tight',pad_inches=0.1 )

    plt.show()    

def plot_PD_stability(input_data, q_u_check, plot_verbose , mjd_align_check, verbose, sv_im):
    #MJD PLOT SHOULD BE ITS OWN FUNCTION AYE
    print("Plot Polarization Degree stability")

    means_arr = []
    means_err_arr = []
    array_des_dates = []
    
    array_des_noms = []
    mjd_strs = []
    
    for dats in input_data:        
        result = re.search('(.*)_', list(dats.keys())[0][:12])
        targ_name_str = re.search('_(.*)', list(dats.keys())[0])

        array_des_noms.append(targ_name_str.group(1))
        mjd_strs.append(result.group(1)+'T00:00:00.000000000')
        array_des_dates.append(result.group(1)) #This is 
        
        
        mean_q_squared = np.mean(dats[list(dats.keys())[0]][0][1:])**2
        mean_u_squared = np.mean(dats[list(dats.keys())[0]][2][1:])**2
        
        mean_q_err_squared = np.mean(dats[list(dats.keys())[0]][1][1:])**2
        mean_u_err_squared = np.mean(dats[list(dats.keys())[0]][3][1:])**2
        
        sum_o_squares = mean_q_squared + mean_u_squared
        
        pol_d = math.sqrt(sum_o_squares)
        #From Slowikowska et al 2016
        pol_d_err  = math.sqrt( ((mean_q_squared)/(sum_o_squares))*(mean_q_err_squared)+ ((mean_u_squared)/(sum_o_squares))*(mean_u_err_squared))

        if(verbose and mjd_align_check):
            print(targ_name_str.group(1), Time([result.group(1)+'T00:00:00.000000000'])[0].mjd  )
        elif(verbose):
            print(targ_name_str.group(1), result.group(1))
            
        if(q_u_check=='perc'):
            if(verbose):
                print("Pol D!!!", pol_d*100,u"\u00B1",pol_d_err*100 )
                #print(np.mean(dats[list(dats.keys())[0]][0][1:]),u"\u00B1",np.std(dats[list(dats.keys())[0]][0][1:]))
                means_arr.append(  pol_d*100   ) #mean of qs. Calculate PD
                means_err_arr.append(pol_d_err*100) #std of list. 
        elif(q_u_check=='nperc'):
            if(verbose):
                print("Pol D!!!", pol_d,u"\u00B1",pol_d_err)
                #print(np.mean(dats[list(dats.keys())[0]][2][1:]),u"\u00B1",np.std(dats[list(dats.keys())[0]][2][1:]))
                means_arr.append(pol_d)
                means_err_arr.append(pol_d_err)
            
    t = Time(mjd_strs, format='isot', scale='utc')
    
    if(mjd_align_check):
        for a in range(0, len(input_data)):
            print(list(input_data[a].keys())[0], t[a], t.mjd[a])

    fig, ax = plt.subplots()
    markers, caps, bars = ax.errorbar(t.mjd, means_arr, yerr=means_err_arr, xerr =[0]*len(means_arr),
            fmt='o', ecolor='blue',capsize=2, capthick=2)
    plt.title(array_des_noms[0]+" Polarization Degree stability", fontsize=24)
    if(plot_verbose):
        for l in range(0, len(t.mjd)):
            plt.text(t.mjd[l], means_arr[l], str(round(means_arr[l], 2))+u"\u00B1"+str(np.round(means_err_arr[l],4)), fontsize=24)
    plt.yticks(fontsize = 22)
    plt.xticks(fontsize = 22)
    if(q_u_check=='q'):
        plt.ylabel('q', fontsize=24)
    elif(q_u_check=='u'):
        plt.ylabel('u', fontsize=24)
    plt.axhline(y=0, color = 'black')
    plt.xlabel('MJD',fontsize=24)
    plt.grid()

    [bar.set_alpha(0.22) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]
    
    if(sv_im != ''):
        plt.savefig(sv_im,bbox_inches='tight',pad_inches=0.1 )

    plt.show()    

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


def plot_q_u_stability(input_data, q_u_check, plot_verbose , mjd_align_check, verbose, sv_im):
    #MJD PLOT SHOULD BE ITS OWN FUNCTION AYE
    print("Plot q or u stability")

    means_arr = []
    means_err_arr = []
    array_des_dates = []
    
    array_des_noms = []
    mjd_strs = []
    
    for dats in input_data:        
        result = re.search('(.*)_', list(dats.keys())[0][:12])
        targ_name_str = re.search('_(.*)', list(dats.keys())[0])

        array_des_noms.append(targ_name_str.group(1))
        mjd_strs.append(result.group(1)+'T00:00:00.000000000')
        array_des_dates.append(result.group(1)) #This is 
        
        if(verbose and mjd_align_check):
            print(targ_name_str.group(1), Time([result.group(1)+'T00:00:00.000000000'])[0].mjd  )
        elif(verbose):
            print(targ_name_str.group(1), result.group(1))
        
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
            
    t = Time(mjd_strs, format='isot', scale='utc')
    
    if(mjd_align_check):
        for a in range(0, len(input_data)):
            print(list(input_data[a].keys())[0], t[a], t.mjd[a])
            
    fig, ax = plt.subplots()
    markers, caps, bars = ax.errorbar(t.mjd, means_arr, yerr=means_err_arr, xerr =[0]*len(means_arr),
            fmt='o', ecolor='blue',capsize=2, capthick=2)
    
    plt.title(array_des_noms[0]+" "+q_u_check + " stability", fontsize=24)
    if(plot_verbose):
        for l in range(0, len(t.mjd)):
            plt.text(t.mjd[l], means_arr[l], str(round(means_arr[l], 6))+u"\u00B1"+str(np.round(means_err_arr[l],4)), fontsize=24)
            
    plt.yticks(fontsize = 22)
    plt.xticks(fontsize = 22)
    if(q_u_check=='q'):
        plt.ylabel('q', fontsize=24)
        plt.axhline(np.mean(means_arr), linestyle='dashed', alpha = 0.45 )
        locs, labels = plt.yticks() 

        plt.text(t.mjd[0],
                 np.mean(means_arr) - (locs[1]-locs[0])/6 ,
                 'Stokes q Instrumental Polarization:'+ str(np.round(np.mean(means_arr),4))+u"\u00B1"+str(np.round(np.mean(means_err_arr), 4  ) ),
                 fontsize=24)
    elif(q_u_check=='u'):
        plt.ylabel('u', fontsize=24)
        plt.axhline(np.mean(means_arr), linestyle='dashed', alpha = 0.45)
        locs, labels = plt.yticks() 

        plt.text(t.mjd[0],
                 np.mean(means_arr) - (locs[1]-locs[0])/6,
                 'Stokes u Instrumental Polarization:'+ str(np.round(np.mean(means_arr),4))+u"\u00B1"+str(np.round(np.mean(means_err_arr), 4)), 
                 fontsize=24)
        
    plt.axhline(y=0, color = 'black')
    plt.xlabel('MJD',fontsize=24)
    plt.grid()

    [bar.set_alpha(0.22) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]
    
    if(sv_im != ''):
        plt.savefig(sv_im,bbox_inches='tight',pad_inches=0.1 )
    plt.show()
    
    #Here is where you do the interpolate
    
    #bivariate spline interpolation
    tck = interpolate.splrep(t.mjd, means_arr, s=0)
    print("tck:", tck)
    print("tck[0]:", tck[0])
    print("tck[1]:", tck[1])
    
    #simple linear interpolation
    f = interp1d(t.mjd, means_arr)
    f2 = interp1d(t.mjd, means_arr, kind='cubic')
    
    #lets
    #
    plt.plot(tck[0], tck[1])
    plt.errorbar(t.mjd, means_arr, yerr=means_err_arr, xerr =[0]*len(means_arr),
            fmt='o', ecolor='blue',capsize=2, capthick=2)
    plt.grid()
    plt.show()
    
    print(f)
    print(f2)
    
    plt.plot(f(t.mjd))
    plt.plot(f2(t.mjd))
    #plt.errorbar(t.mjd, means_arr, yerr=means_err_arr, xerr =[0]*len(means_arr),
    #        fmt='o', ecolor='blue',capsize=2, capthick=2)
    plt.grid()
    plt.show()
    
    return(  [np.mean(means_arr) , np.mean(means_err_arr)]  )

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

def q_n_u_single_plot_v1(pol_data, plot_c, sv_im_str,  sv_im,MJD_arg ,pol_deg, only_means ,key_verb):
    """
    #please document. Thanks past me...
    """
    #data input, plot color, save image string, arguement for MJD (boolean), polaraization degree draw a line, launch verb (redundant)
    #key_verb is for opennig text
    
    #What is not needed in v2 is the plot_c
    #what is needed in v2 is that sv_im_str
    result = re.search('./img_out/(.*)', sv_im_str)
    
    if(key_verb):
        print("N obs:", len(pol_data))
    #something that I notice that you havent done... Calculate q means std and u means std... easy take it in quadrature.
    
    targ_qmeans = []
    targ_umeans = []
    
    targ_qmeans_err = []
    targ_umeans_err = []
    
    target_qs = []
    target_us = []
    targ_qstds = []
    targ_ustds = []
    
    targ_date_strs = [] #
    mjd_strs = []
    
    for things in pol_data:
        targ_qmeans.append(np.mean(things[list(things.keys())[0]][0][1:])) #means
        targ_umeans.append(np.mean(things[list(things.keys())[0]][2][1:])) #ummeans
        
        targ_qmeans_err.append(np.std(things[list(things.keys())[0]][0][1:])) #means
        targ_umeans_err.append(np.std(things[list(things.keys())[0]][2][1:])) #ummeans
                                                                           
        #yo we can do it here
        date = re.search('(.*)_', list(things.keys())[0][:12]) #suppress the target name, just give me the 
        targ_date_strs.append(date.group(1)) #just the keys ay
        mjd_strs.append(date.group(1)+'T00:00:00.000000000')

        target_qs = target_qs + things[list(things.keys())[0]][0][1:]
        targ_qstds = targ_qstds + things[list(things.keys())[0]][1][1:]
        target_us = target_us + things[list(things.keys())[0]][2][1:] 
        targ_ustds = targ_ustds + things[list(things.keys())[0]][3][1:]
        
    t = Time(mjd_strs, format='isot', scale='utc')
    
    if(only_means):
        plt.scatter(targ_qmeans, targ_umeans, alpha=0.9, color = plot_c,)
        plt.errorbar(targ_qmeans, targ_umeans, xerr=targ_qmeans_err, yerr=targ_umeans_err, color = plot_c, lw=0.75, fmt="o", alpha=0.9)
    else:
        plt.scatter(target_qs, target_us, color = plot_c, alpha=0.11)
        plt.errorbar(target_qs, target_us, xerr=targ_qstds, yerr=targ_ustds, lw=0.75, fmt="o", color=plot_c, alpha=0.1)

        plt.scatter(targ_qmeans, targ_umeans, alpha=0.9)
        plt.errorbar(targ_qmeans, targ_umeans, xerr=targ_qmeans_err, yerr=targ_umeans_err, lw=0.75, fmt="o", alpha=0.9)

    if(pol_deg):
        for z in range(0, len(targ_qmeans)):
            plt.plot([0,targ_qmeans[z]], [0, targ_umeans[z]], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
            if(MJD_arg):
                plt.text(targ_qmeans[z], targ_umeans[z], int(t[z].mjd), rotation=-45, fontsize=10) #suppress the title
            else:
                plt.text(targ_qmeans[z], targ_umeans[z], targ_date_strs[z], rotation=-45, fontsize=10) #suppress the title
    
    plt.grid()
    if(sv_im_str != ''):
        plt.title(result.group(1)+" Pol Scatter")
        plt.savefig(sv_im_str,bbox_inches='tight',pad_inches=0.1 )
    else:
        plt.title("Pol Scatter")        
        
    plt.yticks(fontsize = 22)
    plt.xticks(fontsize = 22)        
    plt.ylabel('u', fontsize = 24)
    plt.xlabel('q', fontsize = 24)
    plt.axhline(y=0, color = 'black')
    plt.axvline(x=0, color = 'black')
    
    if(sv_im != ''):
        plt.savefig(sv_im,bbox_inches='tight',pad_inches=0.1 )
    
    plt.show()
    
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
    
    for things in pol_data:
        print(list(things.keys())[0])#, ,things[list(things.keys())[0]], type(things[list(things.keys())[0]]))

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
    
    #combined_q = target_data[0][1:] + zero_pol_std[0][1:] +  high_pol_std[0][1:]
    #combined_u = target_data[2][1:] + zero_pol_std[2][1:] +  high_pol_std[2][1:]

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
        
def calib_data(inp_data, instrumental_pol, plt_show,verbose):
    """
    #Function that scatter plots q nad u
    #This thing just plots. It does not do anything fance such as compute mean blah
    #Error should be taken in quadrature
    
    #Defunct
    #Not built you
    
    if(stats):
        mean_q, mean_u, median_q, median_u = q_u_stats(q, q_err, u, u_err)
        
    #parse the target name from the data_name
    result = re.search(MJD+'_(.*)_P', data_name)
    #print()
    
    plt.scatter(q, u) #u is y (vertical), q is x (horizontal) 
    plt.title(result.group(1)+" q and u scatter plot")
    plt.xlabel("q")
    plt.ylabel("u")
    
    if(mean_med == 'mean'):
        plt.plot([0,mean_q], [0, mean_u], 'k-', lw=1.75, alpha=0.4, linestyle = '--') #vertical
    elif(mean_med == 'median'):
        plt.plot([0,median_q], [0, median_u], 'k-', lw=1.75, alpha=0.4, linestyle = '--') #vertical
    elif(mean_med == ''):
        plt.plot([0,0], [0, 0], 'k-', lw=1.75, alpha=0.4, linestyle = '--') #vertical
        
    plt.grid()
    plt.show()
    
    if(mean_med == 'mean'):
        zero_pol_calib_factor=(mean_q, mean_u)
    elif(mean_med == 'median'):
        zero_pol_calib_factor=(median_q, median_u    )
    elif(mean_med == ''):
        zero_pol_calib_factor = (0,0)
    
    return zero_pol_calib_factor
    #Ideally you have to recreate the data just as it was on the other end.
    #Work on the artefact
    """
    cal_prod = cp(inp_data) #I have copied the data
    
    print("Function that takes in a dataset and a calibration points and subtracts the data") #calibration point an
    print("I am here!")
    if(verbose):
        print("Data (pre cal):", inp_data)
        print("Instrumental Polarization:", instrumental_pol) 
        
    qs_uncal = []
    us_uncal = []
    
    qs_err = []
    us_err = []
        
    qs_cal = []
    us_cal = []
    
    print("Iterate through calibrated product:")
    
    for k in range(0, len(cal_prod)):
        print(cal_prod[k].keys())
        print(cal_prod[k][list(cal_prod[k].keys())[0]]) #
        print("Len of that things. Expectation 4:", len(cal_prod[k][list(cal_prod[k].keys())[0]]    ))
        for z in range(1, len(cal_prod[k][list(cal_prod[k].keys())[0]][0][:])):
            print(cal_prod[k][list(cal_prod[k].keys())[0]][0][z])
            
    
    #rewrite
    #for j in range(0, len(inp_data)):
        #inp_data[j]
    #    for k in range(0, len(       inp_data[list(inp_data[j].keys())[0]]      )):
    #        print("Replace Values")
    
    #for j in range(0, len(inp_data)):
    #    if(verbose):
    #        print(inp_data[j])
    #    print(  )
        #for h in range(0, len(     inp_data[j][list(inp_data[j].keys())[0]][0][1:]   )):
    
    
    
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
    """
        
        
        
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
            

#Per filter
    
#First step is to calculate q and u for zero polarization standard stars
    
#Second is to calculate the instrumental q and u which is the average q and u for zero pol standard

#calulate q and u for high pol, subtract high pol q and u with instrumental q and u

#calculate position angle. This gives you the PA correction. #That PA offset should be the original - new.
#we expect the atan to change very little between calibration steps

#Apply the calibration factors to the target data. PA correction to the target data.

#Done
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