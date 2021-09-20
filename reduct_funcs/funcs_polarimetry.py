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

def source_peak_finder_pol_std(fits_data_1, siegma, search_array, trim, plot_peaks, verbose):
    """
    Function that finds coordinates of stars.
    Star Finder used for the really faint stars of the pol std group
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
    
def load_pol_data(file_dir):
    loc = (file_dir)

    wb = xlrd.open_workbook(loc)
    sheet = wb.sheet_by_index(0)

    q_arr = []
    q_err_arr = []
    u_arr = []
    u_err_arr = []

    for j in range(0, sheet.nrows ):
        q_arr.append(sheet.cell_value(j, 20))
        q_err_arr.append(sheet.cell_value(j, 21))
        u_arr.append(sheet.cell_value(j, 22))
        u_err_arr.append(sheet.cell_value(j, 23))
        
    return(q_arr,  q_err_arr, u_arr, u_err_arr)
    
def calc_stokes_q_n_u(counts, counts_errors, verbose):
    """
    #Function that calculates Stokes q n u
    
    #I assune that the first and the second are the ordinary rays. The bright one.
    #Third and Fourth are the extraordinary rays.
    """
    print("Function that calculates Stokes q n u")
    print("Input Counts and Error counts in order of target for simplicity")
    print("Does the same thing as the excel file. The order of the input determines the order of target 1, 2, 3 or 4.")
    
    q=(counts[0]-count[3])/(counts[0]+count[3])
    u=(counts[1]-counts[2])/(counts[1]+counts[2]) #Have a think about this matey
    
    q_err = 2
    u_err = 3
    
    print("Norm Q (q)", q, "Norm U (u)", u)

    if(verbose):
        plt.plot(counts)
        plt.grid()
        plt.show()
        
    return [q, q_err, u, u_err] #returns an array of q and u

def plot_PD(pds):
    print("Takes in an array and plots the polarization degree")
    
    plt.plot(pds)
    plt.title("Polarization Degree")
    plt.axhline(y=np.mean(pds))
    plt.grid()
    plt.show()
    
def q_u_stats(q, q_err, u, u_err):
    mean_q = np.mean(q)
    mean_u = np.mean(u)
    
    median_q = np.median(q)
    median_u = np.median(u)
    
    return(mean_q, mean_u, median_q, median_u)

def correct_q_u(target_data, zero_pol_std, high_pol_std, zero_pol_offset):   
    target_data[0][1:] = target_data[0][1:] - zero_pol_offset[0]
    target_data[2][1:] = target_data[2][1:] - zero_pol_offset[1]
    
    zero_pol_std[0][1:] = zero_pol_std[0][1:] - zero_pol_offset[0]
    zero_pol_std[2][1:] = zero_pol_std[2][1:] - zero_pol_offset[1]
    
    high_pol_std[0][1:] = high_pol_std[0][1:] - zero_pol_offset[0]
    high_pol_std[2][1:] = high_pol_std[2][1:] - zero_pol_offset[1]
    
    return(target_data, zero_pol_std, high_pol_std)

def calc_pd(target_data, zero_pol_std, high_pol_std):
    """
    #Just do all?
    """
    target_pd = math.sqrt( np.mean(target_data[0][1:])**2 + np.mean(target_data[2][1:])**2) #This works using scatter bro
    zero_pol_pd = math.sqrt(np.mean(zero_pol_std[0][1:])**2 + np.mean(zero_pol_std[2][1:])**2)
    high_pol_pd = math.sqrt(np.mean(high_pol_std[0][1:])**2 + np.mean(high_pol_std[2][1:])**2)
    
    print("Target PD:", target_pd)
    print("Zero pol PD:", zero_pol_pd)
    print("Target PD:", high_pol_pd)
    
    return (target_pd, zero_pol_pd, high_pol_pd )

def calc_pa(target_data, zero_pol_std, high_pol_std):
    
    target_pa = 0.5*math.atan(np.mean(target_data[0][1:])/np.mean(target_data[2][1:]))
    zero_pol_pa = 0.5*math.atan(np.mean(zero_pol_std[0][1:])/np.mean(zero_pol_std[2][1:]))
    high_pol_pa = 0.5*math.atan(np.mean(high_pol_std[0][1:])/np.mean(high_pol_std[2][1:]))
        
    print("Target PA:", target_pa)
    print("Zero pol std PA:", zero_pol_pa)
    print("High pol std PA:", high_pol_pa)
    
    return (target_pa, zero_pol_pa, high_pol_pa)
    
def q_n_u_stack_plot(target_data, zero_pol_std, high_pol_std, title , pol_deg ):
    """
    #Function that scatter plots q nad u
    #This thing just plots. It does not do anything fance such as compute mean blah
    #Error should be taken in quadrature
    
    #Plot all in 4
    
    """
    #print("Plot all 4!")

    combined_q = target_data[0][1:] + zero_pol_std[0][1:] +  high_pol_std[0][1:]
    combined_u = target_data[2][1:] + zero_pol_std[2][1:] +  high_pol_std[2][1:]

    targ_mean_q, targ_mean_u, median_q, median_u = q_u_stats(target_data[0][1:] ,target_data[1][1:] ,target_data[2][1:] ,target_data[3][1:] )
    zero_pol_mean_q, zero_pol_mean_u, median_q, median_u = q_u_stats(zero_pol_std[0][1:] ,zero_pol_std[1][1:] ,zero_pol_std[2][1:] ,zero_pol_std[3][1:] )
    high_pol_mean_q, high_pol_mean_u, median_q, median_u = q_u_stats(high_pol_std[0][1:] ,high_pol_std[1][1:] ,high_pol_std[2][1:] ,high_pol_std[3][1:] )

    fig, axs = plt.subplots(2, 2)
    fig.suptitle("Polarimetry Scatter Plot "+title)
    axs[0, 0].scatter(target_data[0][1:], target_data[2][1:])
    axs[0, 0].set_title("Target")
    axs[0, 0].grid()
    if(pol_deg):
        axs[0, 0].plot([0,targ_mean_q], [0, targ_mean_u], 'k-', lw=1.75, alpha=0.4, linestyle = '--') #vertical
    axs[1, 0].scatter(zero_pol_std[0][1:], zero_pol_std[2][1:])
    axs[1, 0].set_title("Zero Polarization Standard Star")
    if(pol_deg):
        axs[1, 0].plot([0,zero_pol_mean_q], [0, zero_pol_mean_u], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
    axs[1, 0].grid()
    axs[0, 1].scatter(high_pol_std[0][1:], high_pol_std[2][1:])
    axs[0, 1].set_title("High Polarization Standard Star")
    if(pol_deg):
        axs[0, 1].plot([0,high_pol_mean_q], [0, high_pol_mean_u], 'k-', lw=1.75, alpha=0.4, linestyle = '--')
    axs[0, 1].grid()
    axs[1, 1].scatter(combined_q, combined_u)
    axs[1, 1].set_title("Combined")
    axs[1, 1].grid()
    fig.tight_layout()


def calib_pipe():
    """
    #Function that scatter plots q nad u
    #This thing just plots. It does not do anything fance such as compute mean blah
    #Error should be taken in quadrature
    
    #Defunct
    """
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

#Per filter
    
#First step is to calculate q and u for zero polarization standard stars
    
#Second is to calculate the instrumental q and u which is the average q and u for zero pol standard

#calulate q and u for high pol, subtract high pol q and u with instrumental q and u

#calculate position angle. This gives you the PA correction. #That PA offset should be the original - new.
#we expect the atan to change very little between calibration steps

#Apply the calibration factors to the target data. PA correction to the target data.

#Done