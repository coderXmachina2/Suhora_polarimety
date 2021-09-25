import matplotlib.pyplot as plt
import numpy as np
import astropy
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

import matplotlib.pyplot as plt

# Some style for better looking plots
from pylab import rcParams
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Verdana']
plt.rcParams['font.size'] = 12
plt.rcParams['font.size'] = 12
plt.rcParams['lines.linewidth'] = 4.
plt.rcParams['axes.labelsize'] = 'medium'
plt.rcParams['grid.linewidth'] = 1.0
plt.rcParams['grid.linestyle'] = '-'
plt.rcParams['xtick.minor.size']=4
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.minor.size']=4
plt.rcParams['ytick.major.size']=8
plt.rcParams['figure.figsize'] = 14,14
plt.rcParams['figure.subplot.bottom'] = 0.15
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['xtick.labelsize'] = 10

def make_calib_files(bias, dark, flat):
    """
    Takes in bias, dark, and flat. Returns master bias, master dark, and master flat
    """
    b_frames = []
    drk_frames = []
    normed_bs_fl_frames = []

    for bias_file in bias:
        bias_frame = astropy.io.fits.open(bias_file)
        b_frames.append(bias_frame[0].data)
    
    master_bias = np.median(b_frames, axis= 0) #I am sure that my master bias method is corrct

    for dark_file in dark:
        drk_frame = astropy.io.fits.open(dark_file)
        bias_subtracted_dark = drk_frame[0].data - master_bias
        drk_frames.append(bias_subtracted_dark)
    
    master_drk =  np.median(drk_frames, axis= 0)
    
    #Here is an idea, depending on which pol filter you want to use, only take 
    #the corresponding p filter (P3 - R) or (P1 - R)
    #
    
    for flat_f in flat:
        flat_frame = astropy.io.fits.open(flat_f)
        bias_subtracted_flat = np.subtract(flat_frame[0].data,master_bias)
        normed_bias_subtracted_flat = np.divide(bias_subtracted_flat, np.median(flat_frame[0].data )    )
        normed_bs_fl_frames.append(normed_bias_subtracted_flat)

    master_flat = np.median(normed_bs_fl_frames, axis= 0)
    
    return(master_bias, master_drk, master_flat)

def file_splits(list_data):
    filts= []
    objects = []
    exp_times = []
    
    for k in range(0, len(list_data)):
        if(astropy.io.fits.open(list_data[0])[0].header['IMAGETYP'] == 'object'):
            exp_times.append(  astropy.io.fits.open(list_data[k])[0].header['EXPTIME']   ) 
        filts.append(astropy.io.fits.open(list_data[k])[0].header['FILTER'])
        objects.append(astropy.io.fits.open(list_data[k])[0].header['OBJECT'])

    print(astropy.io.fits.open(list_data[0])[0].header['IMAGETYP'],
          set(objects), 
          "Filters:", [(filters, filts.count(filters)) for filters in set(filts)], #Wow thats some bullshit
          "total " + str(np.round(sum(exp_times)/60)) + " min exposure time" if astropy.io.fits.open(list_data[0])[0].header['IMAGETYP'] == 'object' else '' )

def reduction(raw_data, calib_files):
    """
    Takes in the calibration files and raw data. Takes in two input. Returns reduced fits image.
    
    """
    bs_data = np.subtract(raw_data, calib_files[0]) #subtract the bias
    bs_ds_data = np.subtract(bs_data, calib_files[1]) #subtract the dark
    bs_ds_ff_data = np.divide(bs_ds_data, calib_files[2] ) #divide the flat
                          
    return(bs_ds_ff_data)

def plot_raw_double_compare(fits_data_1, fits_data_2, scale_arr, sigma, plot):
    """
    Plots 2 data side by side. Without reduction. Sigma Clipped stats is an option.
    
    Take is a whole fits file. Does not return anything.
    
    """
                          
    m1 = fits_data_1[0].data[scale_arr[0]-scale_arr[2]:scale_arr[1]+scale_arr[2], 
                          scale_arr[1]-scale_arr[3]:scale_arr[1]+scale_arr[3]]
    m2 = fits_data_2[0].data[scale_arr[0]-scale_arr[2]:scale_arr[1]+scale_arr[2], 
                          scale_arr[1]-scale_arr[3]:scale_arr[1]+scale_arr[3]]

    fig = plt.figure(figsize=(16, 12))
    ax1 = fig.add_subplot(121)
    
    if(sigma):
        mean, median, std = sigma_clipped_stats(m1)
        im1 = ax1.imshow(m1, 
                         vmin = median - 3*std,
                         vmax = median + 3*std,
                         interpolation='None')
    else:
        im1 = ax1.imshow(m1)

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')
    
    ax1.set_title(fits_data_1[0].header['OBJECT']+ " " +fits_data_1[0].header['FILTER'] + " raw " + fits_data_1[0].header['TIME-OBS'])
    
    ax1.axvline(x=scale_arr[2],alpha = 0.15, color='red')
    ax1.axhline(y=scale_arr[2],alpha = 0.15, color='red')
                          
    ax2 = fig.add_subplot(122)
    
    if(sigma):
        mean, median, std = sigma_clipped_stats(m2)
        im2 = ax2.imshow(m2,
                         vmin = median - 3*std,
                         vmax = median + 3*std,
                         interpolation='None')
    else:
        im2 = ax2.imshow(m2)
    
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im2, cax=cax, orientation='vertical')
    
    ax2.set_title(fits_data_2[0].header['OBJECT']+ " " +fits_data_2[0].header['FILTER'] + " raw " + fits_data_2[0].header['TIME-OBS'])
    ax2.axvline(x=scale_arr[2],alpha = 0.15, color='red')
    ax2.axhline(y=scale_arr[2],alpha = 0.15, color='red')
    
    if (plot):
        if(fits_data_1[0].header['OBJECT'] =='eecep'):
            string_out='./img_out/target/'+ fits_data_1[0].header['OBJECT']+'_'+ fits_data_1[0].header['FILTER'] +'_'+fits_data_1[0].header['TIME-OBS'].replace(":", "-")
        else:
            string_out='./img_out/pol_std/'+ fits_data_1[0].header['OBJECT']+'_'+fits_data_1[0].header['FILTER'] +'_'+fits_data_1[0].header['TIME-OBS'].replace(":", "-")
        if(sigma):
            string_out += '_sigma.png'
        else:
            string_out += '.png'
        plt.savefig(string_out ,bbox_inches='tight',
                        pad_inches=0.1)

def plot_double_raw_v_reduced(fits_data_1, scale_arr, calib_files, sigma, file_out, plot):
    """
    Plots 2 data side by side. One data is reduced/ bias, dark subtracted
    Only needs one input. Sigma Clipped stats is an option.
    
    Instead of returning a numpy array image I need it to return a fits_object
    
    """    
    reduced = reduction(fits_data_1[0].data, calib_files)
    
    reduced_fits_obj = cp(fits_data_1)
    reduced_fits_obj[0].data = reduced #so I can overwrite it just like this
                          
    m1 = fits_data_1[0].data[scale_arr[0]-scale_arr[2]:scale_arr[1]+scale_arr[2], 
                             scale_arr[1]-scale_arr[3]:scale_arr[1]+scale_arr[3]]
    m2 = reduced[scale_arr[0]-scale_arr[2]:scale_arr[1]+scale_arr[2], 
                 scale_arr[1]-scale_arr[3]:scale_arr[1]+scale_arr[3]]

    if(plot):
        fig = plt.figure(figsize=(16, 12))

        ax1 = fig.add_subplot(121)

        if(sigma):
            mean, median, std = sigma_clipped_stats(m1)
            im1 = ax1.imshow(m1, 
                             vmin = median - 3*std,
                             vmax = median + 3*std,
                             interpolation='None')
        else:
            im1 = ax1.imshow(m1)

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im1, cax=cax, orientation='vertical')

        ax1.set_title(fits_data_1[0].header['OBJECT']+ " " +fits_data_1[0].header['FILTER'] + " raw " + fits_data_1[0].header['TIME-OBS'])
        ax1.axvline(x=scale_arr[2],alpha = 0.15, color='red')
        ax1.axhline(y=scale_arr[2],alpha = 0.15, color='red')

        ax2 = fig.add_subplot(122)

        if(sigma):
            mean, median, std = sigma_clipped_stats(m2)
            im2 = ax2.imshow(m2, 
                             vmin = median - 3*std,
                             vmax = median + 3*std,
                             interpolation='None')
        else:
            im2 = ax2.imshow(m2)

        divider = make_axes_locatable(ax2)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        fig.colorbar(im2, cax=cax, orientation='vertical')
        ax2.set_title(fits_data_1[0].header['OBJECT']+ " " +fits_data_1[0].header['FILTER'] + " reduced " + fits_data_1[0].header['TIME-OBS'])
        ax2.axvline(x=scale_arr[2],alpha = 0.15, color='red')
        ax2.axhline(y=scale_arr[2],alpha = 0.15, color='red')

        plt.show()
    
    if (file_out):
        if(fits_data_1[0].header['OBJECT'] =='eecep'):
            string_out='./img_out/target/'+ fits_data_1[0].header['OBJECT']+'_'+ fits_data_1[0].header['FILTER'] +'_'+fits_data_1[0].header['TIME-OBS'].replace(":", "-")
        else:
            string_out='./img_out/pol_std/'+ fits_data_1[0].header['OBJECT']+'_'+fits_data_1[0].header['FILTER'] +'_'+fits_data_1[0].header['TIME-OBS'].replace(":", "-")
        if(sigma):
            string_out += '_sigma.png'
        else:
            string_out += '.png'
        plt.savefig(string_out ,bbox_inches='tight',
                        pad_inches=0.1)
    
    
    return(reduced_fits_obj)