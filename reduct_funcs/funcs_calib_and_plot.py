import matplotlib.pyplot as plt
import numpy as np
import astropy
import ccdproc

from copy import deepcopy as cp
from scipy import stats
from astropy import units as u
from astropy.nddata import CCDData

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
from reduct_funcs import funcs_utils
from reduct_funcs import funcs_polarimetry

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

def make_calib_files(bias, dark, flat_p1, flat_p3):
    """
    Takes in bias, dark, and flat. Returns master bias, master dark, and master flat
    """
    b_frames = []
    drk_frames = []
    normed_bs_p1fl_frames = []
    normed_bs_p3fl_frames = []

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
    #the corresponding p filter (P3 - R) or (P1 - R). The flats also subtract master_bias
    
    for flat_x in flat_p1:
        flat_frame = astropy.io.fits.open(flat_x)
        bias_subtracted_flat = np.subtract(flat_frame[0].data, master_bias)
        normed_bias_subtracted_flat = np.divide(bias_subtracted_flat, np.median(flat_frame[0].data))
        normed_bs_p1fl_frames.append(normed_bias_subtracted_flat)
        
    for flat_y in flat_p3:
        flat_frame = astropy.io.fits.open(flat_y)
        bias_subtracted_flat = np.subtract(flat_frame[0].data,master_bias)
        normed_bias_subtracted_flat = np.divide(bias_subtracted_flat, np.median(flat_frame[0].data))
        normed_bs_p3fl_frames.append(normed_bias_subtracted_flat)

    masterf_p1 = np.median(normed_bs_p1fl_frames, axis= 0)
    masterf_p3 = np.median(normed_bs_p3fl_frames, axis= 0)    
    
    return(master_bias, master_drk, masterf_p1, masterf_p3)

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
    There should be more steps need to be taken with CCD proc.
    """
    
    #Everything is here. subtract bias, subtract dark, and flat fielding.
    #what is raw data
    #Its not a numpy thing. It is an image.
    #print(type(raw_data))
    #print("Raw data:", raw_data)
    
    bs_data = np.subtract(raw_data, calib_files[0]) #subtract the bias
    bs_ds_data = np.subtract(bs_data, calib_files[1]) #subtract the dark
    bs_ds_ff_data = np.divide(bs_ds_data, calib_files[2] ) #divide the flat
                          
    return(bs_ds_ff_data)

def reduction_ccd_proc(unredu_fits_data, calib_files, key):
    """
    Takes in the calibration files and raw data. Takes in two input. Returns reduced fits image.
    There should be more steps need to be taken with CCD proc.
    """
    
    #It is the same image that pops out in the end. I acknowledge that the change has been applied
    #
    shad_cop = cp(astropy.io.fits.open(unredu_fits_data)) #Shadow copy
    #Everything is here. subtract bias, subtract dark, and flat fielding.
    #what is raw data
    #Its not a numpy thing. It is an image.
    #calib_files = (bias, dark, p_flat)
    
    masterf = CCDData(calib_files[2], unit=u.adu)
    masterf.header = astropy.io.fits.open(unredu_fits_data)[0].header 

    master_b = CCDData(calib_files[0], unit=u.adu)
    master_b.header['exposure'] = astropy.io.fits.open(unredu_fits_data)[0].header['EXPTIME']

    master_d = CCDData(calib_files[1], unit=u.adu)
    master_d.header['exposure'] = astropy.io.fits.open(unredu_fits_data)[0].header['EXPTIME'] 
    
    ccd_data = CCDData(astropy.io.fits.open(unredu_fits_data)[0].data, unit=u.adu) 
    
    ccd_data.header = astropy.io.fits.open(unredu_fits_data)[0].header
    ccd_data.header['exposure'] = ccd_data.header['EXPTIME']
    
    if(key=='2020'):
        #Aspen
        #AspenCG47:
        #gain: 1.15 e/ADU   
        #readout noise: 48.9 e
        
        #Gain correction
        print("Apply Gain Correction here. New Stuff LoL XD. Today. Do target again.")
        
        #Method 1 subtract biast, subtract dark, gain correct and flat field together. no night and day difference
        #bias_subtracted = ccdproc.subtract_bias(ccd_data, master_b)
        #dark_subtracted = ccdproc.subtract_dark(bias_subtracted, master_d,
        #                                        exposure_time='exposure',
        #                                        exposure_unit=u.second)
        #nccd = ccdproc.ccd_process(dark_subtracted,master_flat = masterf,
        #                           gain= 1.15*u.electron/u.adu, 
        #                           readnoise= 48.9*u.electron,                              
        #                           exposure_unit=u.second, 
        #                           exposure_key='exposure')
        
        #Method 2 subtract biast, subtract dark, flat field, and gain correct. no night and day difference
        bias_subtracted = ccdproc.subtract_bias(ccd_data, master_b)
        dark_subtracted = ccdproc.subtract_dark(bias_subtracted, master_d,
                                                exposure_time='exposure',
                                                exposure_unit=u.second)
        
        reduced_image = ccdproc.flat_correct(dark_subtracted, masterf)
        nccd = ccdproc.ccd_process(dark_subtracted,
                                   gain= 1.15*u.electron/u.adu, 
                                   readnoise= 48.9*u.electron,                              
                                   exposure_unit=u.second, 
                                   exposure_key='exposure')
        
        #Method 3 gain correct First to produce with data and deviation, then subtract bias, and dark and flat field.
        
        #data_with_deviation = ccdproc.create_deviation(ccd_data, 
        #                          gain=1.15 * u.electron/u.adu,
        #                          readnoise=48.9 * u.electron)

        #gain_corrected = ccdproc.gain_correct(data_with_deviation, 1.15*u.electron/u.adu) #there is an uncertainty propagated.
        
        #bs_data = np.subtract(gain_corrected.data, calib_files[0]) #subtract the bias
        #bs_ds_data = np.subtract(bs_data, calib_files[1]) #subtract the dark
        #bs_ds_ff_data = np.divide(bs_ds_data, calib_files[2] ) #divide the flat
        
    elif(key=='2021'):
        #Altau
        bias_subtracted = ccdproc.subtract_bias(ccd_data, master_b)
        dark_subtracted = ccdproc.subtract_dark(bias_subtracted, master_d,
                                                exposure_time='exposure',
                                                exposure_unit=u.second)
        #Flat fielding
        
        print("Commencing gain correction:")

        nccd = ccdproc.ccd_process(ccd_data, 
                                   gain= 1.5*u.electron/u.adu, 
                                   readnoise= 8.5*u.electron, 
                                   exposure_unit=u.second, 
                                   exposure_key='exposure') #and put in the values 
        
    #make the new one become the shadcop. return the shadcop
    shad_cop.data = nccd.data
                          
    return(shad_cop)

def plot_raw_double_compare(fits_data_1, fits_data_2, scale_arr, comp_what, sup_tit ,sigma ,plot):
    """
    Plots 2 data side by side. Without reduction. Sigma Clipped stats is an option.
    
    Take is a whole fits file. Does not return anything.
    
    """
    
    #So which one is whcih? Compare with your old one ish!
    #It is with the [0] as per the original
    m1 = fits_data_1.data[scale_arr[0]-scale_arr[2]:scale_arr[1]+scale_arr[2],  #So the common syntax is file[0].data[blah]
                          scale_arr[1]-scale_arr[3]:scale_arr[1]+scale_arr[3]] #
    m2 = fits_data_2.data[scale_arr[0]-scale_arr[2]:scale_arr[1]+scale_arr[2], 
                          scale_arr[1]-scale_arr[3]:scale_arr[1]+scale_arr[3]]

    fig = plt.figure(figsize=(16, 12))
    ax1 = fig.add_subplot(121)
    fig.suptitle(sup_tit)
    
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
    
    ax1.set_title(fits_data_1[0].header['OBJECT']+ " " +fits_data_1[0].header['FILTER']+ ' ' + comp_what[0]+ ' ' + fits_data_1[0].header['TIME-OBS'])
    
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
    
    ax2.set_title(fits_data_2[0].header['OBJECT']+ " " +fits_data_2[0].header['FILTER']+ ' ' + comp_what[1] + ' ' + fits_data_2[0].header['TIME-OBS'])
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
        
    fig.tight_layout()
    plt.show()
    #Does not return anything
        
def plot_double_raw_v_reduced(fits_data_1, scale_arr, calib_files, sup_tit, sigma, file_out, plot):
    """
    Plots 2 data side by side. One data is reduced/ bias, dark subtracted
    Only needs one input. Sigma Clipped stats is an option.
    
    Instead of returning a numpy array image I need it to return a fits_object
    
    """    

    reduced = reduction(fits_data_1[0].data, calib_files) #does the redictopm
    
    #you cloned it and you made the 
    reduced_fits_obj = cp(fits_data_1)
    reduced_fits_obj[0].data = reduced #so I can overwrite it just like this
                          
    m1 = fits_data_1[0].data[scale_arr[0]-scale_arr[2]:scale_arr[1]+scale_arr[2], 
                             scale_arr[1]-scale_arr[3]:scale_arr[1]+scale_arr[3]]
    m2 = reduced[scale_arr[0]-scale_arr[2]:scale_arr[1]+scale_arr[2], 
                 scale_arr[1]-scale_arr[3]:scale_arr[1]+scale_arr[3]]

    if(plot):
        fig = plt.figure(figsize=(16, 12))
        if(sup_tit != ''):
            fig.suptitle(sup_tit)
            
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

        fig.tight_layout()
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
    
    
    return(reduced_fits_obj) #returns this


def calib_pipe(target_data, zero_pol_data, plot_me, key_verb_t):
    """
    #Copy and paste the latest applicable calibration process here.
    #Takes in target data and zero pol data and returns calibrated data
    """
    ##Just get g191. That is not very robust isn't it...
    cal_prod = cp(target_data)
    data_strs = ['g191b2b', 'G191B2B']
    G191_low_pol = funcs_utils.filter_data(zero_pol_data, data_strs, False)

    if(plot_me):
        funcs_polarimetry.q_n_u_single_plot_v1(G191_low_pol, 'blue','./img_out/all_zero_pols','zero_pol',True, True, True, False, False) 

    #this guy is suppressed. For calculation only
    q_cal = funcs_polarimetry.plot_q_u_stability(G191_low_pol, 'q', './img_out/stability/EE Cephei_q_stab', True,False, False, False)
    u_cal = funcs_polarimetry.plot_q_u_stability(G191_low_pol, 'u', './img_out/stability/EE Cephei_u_stab', True,False, False, False)
    
    #This does that removal of points... how about that
    del G191_low_pol[1]
    del G191_low_pol[-1]
    
    if(plot_me):
        funcs_polarimetry.q_n_u_single_plot_v1(G191_low_pol[:],'blue','./img_out/all_zero_pols','zero_pol',True,True,True,False,False)
        
    #this guy is suppressed. For calculation only
    q_cal = funcs_polarimetry.plot_q_u_stability(G191_low_pol[:], 'q', './img_out/stability/EE Cephei_q_stab', True,False, False,False)
    u_cal = funcs_polarimetry.plot_q_u_stability(G191_low_pol[:], 'u', './img_out/stability/EE Cephei_u_stab', True,False, False,False)

    #print("For all 0 pols. q n u instrumental points:")
    #print("q inst:", q_cal[0], u"\u00B1",q_cal[1])
    #print("u inst:", u_cal[0], u"\u00B1",u_cal[1])
    #print("\n")    
    
    #0-2: 1
    #2-4: 2
    #5-13: 3
    #13-22: 4
    #22-len(target_data): 5
    
    c = 1
    cal_c = 0
    arr_qcal = []
    cal_targ = []
    targ_data_arr = [cal_prod[0:2],  cal_prod[2:4], cal_prod[4:13], cal_prod[13:22], cal_prod[22:]]
    targ_slice = ['0-2', '2-4','5-13', '13-22', '22-len(target_data)']
    
    for k in range(0, len(G191_low_pol)):
        q_cal = funcs_polarimetry.plot_q_u_stability(G191_low_pol[0:c], 'q', './img_out/stability/EE Cephei_q_stab',  True,False, False, False)
        u_cal = funcs_polarimetry.plot_q_u_stability(G191_low_pol[0:c], 'u', './img_out/stability/EE Cephei_u_stab', True,False, False, False)

        cal_section = funcs_polarimetry.calib_data(targ_data_arr[k], (q_cal, u_cal), False, False)
        cal_targ = cal_targ + cal_section
        the_slice = [list(x.keys())[0] for x in G191_low_pol[0:c]]
        if(key_verb_t):
            print("For 0 pols:", the_slice,the_slice[0])
            print("And targets:", targ_slice[k])
            print("q cal:", q_cal[0], u"\u00B1", q_cal[1])
            print("u cal:", u_cal[0], u"\u00B1", u_cal[1])
            print("\n")
        c =c + 1
                          
    return(cal_targ)