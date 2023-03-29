import matplotlib.pyplot as plt
import numpy as np
import astropy
import ccdproc
import glob

from copy import deepcopy as cp
from astropy import units as u
from astropy.nddata import CCDData

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

def make_calib_files(bias, 
                     dark, 
                     flat_p1, 
                     flat_p3,
                     MJD):
    """
    A function that takes in lists of bias, dark, p1 and p3 FITS files. Creates master bias, master dark, and master flat by taking the median of the data and returning result.

    Parameters
    ----------
        bias : list
            List of bias FITS file names with their relative directory path. 
        dark : list
            List of dark FITS files names with their relative directory path.
        flat_p1 : list
            List of P1 flat FITS files names with  with their relative directory path.
        flat_p3 : list
            List of P3 flat FITS files names with their relative directory path.
        flat_dir : str
            String that specifies flat
    """

    b_frames = []
    drk_frames = []
    normed_bs_p1fl_frames = []
    normed_bs_p3fl_frames = []
    
    for bias_file in bias:
        bias_frame = astropy.io.fits.open(bias_file)
        b_frames.append(bias_frame[0].data)
    
    master_bias = np.median(b_frames, axis= 0)

    for dark_file in dark:
        drk_frame = astropy.io.fits.open(dark_file)
        bias_subtracted_dark = drk_frame[0].data - master_bias
        drk_frames.append(bias_subtracted_dark)
    
    master_drk =  np.median(drk_frames, axis= 0)
    
    #Get all flats within ± 2 the directories
    
    x = glob.glob('./files_sorted/*')
    data_dirs = []
    print("Search index:", x.index('./files_sorted\\'+MJD))
    for l in range(0, len(x)):
        if(l == x.index('./files_sorted\\'+MJD)):
            print(x[l], "<--- Main")
            data_dirs.append(x[l])
        elif(l == x.index('./files_sorted\\'+MJD) + 1 or 
             l == x.index('./files_sorted\\'+MJD) - 1 or
             l == x.index('./files_sorted\\'+MJD) + 2 or
             l == x.index('./files_sorted\\'+MJD) - 2):
            print(x[l], "<-")
            data_dirs.append(x[l])
        else:
            print(x[l])
            
    flat_p1 = []
    flat_p3 = []
    
    for dirs in data_dirs:
        print(dirs+'flat/*'   )
        gdirs = glob.glob(dirs+'/flat/*')
        print(gdirs)
        for flatf in gdirs:
            if 'p1' in flatf:
                print(flatf, "<- p1 flat")
                flat_p1.append(flatf)
            elif 'p3' in flatf:
                print(flatf, "<- p3 flat")
                flat_p3.append(flatf)        

           
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
    """
    A function that takes in a list of any data type (bias, dark, flat or target) and iterates through the list to get basic information on the image type, filters, and exposure time if the object type not a reduction file. Used for initial checking and was initially designed for EDA. Expects single type in the input.

    Parameters
    ----------
        list_data : list
            List of either bias, dark, flat, or data FITS files with their relative directory path. Loops through the list and adds FITS metadata from the header into separate list. Prints the information with some formatting.
    """
    img_type= []
    filts= []
    objects = []
    exp_times = []
    
    for k in range(0, len(list_data)):
        if(astropy.io.fits.open(list_data[0])[0].header['IMAGETYP'] == 'object'):
            exp_times.append(astropy.io.fits.open(list_data[k])[0].header['EXPTIME']   ) 
        filts.append(astropy.io.fits.open(list_data[k])[0].header['FILTER'])
        objects.append(astropy.io.fits.open(list_data[k])[0].header['OBJECT'])
        img_type.append(astropy.io.fits.open(list_data[k])[0].header['IMAGETYP'])

    print("Img type:", [x for x in list(set(img_type))],
          "Obj:", [t for t in list(set(objects))], 
          "Filters:", [(filters, filts.count(filters)) for filters in set(filts)],
          "\nTotal " + str(np.round(sum(exp_times)/60)) + " min exposure time" if astropy.io.fits.open(list_data[0])[0].header['IMAGETYP'] == 'object' else '' )

def reduction(raw_data, 
              calib_files):
    """
    A function that takes in openned FITS image data and performs reduction with raw numpy routines (direct quick and dirty computation). Image subtract bias, subtract dark, and divide flat. Returns 2d nd array data.

    Parameters
    ----------
        raw_data : numpy ndarray
            Single numpy ndarray of image intensities. 2D image array.
        calib_files : tuple
            Tuple of length 3 comprising of bias (calib_files[0]), dark (calib_files[1]), and flat (calib_files[1]).
    """
    
    bs_data = np.subtract(raw_data, calib_files[0]) #subtract the bias
    bs_ds_data = np.subtract(bs_data, calib_files[1]) #subtract the dark
    bs_ds_ff_data = np.divide(bs_ds_data, calib_files[2] ) #divide the flat
                          
    return(bs_ds_ff_data)

def reduction_ccd_proc(unredu_fits_data, 
                       calib_files, 
                       key=''):
    """
    A function that takes in FITS file data and performs reduction with ccd proc routines. Image subtract bias, subtract dark, and divide flat. Was still in testing and not used in deployment for EAS 2022. As of 21/12/2022 still not deployed and in testing but last I recall the results were approximately similar with quick and dirty method. Followup to do another round of verification and deploy.Returns 2D nd array data.

    Parameters
    ----------
        unredu_fits_data : str
            Single FITS file name with their relative directory path
        calib_files : tuple,
            Tuple of length 3 comprising of bias (calib_files[0]), dark (calib_files[1]), and flat (calib_files[1]).
        key: str
            Year '2020' or '2021'. Conducts different reduction based on camera used to capture data assuming the correct input files are given.
    """
    
    shad_cop = cp(astropy.io.fits.open(unredu_fits_data))
    
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
               
        #Method 1 subtract bias, subtract dark, gain correct and flat field together. No night and day difference
        #bias_subtracted = ccdproc.subtract_bias(ccd_data, master_b)
        #dark_subtracted = ccdproc.subtract_dark(bias_subtracted, master_d,
        #                                        exposure_time='exposure',
        #                                        exposure_unit=u.second)
        #nccd = ccdproc.ccd_process(dark_subtracted,master_flat = masterf,
        #                           gain= 1.15*u.electron/u.adu, 
        #                           readnoise= 48.9*u.electron,                              
        #                           exposure_unit=u.second, 
        #                           exposure_key='exposure')
        
        #Method 2 subtract bias, subtract dark, flat field, and gain correct. no night and day difference
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
        reduced_image = ccdproc.flat_correct(dark_subtracted, masterf)
        nccd = ccdproc.ccd_process(ccd_data, 
                                   gain= 1.5*u.electron/u.adu, 
                                   readnoise= 8.5*u.electron, 
                                   exposure_unit=u.second, 
                                   exposure_key='exposure') #and put in the values 
    else:
        print("Invalid Camera Key")
        
    #make the new one become the shadcop. return the shadcop
    shad_cop.data = nccd.data
                          
    return(shad_cop)

def plot_raw_double_compare(fits_data_1, 
                            scale_arr, 
                            comp_what = ["plot A", "plot B"], 
                            sv_img = False):
    """
    A function that plots FITS image data. Does not return anything

    Parameters
    ----------
        fits_data_1 : str
            Single numpy ndarray of image intensities. 2D image array.
        scale_arr : list
            List of integers. Array that scales both image data to act as zoom. 
        comp_what :  list
            A list that contains image sub titles.
        sv_img : bool, optional
            Saves image to file. Expects the correct directory to already exist.
    """
    
    op = astropy.io.fits.open(fits_data_1)[0] #header abd data
    
    m1 = op.data[scale_arr[0]-scale_arr[2]:scale_arr[1]+scale_arr[2],  
                          scale_arr[1]-scale_arr[3]:scale_arr[1]+scale_arr[3]]

    fig = plt.figure(figsize=(16, 12))
    ax1 = fig.add_subplot(121)

    im1 = ax1.imshow(m1)

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')
    
    ax1.set_title(op.header['OBJECT']+ " " +op.header['FILTER']+ ' ' + comp_what[0]+ ' ' + op.header['TIME-OBS'])
    
    ax1.axvline(x=scale_arr[2],alpha = 0.15, color='red')
    ax1.axhline(y=scale_arr[2],alpha = 0.15, color='red')
                          
    ax2 = fig.add_subplot(122)  

    mean, median, std = sigma_clipped_stats(m1)
    m1 = ax2.imshow(m1, vmin = median - 3*std,
                         vmax = median + 3*std,
                         interpolation='None')
    
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')
    
    ax2.set_title(op.header['OBJECT']+ " " +op.header['FILTER']+ ' ' + comp_what[1] + ' ' + op.header['TIME-OBS'])
    ax2.axvline(x=scale_arr[2],alpha = 0.15, color='red')
    ax2.axhline(y=scale_arr[2],alpha = 0.15, color='red')
    
    if (sv_img):
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
        
    fig.tight_layout(pad=0.6, w_pad=0.5, h_pad=1.0)
    plt.show()
        
def plot_double_raw_v_reduced(fits_data_1, 
                              calib_files,
                              scale_arr, 
                              sigma=False, 
                              plot=True, 
                              sv_img=False):
    """
    A function that plots fits image data. Does not return anything

    Parameters
    ----------
        fits_data_1 : str
            Single FITS file name with their relative directory path
        calib_files : tuple
            Tuple of length 3 comprising of bias (calib_files[0]), dark (calib_files[1]), and flat (calib_files[1])
        scale_arr : list
            List of integers. Array that scales both image data to act as zoom. 
        sigma : bool, optional
            Applies sigma_clipped_stats to image. Sigma is 3 by default
        plot : bool, optional
            Plot image. True by default.
        sv_img : bool, optional
            Saves image to file.
    """
    
    op = astropy.io.fits.open(fits_data_1)[0] 
    reduced = reduction(op.data, calib_files)     
    reduced_fits_data = cp(reduced)
    
    m1 = op.data[scale_arr[0]-scale_arr[2]:scale_arr[1]+scale_arr[2], 
                             scale_arr[1]-scale_arr[3]:scale_arr[1]+scale_arr[3]]
    m2 = reduced_fits_data [scale_arr[0]-scale_arr[2]:scale_arr[1]+scale_arr[2], 
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

        ax1.set_title(op.header['OBJECT']+ " " +op.header['FILTER'] + " raw " + op.header['TIME-OBS'])
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
        
        ax2.set_title(op.header['OBJECT']+ " " +op.header['FILTER'] + " reduced " + op.header['TIME-OBS'])
        ax2.axvline(x=scale_arr[2],alpha = 0.15, color='red')
        ax2.axhline(y=scale_arr[2],alpha = 0.15, color='red')

        fig.tight_layout()
        plt.show()
    
    if (sv_img):
        if(op.header['OBJECT'] =='eecep'): #What if it is spelled differently
            string_out='./img_out/target/'+ op.header['OBJECT']+'_'+ op.header['FILTER'] +'_'+fits_data_1[0].header['TIME-OBS'].replace(":", "-")
        else:
            string_out='./img_out/pol_std/'+ op.header['OBJECT']+'_'+op.header['FILTER'] +'_'+op.header['TIME-OBS'].replace(":", "-")
        if(sigma):
            string_out += '_sigma.png'
        else:
            string_out += '.png'
        plt.savefig(string_out,
                    bbox_inches='tight',
                    pad_inches=0.1)
       
    return(reduced_fits_data)

def calib_data(inp_data, 
               instrumental_pol, 
               plt_show = False,
               verbose=False):
    """
    A function that takes in data and value for instrumental polarization and computes calibration shift. Returns standard data artefact.

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
    
    if(verbose):
        print("Calibrating data...") #calibration point an
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


def calib_pipe(input_data, 
               zero_pol_data, 
               list_index_unstable_data=[],
               verbose_plot_zpol=False,
               verbose_plot_points=False,
               verbose_text=False):
    """
    A function that implements an experimental draft pipeline. Copy and pasted

    Parameters
    ----------
        input_data : list
            Tuple containing a list of dictionaries with q, q error, u, u error data, and a list of date time objects (data timestamps)
        zero_pol_data : list
            Tuple of length 3 comprising of bias (calib_files[0]), dark (calib_files[1]), and flat (calib_files[1]).
        list_index_unstable_data : list
            List of indexes of the zero polarization points to be removed. Be warned, remove one and the index shifts. Practice in notebook first then proceed.
        verbose_plot_zpol : bool, optional
            List of integers. Array that scales both image data to act as zoom. 
        verbose_plot_points : bool, optional
            List of integers. Array that scales both image data to act as zoom. 
        key_verb_t : bool, optional
            Applies sigma_clipped_stats to image. Sigma is 3 by default.
    """
    
    cal_prod = cp(input_data)
    zero_cal_point = cp(zero_pol_data)
    
    #This could technically still be useful for other cases of data
    #data_strs = ['g191b2b', 'G191B2B']
    #G191_low_pol = funcs_utils.filter_data(zero_pol_data, data_strs, False)

    if(verbose_plot_zpol):
        means = funcs_polarimetry.q_n_u_single_plot_v1(zero_cal_point,
                                       plot_c='blue',
                                       only_means=True,
                                       pol_deg=True)
    
    #Redundant
    mean_q, mean_q_err = funcs_polarimetry.plot_q_u_stability(zero_cal_point, q_u_check='q', 
                                                          m_plot=verbose_plot_points, 
                                                          plot_verbose=verbose_plot_points)
    mean_u, mean_u_err = funcs_polarimetry.plot_q_u_stability(zero_cal_point, q_u_check='u', 
                                                          m_plot=verbose_plot_points, 
                                                          plot_verbose=verbose_plot_points)

    #This does that removal of points... how about that
    if(verbose_text):
        print("Removing unstable points")

    for points in list_index_unstable_data:
        zero_cal_point[0] = np.delete(zero_cal_point[0], points)
        zero_cal_point[1] = np.delete(zero_cal_point[1], points)

    if(verbose_plot_zpol):
        means = funcs_polarimetry.q_n_u_single_plot_v1(zero_cal_point,
                                       plot_c='blue',
                                       only_means=True,
                                       pol_deg=True)
        
    mean_q, mean_q_err = funcs_polarimetry.plot_q_u_stability(zero_cal_point, q_u_check='q', 
                                                          m_plot=verbose_plot_points, 
                                                          plot_verbose=verbose_plot_points)
    mean_u, mean_u_err = funcs_polarimetry.plot_q_u_stability(zero_cal_point, q_u_check='u', 
                                                          m_plot=verbose_plot_points, 
                                                          plot_verbose=verbose_plot_points)

    if(verbose_text):
        print("For all 0 pols. q n u instrumental points:")
        print("q inst:", mean_q, u"\u00B1",mean_q_err)
        print("u inst:", mean_u, u"\u00B1",mean_u_err)
        print("\n")    
    
    c = 1
    
    cal_targ = [[],
                []]
    
    #Sections defined arbitrarily
    targ_data_arr = [np.array(cal_prod)[:,:2],
                     np.array(cal_prod)[:,2:4],
                     np.array(cal_prod)[:,4:13], 
                     np.array(cal_prod)[:,13:22], 
                     np.array(cal_prod)[:,22:]]   
        
    targ_slice = ['0-2', 
                  '2-4',
                  '5-13',
                  '13-22',
                  '22-len(target_data)']
   
    for k in range(0, len(zero_cal_point[0])):
        mean_q, mean_q_err = funcs_polarimetry.plot_q_u_stability(np.array(zero_cal_point)[:,:c],
                                                                  q_u_check='q', 
                                                                  m_plot=False, 
                                                                  plot_verbose=False)
        mean_u, mean_u_err = funcs_polarimetry.plot_q_u_stability(np.array(zero_cal_point)[:,:c], 
                                                                  q_u_check='u', 
                                                                  m_plot=False)        

        q_cal = [np.mean(mean_q), np.mean(mean_q_err)]
        u_cal = [np.mean(mean_u), np.mean(mean_u_err)]

        cal_section = calib_data(targ_data_arr[k], (q_cal, u_cal))
        #cal_section = funcs_polarimetry.calib_data(targ_data_arr[k], (q_cal, u_cal))

        for di in range(0, len(cal_section[0])):
            cal_targ[0].append(cal_section[0][di])
            cal_targ[1].append(cal_section[1][di])
                  
        #compute this ting called the slice
        #the_slice = [list(x.keys())[0] for x in G191_low_pol[0:c]]
        #if(verbose_text):
        #    print("For 0 pols:", the_slice,the_slice[0])
        #    print("And targets:", targ_slice[k])
        #    print("q cal:", q_cal[0], u"\u00B1", q_cal[1])
        #    print("u cal:", u_cal[0], u"\u00B1", u_cal[1])
        #    print("\n")

        c =c + 1

    return(cal_targ)