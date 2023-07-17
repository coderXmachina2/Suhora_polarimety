import matplotlib.pyplot as plt
import numpy as np
import astropy
import importlib
import glob
import xlsxwriter
import re

from copy import deepcopy as cp

from reduct_funcs import funcs_utils
from reduct_funcs import funcs_calib_and_plot
from reduct_funcs import funcs_star_finder
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.aperture import aperture_photometry
from astropy.visualization import simple_norm
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from photutils import find_peaks
from photutils.detection import DAOStarFinder
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from kneed import KneeLocator

importlib.reload(funcs_star_finder)
importlib.reload(funcs_calib_and_plot)
importlib.reload(funcs_utils)

def apt_phot_global_bkg_sub(fits_data_1, 
                            search_offset, 
                            positions, 
                            apt_rad, 
                            plot_phot_tab=False,
                            plot_sky_med=False): #, ann_in_rad, ann_out_rad   ):
    """
    A function that does apeture photometry with global background subtraction. Returns a table of photometry fluxes with errors. Was verified to yield similarity to local background subtraction.

    Parameters
    ----------
        fits_data_1 : numpy ndarray
            Single numpy ndarray of image intensities. 2D image array.
        search_offset : list
            List of length 4 integers that determine a region to search within image
        positions : list
            List of integers. Array that scales both image data to act as zoom. 
        apt_rad : float
            Float defining inner aperture radius
        plot_phot_tab : bool, optional
            Plot image. False by default
        plot_sky_med : bool, optional
            Saves image to file. False by default
    """    
    search_this = fits_data_1[0].data[512-search_offset:512+search_offset,
                                       512-search_offset:512+search_offset]
    
    apertures = CircularAperture(positions, r=apt_rad) #Its just standard 4
    
    median     = np.median(fits_data_1[0].data)
    median_sky = np.median(fits_data_1[0].data[0:390, 200:400])
    bkg = (median_sky + median)/2

    errs = 0.1 * search_this

    if(plot_sky_med):
        print("Median data:", median, "Median sky", median_sky, "bkg:", bkg)
        plt.imshow(fits_data_1[0].data[0:390, 200:400])
        plt.colorbar()
        plt.title("Background region")
        plt.show()
    
    phot_table = aperture_photometry(search_this - bkg, apertures, error = errs ) 
   
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
        
    if (plot_phot_tab):
        print(phot_table)
        #print("\n Aperture are:", apertures.area, "Annulus Area:" , annulus_aperture.area, "\n")
    
    #I suppose getting the residual error is the same as getting the residual but with the error.
    
    #bkg_mean_err = phot_table['aperture_sum_err_1'] / annulus_aperture.area
    #bkg_sum_err = bkg_mean_err * apertures.area
    #final_sum_err = phot_table['aperture_sum_err_0'] - bkg_sum_err
    #phot_table['residual_aperture_sum_err'] = final_sum_err
    #phot_table['residual_aperture_sum_err'].info.format = '%.8g'
    
    #bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
    #bkg_sum  = bkg_mean * apertures.area
    #final_sum = phot_table['aperture_sum_0'] - bkg_sum
    #phot_table['residual_aperture_sum'] = final_sum
    #phot_table['residual_aperture_sum'].info.format = '%.8g'  # for consistent table output 
    
    #for col in phot_table.colnames:
    #    phot_table[col].info.format = '%.8g'  # for consistent table output
        
    #if (plot):
    #    print(phot_table)
    
    return phot_table
    
def apt_phot_local_bkg_sub(fits_data_1,
                           positions,
                           search_array,
                           img_offset,  
                           apt_rad,
                           ann_radii,
                           plot = False):
    """   
    A function that does apterture photometry with local background subtraction. Returns photometry table

    Parameters
    ----------
        fits_data_1 : numpy ndarray
            Single numpy ndarray of image intensities. 2D image array.
        positions : list
            List of lists. Should be a list of 2 lists for Suhora data polarimeter. The nested list contain x and y coordinates. 
        search_array :  list
            A list of length 4 containing integer limits that determine the area of interest.
        img_offset : int
            Integer that determines new bounds of sub image.
        apt_rad : float
            Aperture radius.
        plot : bool, optional
            Plots image. Prints photometry table.  
    """
    search_this = fits_data_1[512-img_offset:512+img_offset,
                                       512-img_offset:512+img_offset]
    
    apertures = CircularAperture(positions, r=apt_rad) #Its just standard 4
    annulus_aperture = CircularAnnulus(positions, r_in= ann_radii[0], r_out=ann_radii[1])
    errs = 0.1 * search_this
    
    if (plot):
        plt.imshow(search_this)#, cmap='Greys', origin='lower', norm=norm,interpolation='nearest')
        #plt.title(fits_data_1[0].header['OBJECT']+ " " +fits_data_1[0].header['FILTER'] + " reduced " + fits_data_1[0].header['TIME-OBS']+" DAOphot.")
        plt.axvline(x=len(search_this)/2, c='black', alpha=0.2)
        plt.axhline(y=len(search_this[0])/2, c='black', alpha=0.2)
        
        for z in range(0, len(positions)):
            plt.text(  positions[z][0],
                       positions[z][1], 
                       "Target " +str(z+1)+"\n" +str((round(positions[z][0], 2), round(positions[z][1],2))), 
                     fontsize=16, alpha=0.8)
        
        #should have the star tracker brackets
            plt.plot([img_offset+search_array[3],img_offset+search_array[3]], 
                     [img_offset-search_array[0],img_offset+search_array[1]], 
                     'k-', lw=1.75, alpha=0.4, linestyle = '--') #vertical
            plt.plot([img_offset-search_array[2],img_offset-search_array[2]], 
                     [img_offset-search_array[0],img_offset+search_array[1]], 
                     'k-', lw=1.75, alpha=0.4, linestyle ='--') #vertical 
            plt.plot([img_offset-search_array[2],img_offset+search_array[3]], 
                     [img_offset+search_array[1],img_offset+search_array[1]], 
                     'k-', lw=1.75, alpha=0.4, linestyle ='--') #horizontal
            plt.plot([img_offset+search_array[3],img_offset-search_array[2]], 
                     [img_offset-search_array[0],img_offset-search_array[0]],
                     'k-', lw=1.75, alpha=0.4, linestyle = '--') #horizontal 
               
        ap_patches = apertures.plot(color='white', lw=2,label='Photometry aperture')
        ann_patches = annulus_aperture.plot(color='red', lw=2,label='Background annulus')
        handles = (ap_patches[0], ann_patches[0])
        plt.legend(loc=(0.17, 0.05), 
                   facecolor='#458989', 
                   labelcolor='white', 
                   handles=handles, 
                   prop={'weight': 'bold', 'size': 11})
        plt.colorbar()
        plt.grid(lw='0.15')
        plt.show()
    
    apers = [apertures, annulus_aperture]
    phot_table = aperture_photometry(search_this, apers, error = errs )
    
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
        
    if (plot):
        print(phot_table)
        print("\n Aperture are:", apertures.area, "Annulus Area:" , annulus_aperture.area, "\n")
    
    #I suppose getting the residual error is the same as getting the residual but with the error.
    
    bkg_mean_err = phot_table['aperture_sum_err_1'] / annulus_aperture.area
    bkg_sum_err = bkg_mean_err * apertures.area
    final_sum_err = phot_table['aperture_sum_err_0'] - bkg_sum_err
    phot_table['residual_aperture_sum_err'] = final_sum_err
    phot_table['residual_aperture_sum_err'].info.format = '%.8g'
    
    bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
    bkg_sum  = bkg_mean * apertures.area
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    phot_table['residual_aperture_sum'] = final_sum
    phot_table['residual_aperture_sum'].info.format = '%.8g'  # for consistent table output 
    
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
        
    if (plot):
        print(phot_table)
    
    return phot_table

def process_obs_track(list_data_files,  
                      cal_files,
                      search_array,
                      radii_range,
                      track_ints,
                      ann_radii,
                      verbose_print=False):

    MJD = re.search('./files_sorted/(.*)/',  list_data_files[0]).group(1)
    
    if(verbose_print):
        print("Track inputs:",    
               "Search arrray:", search_array,
               "Start int:", track_ints[0], 
               "Stop int:", track_ints[1],
               "Flat field:")
        fig, ax = plt.subplots(figsize=(12, 12))
        plt.imshow(cal_files[2])
        plt.show()
        
        datadf = funcs_utils.print_list(MJD, list_data_files)
        print("The data frame")
        print(datadf)

    filename = 'ref' + astropy.io.fits.open(list_data_files[track_ints[0]])[0].header['OBJECT'] +"_" + astropy.io.fits.open(list_data_files[track_ints[0]])[0].header['FILTER']
    
    workbook = xlsxwriter.Workbook('./stats/'+MJD+"_"+ filename+str(track_ints[0])+"-" +str(track_ints[1])  +'.xlsx')
    worksheet = workbook.add_worksheet()
    worksheet.write('A1', 'int')
    worksheet.write('B1', 'time obs')
    worksheet.write('C1', 'filter')
    worksheet.write('D1', 'exptime')
    worksheet.write('E1', 'aperture radius')
    worksheet.write('F1', 'target 1 x center')
    worksheet.write('G1', 'target 1 y center')
    worksheet.write('H1', 'target 1 counts')
    worksheet.write('I1', 'target 1 error')
    worksheet.write('J1', 'target 2 x center')
    worksheet.write('K1', 'target 2 y center')
    worksheet.write('L1', 'target 2 counts')
    worksheet.write('M1', 'target 2 error')
    worksheet.write('N1', 'q')
    worksheet.write('O1', 'q error')
    worksheet.write('P1', 'u')
    worksheet.write('Q1', 'u error')
    worksheet.write('R1', 'PD')
    worksheet.write('S1', 'PD error')
    worksheet.write('T1', 'PA')
    worksheet.write('U1', 'PA error')
    
    row = 1
    column = 0
    search_off = 150
    sigma_search = 3

    #This is the processing loop
    for k in range(track_ints[0], track_ints[1]):
        if(verbose_print):
            print(k, list_data_files[k])
            
        #Key processing step make a reduced object
        #Discrepancy of structure with laptop. OMG!
        reduced_obj = funcs_calib_and_plot.plot_double_raw_v_reduced(list_data_files[k],
                                            (cal_files[0], cal_files[1], cal_files[2]),
                                            [512, 512, 512, 512], 
                                            sigma=True, 
                                            plot=True, 
                                            sv_img=False)
        
        trial_radii = np.linspace(radii_range[0], radii_range[1], num=60)
        
        target_a = []
        target_b = []
        
        for radii in trial_radii:            
            
            #Main Detector. This order is similar to laptop version
            DAO_positions = funcs_star_finder.dao_star_finder(fits_data_1=reduced_obj,
                                                              search_array=search_array,
                                                              siegma=sigma_search,   
                                                              trim = 1000,  #Hardcode this for now
                                                              img_offset=search_off, 
                                                              apt_rad=radii, 
                                                              ann_radii=(ann_radii[0],ann_radii[1]), 
                                                              plot=False, 
                                                              plot_tab=False)
            if (len(DAO_positions) != 2):
                #what do I return from this? 
                print("Second Detector Triggered")
                x_targ, y_targ, peak_targ = funcs_star_finder.source_peak_finder(reduced_obj, 
                                                                                 search_array,
                                                                                 3, 
                                                                                 1000, 
                                                                                 plot_peaks=False, 
                                                                                 verbose=False)

                apt_pos = funcs_star_finder.plot_spotted(reduced_obj, 
                                                         search_array, 
                                                         x_targ, 
                                                         y_targ, 
                                                         peak_targ, 
                                                         search_off, 
                                                         plot_im=False)
                
                DAO_positions = funcs_star_finder.peak_to_DAO(apt_pos)
            
            phot_tab = funcs_apt_phot.apt_phot_local_bkg_sub(reduced_obj,
                                                             DAO_positions,
                                                             search_bracket, 
                                                             search_off,  
                                                             radii, 
                                                             (12, 22), 
                                                             False)
            
            target_a.append(phot_tab['residual_aperture_sum'][0])
            target_b.append(phot_tab['residual_aperture_sum'][1])
            
        combine_target = [target_a, target_b]    
        
        good_radii=funcs_apt_phot.solve_apt(combine_target, #
                                            trial_radii, 
                                            verbose=True) #
        
    workbook.close()
    
def solve_apt(combine_target, 
              trial_radii,
              verbose=False):
    """
    A function that takes in a list of lists of combined targets and a list of trial_radii and fits for the knee point. Returns the median of the knee which is supposed to be a recommended radius for aperture photometry.
    
    Parameters
    ----------
        combine_target : list
            List of lists. The nested lists contain photometry fluxes (counts).
        trial_radii : list
            List of trial radii. 
        verbose : bool, optional
            Prints some extrat things     
    """
    
    knees = []
    c = 0
    #print("Expectation len of 2 with 60 elements inside:", len(combine_target), combine_target)
    
    for x in range(0, len(combine_target)):
        #For both targets 
        #plots the knee locator
        kn = KneeLocator(trial_radii, 
                         combine_target[x], 
                         curve='concave', 
                         direction='increasing')
        knees.append(kn.knee)
        if(verbose):
            #what is this x and y fix the plot
            fig, ax1 = plt.subplots(figsize=(9, 9))
            plt.plot(  trial_radii,  combine_target[x]   )

            plt.title("Target " + str(c) + " Photometry vs Aperture Size")
            plt.text(np.median(knees), 
                     np.min(combine_target[x])+(np.max(combine_target[x])-np.min(combine_target[x]))/2, 
                     "Aperture Photometry Radius Solution: \n"+ str( np.median(knees) )  )          
            plt.axvline(x= np.median(knees), linestyle='--' )
            
            plt.ylabel("Photometric Counts")
            plt.xlabel("Aperture Radius")
            plt.grid()
            plt.show()
            print("Target:", x+1 ,kn.knee)
        c+=1
        
    if(verbose):
        print("Median aperture radii:", np.median(knees))
    
    return (np.median(knees))