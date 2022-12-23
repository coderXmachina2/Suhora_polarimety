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
from kneed import KneeLocator

def apt_phot_global_bkg_sub(fits_data_1, search_offset, positions, apt_rad, plot_phot_tab, plot_sky_med): #, ann_in_rad, ann_out_rad   ):
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
         aperture radius
    plot_phot_tab : bool, optional
         plot image. True by default
    plot_sky_med : bool, optional
         Saves image to file.
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
    
def apt_phot_local_bkg_sub(fits_data_1, positions, search_array, img_offset,  apt_rad, ann_in_rad, ann_out_rad, plot = False):
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
    annulus_aperture = CircularAnnulus(positions, r_in= ann_in_rad, r_out=ann_out_rad)
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
        plt.legend(loc=(0.17, 0.05), facecolor='#458989', labelcolor='white', handles=handles, prop={'weight': 'bold', 'size': 11})
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

def solve_apt(combine_target, trial_radii,verbose):
    """
    A function that takes in a list of lists of combined targets and a list of trial_radii and fits for the knee point. Returns the median of the knee which is supposed to be a recommended radius for aperture photometry.
    
    Parameters
    ----------
    combine_target : list
        List of lists. The nested lists contain photometry fluxes.
    trial_radii : list
        List of trial radii. 
    verbose : bool, optional
        Prints some extrat things     
    """
    
    knees = []

    for x in range(0, len(combine_target)):
        #print(x)
        kn = KneeLocator(trial_radii, 
                         combine_target[x], 
                         curve='concave', 
                         direction='increasing')
        knees.append(kn.knee)
        if(verbose):
            plt.plot(combine_target)
            plt.grid()
            plt.show()
            print("Target:", x+1 ,kn.knee)

    if(verbose):
        print("Median aperture radii:", np.median(knees))
    
    return (np.median(knees))