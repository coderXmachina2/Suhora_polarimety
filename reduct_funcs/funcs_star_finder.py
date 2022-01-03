import matplotlib.pyplot as plt
import numpy as np
import astropy
import math
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

def source_peak_finder(fits_data_1, siegma, search_array, trim, plot_peaks, verbose):
    """
    Function that finds coordinates of stars.
    
    Input should be sigma
    """
    x_peak = []
    y_peak = []
    peak_val =[]

    x_interest = []
    y_interest = []
    peak_interest = []
    index_x = []

    mean, median, std = sigma_clipped_stats(fits_data_1[0].data, sigma=siegma) #searches the whole image
    threshold = median + (10. * std)
    tbl = find_peaks(fits_data_1[0].data, threshold, box_size=40)
    tbl['peak_value'].info.format = '%.3g'  # for consistent table output
    
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
        if peaks < trim:
            sub_sample.append(peaks) #But it cannot identify it out of this sample...
            
    threshold = np.mean(sub_sample)+5*np.std(sub_sample) #This becomes nan and shit is fucked
    
    for z in range(0, len(peak_interest)):
        if(peak_interest[z] > threshold):
            x_targ.append(x_interest[z])
            y_targ.append(y_interest[z])
            peak_targ.append(peak_interest[z])
    #Output
    if(plot_peaks):
        plt.plot(peak_val[:])
        plt.title("Flux peaks")
        plt.grid()
        
        plt.text(0, threshold,
                 'Region of interest threshold ', 
                 fontsize=16, alpha=10, color = 'red')
        plt.axhline(y=threshold, 
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
              "with sigma:", siegma, "and threshold:", threshold, "threshold is nan and shit is fucked" if math.isnan(threshold) else "valid" )

        print("Targets within region of interest: ", len(x_targ))

    return(x_targ, y_targ, peak_targ)


def plot_spotted(fits_data_1, img_offset, search_array ,x_targ, y_targ, peak_targ, plot_im, fade_plot):
    """
    Function that takes in two plots and plots them. It just plots
    """    
    search_offset = img_offset
    
    apt_positions = []
    for xint in range(0,len(x_targ)):
        apt_positions.append((x_targ[xint]-(512-search_offset), y_targ[xint]-(512-search_offset)))
    
    if(plot_im):
        plt.imshow(fits_data_1[0].data[512-img_offset:512+img_offset,
                                       512-img_offset:512+img_offset])
        plt.axvline(x=img_offset, c='black', alpha=0.2)
        plt.axhline(y=img_offset, c='black', alpha=0.2)
        plt.grid(lw='0.15')
        plt.colorbar()
        plt.title(fits_data_1[0].header['OBJECT']+ " " +fits_data_1[0].header['FILTER'] + " reduced " + fits_data_1[0].header['TIME-OBS'])

        apt_positions = []
        plt_alpha= 0.2
        
        for xint in range(0,len(x_targ)):
            plt.text(x_targ[xint]-(512-search_offset),
                     y_targ[xint]-(512-search_offset),
                     'Target '+str((x_targ[xint]-(512-search_offset),
                                    y_targ[xint]-(512-search_offset)))+"\nPeak_val:"+str(round(peak_targ[xint],2)) , 
                     fontsize=16, alpha=10)
            
            apt_positions.append((x_targ[xint]-(512-search_offset), y_targ[xint]-(512-search_offset)))

            #This thing is freaking magical!
            #I need to implement a bracketing for my first filter
            plt.plot([search_offset+search_array[3],search_offset+search_array[3]], [search_offset-search_array[0],search_offset+search_array[1]], 'k-', lw=1.75, alpha=0.4, linestyle = '--') #vertical
            plt.plot([search_offset-search_array[2],search_offset-search_array[2]], [search_offset-search_array[0],search_offset+search_array[1]], 'k-', lw=1.75, alpha=0.4, linestyle ='--') #vertical 
            plt.plot([search_offset-search_array[2],search_offset+search_array[3]], [search_offset+search_array[1],search_offset+search_array[1]], 'k-', lw=1.75, alpha=0.4, linestyle ='--') #horizontal
            plt.plot([search_offset+search_array[3],search_offset-search_array[2]], [search_offset-search_array[0],search_offset-search_array[0]], 'k-', lw=1.75, alpha=0.4, linestyle = '--') #horizontal lines
            
            plt.plot([x_targ[xint]-(512-search_offset),
                      x_targ[xint]-(512-search_offset)], 
                     [y_targ[xint]-(512-search_offset)-2,
                      y_targ[xint]-(512-search_offset)+2], 'k-', lw=2) #vertical line
            plt.plot([x_targ[xint]-(512-search_offset),
                      x_targ[xint]-(512-search_offset)], 
                     [y_targ[xint]-(512-search_offset)-2, y_targ[xint]-(512-search_offset)+2], 'k-', lw=2) #vertical line
            plt.plot([x_targ[xint]-(512-search_offset)-2, x_targ[xint]-(512-search_offset)+2], 
                     [y_targ[xint]-(512-search_offset),y_targ[xint]-(512-search_offset) ], 'k-', lw=2) #horizontal line
        plt.show()
    
    return apt_positions

def peak_to_DAO(apt_pos):
    DAO_array = []
    for pos in apt_pos:
        pos_array = [   ]

        pos_array.append(pos[0])
        pos_array.append(pos[1])
        DAO_array.append(pos_array)
        
    return(DAO_array)
    
def dao_star_finder(fits_data_1, search_array, siegma, second_thresh ,search_offset, apt_rad, ann_in_rad, ann_out_rad, plot, plot_tab):
    
    search_this = fits_data_1[0].data[512-search_offset:512+search_offset,
                                   512-search_offset:512+search_offset]
    
    #DAO star finder. Takes in argument sigma 
    mean, median, std = sigma_clipped_stats(search_this, sigma=siegma)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std) #what happens if we change fwhm
    sources = daofind(search_this - median)#, mask=mask)  #sources is a table just like tbl
    
    if plot_tab:
        print("DAO found", len(sources), "objects discovered")
        print(sources)
    
    for col in sources.colnames:  
        sources[col].info.format = '%.8g'  #for consistent table output
        
    x_peak = []
    y_peak = []
    measure_flux = []
    peak_val = []
    
    for things in sources:
        x_peak.append(things['xcentroid'])    #xcentroid #ycentroid
        y_peak.append(things['ycentroid'])
        measure_flux.append(things['flux'])
        peak_val.append(things['peak'])
        
    #second round thresholding (second filter)    
    x_targ = []
    y_targ = []
    flux_targ = []
    peak_targ = []
    sub_sample = []
    
    if plot_tab:
        print("second_thresh", second_thresh)
    
    for peaks in peak_val: 
        if peaks < second_thresh: #was 2 thousand 
            sub_sample.append(peaks)
                
    threshold = np.mean(sub_sample)+5*np.std(sub_sample)
    
    for i in range(0, len(measure_flux)): #looks inside zone
        """
        #TODO replace with a circle
        #Calculate a distance r
        #if distnace to point (d = sqrt( (x-x)^2 + (y-y)^2) is less than some input d#search offset is the center
        """
        if(y_peak[i] > search_offset-search_array[0] and 
           y_peak[i] < search_offset+search_array[1] and 
           x_peak[i] > search_offset-search_array[2] and 
           x_peak[i] < search_offset+search_array[3] and peak_val[i] > threshold):
            x_targ.append(x_peak[i])
            y_targ.append(y_peak[i])
            flux_targ.append(measure_flux[i])
            peak_targ.append(peak_val[i])
    
    positions = np.transpose((x_targ, y_targ)) #positions
    apertures = CircularAperture(positions, r=apt_rad) #Its just standard 4
    annulus_aperture = CircularAnnulus(positions, r_in= ann_in_rad, r_out=ann_out_rad)
    
    if (plot):
        plt.imshow(search_this)#, cmap='Greys', origin='lower', norm=norm,interpolation='nearest')
        plt.title(fits_data_1[0].header['OBJECT']+ " " +fits_data_1[0].header['FILTER'] + " reduced " + fits_data_1[0].header['TIME-OBS']+" DAOphot")
        plt.axvline(x=len(search_this)/2, c='black', alpha=0.2)
        plt.axhline(y=len(search_this[0])/2, c='black', alpha=0.2)

        #This thing is freaking magical!
        plt.plot([search_offset+search_array[3],search_offset+search_array[3]], [search_offset-search_array[0],search_offset+search_array[1]], 'k-', lw=1.75, alpha=0.4, linestyle = '--') #vertical
        plt.plot([search_offset-search_array[2],search_offset-search_array[2]], [search_offset-search_array[0],search_offset+search_array[1]], 'k-', lw=1.75, alpha=0.4, linestyle ='--') #vertical 
        plt.plot([search_offset-search_array[2],search_offset+search_array[3]], [search_offset+search_array[1],search_offset+search_array[1]], 'k-', lw=1.75, alpha=0.4, linestyle ='--') #horizontal
        plt.plot([search_offset+search_array[3],search_offset-search_array[2]], [search_offset-search_array[0],search_offset-search_array[0]], 'k-', lw=1.75, alpha=0.4, linestyle = '--') #horizontal lines

        for xint in range(0,len(x_targ)):
            plt.text(x_targ[xint], 
                     y_targ[xint], 
                     "Target "+ str(xint+1) + " " + str((round(x_targ[xint], 2), round(y_targ[xint],2))) + "\nFlux: " +
                     str(round(flux_targ[xint], 4)) +
                     "\nPeak_val: " + str(round(peak_targ[xint], 4)), 
                     fontsize=16, alpha=10)

        ap_patches = apertures.plot(color='white', lw=2,
                                   label='Photometry aperture')
        ann_patches = annulus_aperture.plot(color='red', lw=2,
                                            label='Background annulus')
        handles = (ap_patches[0], ann_patches[0])

        plt.legend(loc=(0.17, 0.05), facecolor='#458989', labelcolor='white',
                   handles=handles, prop={'weight': 'bold', 'size': 11})
        plt.colorbar()
        plt.grid(lw='0.15')
        plt.show()
    
    return positions

def dao_star_finder_HD104860(fits_data_1, search_array, siegma, second_thresh ,search_offset, apt_rad, ann_in_rad, ann_out_rad, plot, plot_tab):
    """
    #Special function for 104860
    """
    
    search_this = fits_data_1[0].data[512-search_offset:512+search_offset,
                                   512-search_offset:512+search_offset]
    
    #DAO star finder. Takes in argument sigma 
    mean, median, std = sigma_clipped_stats(search_this, sigma=siegma)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std) #what happens if we change fwhm
    sources = daofind(search_this - median)#, mask=mask)  #sources is a table just like tbl
    
    if plot_tab:
        print("DAO found", len(sources), "objects discovered")
        print(sources)
    
    for col in sources.colnames:  
        sources[col].info.format = '%.8g'  #for consistent table output
        
    x_peak = []
    y_peak = []
    measure_flux = []
    peak_val = []
    
    for things in sources:
        x_peak.append(things['xcentroid'])    #xcentroid #ycentroid
        y_peak.append(things['ycentroid'])
        measure_flux.append(things['flux'])
        peak_val.append(things['peak'])
        
    #second round thresholding (second filter)    
    x_targ = []
    y_targ = []
    flux_targ = []
    peak_targ = []
    sub_sample = []
    
    if plot_tab:
        print("second_thresh", second_thresh)
    
    for peaks in peak_val: 
        if peaks < second_thresh: #was 2 thousand 
            sub_sample.append(peaks)
                
    threshold = np.mean(sub_sample)+5*np.std(sub_sample)
    
    for i in range(0, len(measure_flux)): #looks inside zone
        """
        #TODO replace with a circle
        #Calculate a distance r
        #if distnace to point (d = sqrt( (x-x)^2 + (y-y)^2) is less than some input d#search offset is the center
        """
        if(y_peak[i] > search_offset-search_array[0] and 
           y_peak[i] < search_offset+search_array[1] and 
           x_peak[i] > search_offset-search_array[2] and 
           x_peak[i] < search_offset+search_array[3] and peak_val[i] > threshold):
            x_targ.append(x_peak[i])
            y_targ.append(y_peak[i])
            flux_targ.append(measure_flux[i])
            peak_targ.append(peak_val[i])
    
    positions = np.transpose((x_targ, y_targ)) #positions
    apertures = CircularAperture(positions, r=apt_rad) #Its just standard 4
    annulus_aperture = CircularAnnulus(positions, r_in= ann_in_rad, r_out=ann_out_rad)
    
    if (plot):
        plt.imshow(search_this)#, cmap='Greys', origin='lower', norm=norm,interpolation='nearest')
        plt.title(fits_data_1[0].header['OBJECT']+ " " +fits_data_1[0].header['FILTER'] + " reduced " + fits_data_1[0].header['TIME-OBS']+" DAOphot")
        plt.axvline(x=len(search_this)/2, c='black', alpha=0.2)
        plt.axhline(y=len(search_this[0])/2, c='black', alpha=0.2)

        #This thing is freaking magical!
        plt.plot([search_offset+search_array[3],search_offset+search_array[3]], [search_offset-search_array[0],search_offset+search_array[1]], 'k-', lw=1.75, alpha=0.4, linestyle = '--') #vertical
        plt.plot([search_offset-search_array[2],search_offset-search_array[2]], [search_offset-search_array[0],search_offset+search_array[1]], 'k-', lw=1.75, alpha=0.4, linestyle ='--') #vertical 
        plt.plot([search_offset-search_array[2],search_offset+search_array[3]], [search_offset+search_array[1],search_offset+search_array[1]], 'k-', lw=1.75, alpha=0.4, linestyle ='--') #horizontal
        plt.plot([search_offset+search_array[3],search_offset-search_array[2]], [search_offset-search_array[0],search_offset-search_array[0]], 'k-', lw=1.75, alpha=0.4, linestyle = '--') #horizontal lines

        for xint in range(0,len(x_targ)):
            plt.text(x_targ[xint], 
                     y_targ[xint], 
                     "Target "+ str(xint+1) + " " + str((round(x_targ[xint], 2), round(y_targ[xint],2))) + "\nFlux: " +
                     str(round(flux_targ[xint], 4)) +
                     "\nPeak_val: " + str(round(peak_targ[xint], 4)), 
                     fontsize=16, alpha=10)

        ap_patches = apertures.plot(color='white', lw=2,
                                   label='Photometry aperture')
        ann_patches = annulus_aperture.plot(color='red', lw=2,
                                            label='Background annulus')
        handles = (ap_patches[0], ann_patches[0])

        plt.legend(loc=(0.17, 0.05), facecolor='#458989', labelcolor='white',
                   handles=handles, prop={'weight': 'bold', 'size': 11})
        plt.colorbar()
        plt.grid(lw='0.15')
        plt.show()
    
    return positions
