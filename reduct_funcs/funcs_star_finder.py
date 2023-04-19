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

def source_peak_finder(fits_data_1, 
                       search_array, 
                       siegma, 
                       trim, 
                       plot_peaks=False, 
                       verbose=False):
    """
    A function that finds location and values of peaks in the image bounded within offsets given in array (that define a region of interest shaped in a rectangle/square). Involves some further trhesholding/ filtering routines. Returns list of x coordinates, list of y coordinates, and list of the peak in image. This function is succeeded with plot_spotted for maximum results.
    
    Parameters
    ----------
    fits_data_1 : numpy ndarray
        Single numpy ndarray of image intensities. 2D image array.
    search_array :  list
        A list of length 4 containing integer limits that determine the area of interest.
    siegma : int
        Integer that determines sigma cutoff for sigma_clipped_stats image search function.
    trim : int
         Integer in the hundreds to thousand. Defines a lower limit used to calculate a rejection threshold. 
    plot_peaks : bool, optional
         Plot peaks for verification
    verbose : bool, optional
         Prints results Plots peaks for verification
    """
    x_peak = []
    y_peak = []
    peak_val =[]

    x_interest = []
    y_interest = []
    peak_interest = []
    index_x = []

    mean, median, std = sigma_clipped_stats(fits_data_1, sigma=siegma) #searches the whole image
    threshold = median + (10. * std)
    tbl = find_peaks(fits_data_1, threshold, box_size=40)
    tbl['peak_value'].info.format = '%.3g'  # for consistent table output
    
    #Gets locations
    for things in tbl:
        x_peak.append(things['x_peak'])
        y_peak.append(things['y_peak'])
        peak_val.append(things['peak_value'])

    #looks inside this zone
    for i in range(0, len(peak_val)): 
        """
        #TODO replace search within circle
        #Calculate a distance r
        #save point if distnace to point (d = sqrt( (x-x)^2 + (y-y)^2) is less than some input search offset
        """
        if(y_peak[i] > 512-search_array[0] and 
           y_peak[i] < 512+search_array[1] and 
           x_peak[i] > 512-search_array[2] and 
           x_peak[i] < 512+search_array[3]):
            x_interest.append(x_peak[i])
            y_interest.append(y_peak[i])
            peak_interest.append(peak_val[i])
            index_x.append(i)
    
    x_targ = []
    y_targ = []
    peak_targ = []
    sub_sample = []
    
    for peaks in peak_val:
        if peaks < trim:
            sub_sample.append(peaks)
    
    #Calculates threshold from a sub sample
    threshold = np.mean(sub_sample)+5*np.std(sub_sample)
    
    for z in range(0, len(peak_interest)):
        if(peak_interest[z] > threshold):
            x_targ.append(x_interest[z])
            y_targ.append(y_interest[z])
            peak_targ.append(peak_interest[z])
    
    if(plot_peaks):
        xaxis = np.linspace(0, len(peak_val[:]), len(peak_val[:]))
        plt.scatter(xaxis, peak_val[:])
        plt.title("Image flux peaks")
        plt.grid()
        
        plt.text(0, threshold, 'Threshold ', fontsize=16, color = 'red')
        plt.text(0, trim,'Trim ',fontsize=16,color = 'red')
        
        plt.axhline(y=threshold, color = 'red', linestyle='--', lw=2.5) 
        plt.axhline(y=trim,color = 'red',linestyle='--',lw=2.5)
        for i in range(0, len(index_x)):
            plt.plot([index_x[i],index_x[i]], [peak_interest[i]-300, peak_interest[i]+300], 'k-', lw=2.5) 
            plt.plot([index_x[i]-2.5,index_x[i]+2.5], [peak_interest[i], peak_interest[i]] , 'k-', lw=2.5) 
        plt.show()

    if(verbose):
        print(len(tbl), 
              "peaks detected from image of size", 
              "HEADER NAXIS1 x NAXIS2",
              #fits_data_1[0].header['NAXIS1'], "x", fits_data_1[0].header['NAXIS2'],
              "with sigma:", siegma, "and threshold:", threshold, "threshold is nan and there is error" if math.isnan(threshold) else 
              "valid" )

        print("Targets within region of interest: ", len(x_targ))

    return(x_targ, y_targ, peak_targ)


def plot_spotted(fits_data_1, 
                 search_array, 
                 x_targ, 
                 y_targ, 
                 peak_targ, 
                 img_offset, 
                 plot_im=False):
    """
    A function that takes in 2D image data and plots them along with search offsets. Creats a list of tuples with x and y coordinates packaged together. Takes in list of x and y coordinates and peak fluxes to put them on the image. Returns a list of tuples which are apeture positions. apeture positions x and y coordinates are calculated based on the search offset. So what if the image is not 512? We may face problems moving forward.
    
    Parameters
    ----------
    fits_data_1 : numpy ndarray
        Single numpy ndarray of image intensities. 2D image array.
    search_array :  list
        A list of length 4 containing integer limits that determine the area of interest.
    x_targ : list
         A list of x positions.
    y_targ : list
         A list of y positions.
    peak_targ : list
         A list of  peak intensities at the x and y coordinates.
    img_offset : int
        Integer that determines new bounds of sub image. 
    plot_im :  bool, optional
         Plot image.
    """    
    search_offset = img_offset
    
    apt_positions = []
    for xint in range(0,len(x_targ)):
        apt_positions.append((x_targ[xint]-(512-search_offset), y_targ[xint]-(512-search_offset)))
    
    if(plot_im):
        plt.imshow(fits_data_1[512-img_offset:512+img_offset,
                                       512-img_offset:512+img_offset])
        plt.axvline(x=img_offset, c='black', alpha=0.2)
        plt.axhline(y=img_offset, c='black', alpha=0.2)
        plt.grid(lw='0.15')
        plt.colorbar()
        #plt.title(fits_data_1[0].header['OBJECT']+ " " +fits_data_1[0].header['FILTER'] + " reduced " + fits_data_1[0].header['TIME-OBS'])

        apt_positions = []
        plt_alpha= 0.2
        
        for xint in range(0,len(x_targ)):
            plt.text(x_targ[xint]-(512-search_offset),
                     y_targ[xint]-(512-search_offset),
                     'Target '+str((x_targ[xint]-(512-search_offset),
                                    y_targ[xint]-(512-search_offset)))+
                     "\nPeak_val:"+str(round(peak_targ[xint],2)),fontsize=16)
            
            apt_positions.append((x_targ[xint]-(512-search_offset), y_targ[xint]-(512-search_offset)))

            plt.plot([search_offset+search_array[3],search_offset+search_array[3]], 
                     [search_offset-search_array[0],search_offset+search_array[1]], 
                     'k-', lw=1.75, alpha=0.4, linestyle = '--') #vertical
            plt.plot([search_offset-search_array[2],search_offset-search_array[2]], 
                     [search_offset-search_array[0],search_offset+search_array[1]], 
                     'k-', lw=1.75, alpha=0.4, linestyle ='--') #vertical 
            plt.plot([search_offset-search_array[2],search_offset+search_array[3]], 
                     [search_offset+search_array[1],search_offset+search_array[1]], 
                     'k-', lw=1.75, alpha=0.4, linestyle ='--') #horizontal
            plt.plot([search_offset+search_array[3],search_offset-search_array[2]],
                     [search_offset-search_array[0],search_offset-search_array[0]],
                     'k-', lw=1.75, alpha=0.4, linestyle = '--') #horizontal
            
            plt.plot([x_targ[xint]-(512-search_offset),
                      x_targ[xint]-(512-search_offset)], 
                     [y_targ[xint]-(512-search_offset)-2,
                      y_targ[xint]-(512-search_offset)+2], 
                     'k-', lw=2) #vertical
            plt.plot([x_targ[xint]-(512-search_offset),
                      x_targ[xint]-(512-search_offset)], 
                     [y_targ[xint]-(512-search_offset)-2, 
                      y_targ[xint]-(512-search_offset)+2], 
                     'k-', lw=2) #vertical
            plt.plot([x_targ[xint]-(512-search_offset)-2, 
                      x_targ[xint]-(512-search_offset)+2], 
                     [y_targ[xint]-(512-search_offset),
                      y_targ[xint]-(512-search_offset)], 
                     'k-', lw=2) #horizontal line
        plt.show()
    
    return apt_positions

def peak_to_DAO(apt_pos):
    """
    A function to that converts aperture positions returned by source_peak_finder to a format acceptable by DAOStarFinder. Returns a list
    
    Parameters
    ----------
    apt_pos : list
        List of tuples containing x and y locations relative to a predefined search offset.
    """    
    
    DAO_array = []
    for pos in apt_pos:
        pos_array = [   ]

        pos_array.append(pos[0])
        pos_array.append(pos[1])
        DAO_array.append(pos_array)
        
    return(DAO_array)
    
def dao_star_finder(fits_data_1, 
                    search_array, 
                    siegma, 
                    trim,
                    img_offset, 
                    apt_rad,
                    ann_in_rad, 
                    ann_out_rad, 
                    plot=False,
                    plot_tab=False):
    """
    A function that finds location and values of peaks in the image bounded within offsets given in array (that define a region of interest shaped in a rectangle/square). Involves some further trhesholding/ filtering routines. Does a little bit of preliminary aperture photometry but does not return the results. Returns list of x coordinates, list of y coordinates, and list of the peak in image. This function is used in pipeline counts extraction cell. Returns positions of target of interest.
    
    Parameters
    ----------
    fits_data_1 : numpy ndarray
        Single numpy ndarray of image intensities. 2D image array.
    search_array :  list
        A list of length 4 containing integer limits that determine the area of interest.
    siegma : int
        Integer that determines sigma cutoff for sigma_clipped_stats image search function   
    trim : int
         Integer in the hundreds to thousand. Defines a lower limit used to calculate a rejection threshold. 
    img_offset : int
         Integer that determines new bounds of sub image.
    apt_rad : float
         aperture radius
    ann_in_rad : int
         Integer that determines inner radius of the circular annulus 
    ann_out_rad :  int
         Integer that determines outer radius of the circular annulus 
    plot :  bool, optional
         Plot aperture photometry image
    plot_tab :  bool, optional
         Plot aperture photometry table,
    """    
    
    search_this = fits_data_1[512-img_offset:512+img_offset,512-img_offset:512+img_offset]

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
        print("second_thresh", trim)
    
    for peaks in peak_val: 
        if peaks < trim:
            sub_sample.append(peaks)
                
    threshold = np.mean(sub_sample)+5*np.std(sub_sample)
    
    for i in range(0, len(measure_flux)): #looks inside zone
        """
        #TODO replace with a circle
        #Calculate a distance r
        #if distnace to point (d = sqrt( (x-x)^2 + (y-y)^2) is less than some input d#search offset is the center
        """
        if(y_peak[i] > img_offset-search_array[0] and 
           y_peak[i] < img_offset+search_array[1] and 
           x_peak[i] > img_offset-search_array[2] and 
           x_peak[i] < img_offset+search_array[3] and peak_val[i] > threshold):
            x_targ.append(x_peak[i])
            y_targ.append(y_peak[i])
            flux_targ.append(measure_flux[i])
            peak_targ.append(peak_val[i])
    
    positions = np.transpose((x_targ, y_targ)) #positions
    apertures = CircularAperture(positions, r=apt_rad) #Its just standard 4
    annulus_aperture = CircularAnnulus(positions, r_in= ann_in_rad, r_out=ann_out_rad)
    
    if (plot):
        plt.imshow(search_this)#, cmap='Greys', origin='lower', norm=norm,interpolation='nearest')
        #plt.title(fits_data_1[0].header['OBJECT']+ " " +fits_data_1[0].header['FILTER'] + " reduced " + fits_data_1[0].header['TIME-OBS']+" DAOphot")
        plt.axvline(x=len(search_this)/2, c='black', alpha=0.2)
        plt.axhline(y=len(search_this[0])/2, c='black', alpha=0.2)

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
                 'k-', lw=1.75, alpha=0.4, linestyle = '--') #horizontal lines

        for xint in range(0,len(x_targ)):
            plt.text(x_targ[xint], 
                     y_targ[xint], 
                     "Target "+ str(xint+1) + " " + str((round(x_targ[xint], 2), round(y_targ[xint],2))) + "\nFlux: " +
                     str(round(flux_targ[xint], 4)) +
                     "\nPeak_val: " + str(round(peak_targ[xint], 4)), 
                     fontsize=16, alpha=0.4)

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
    Experimental. I attempted it a few times in 2021. It works but not refined as of 23-12-2019
    
    Parameters
    ----------
    fits_data_1 : numpy ndarray
        Single numpy ndarray of image intensities. 2D image array.
    search_array :  
    siegma : 
    second_thresh : 
    search_offset : 
    apt_rad : 
    ann_in_rad : 
    ann_out_rad : 
    plot : 
    plot_tab : 
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
        #plt.title(fits_data_1[0].header['OBJECT']+ " " +fits_data_1[0].header['FILTER'] + " reduced " + fits_data_1[0].header['TIME-OBS']+" DAOphot")
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
                     fontsize=16, alpha=0.4)

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
