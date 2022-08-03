import gzip
import shutil
import glob
import astropy

from reduct_funcs import funcs_calib_and_plot
from reduct_funcs import funcs_star_finder
from reduct_funcs import funcs_apt_phot
from reduct_funcs import funcs_utils
from reduct_funcs import funcs_polarimetry

def make_dir(directory_path ):
    """
    #used for extracting and unloading data
    """
    x = glob.glob(directory_path+'/*')
    print("Inside dir:", directory_path, len(x), "Things detected")
    print("Print without .gz 3:", x[0][:-3])
    for archives in x:
        with gzip.open(archives, 'rb') as f_in:
            with open( archives[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

def print_list(MJD, f_list):
    print(MJD)
    for i in range(0, len(f_list)):
        print(i, 
        astropy.io.fits.open(f_list[i])[0].header['TIME-OBS'],
        astropy.io.fits.open(f_list[i])[0].header['FILTER'],
        astropy.io.fits.open(f_list[i])[0].header['EXPTIME'],
        astropy.io.fits.open(f_list[i])[0].header['OBJECT'])
        
        
def home_light(reduced_obj, search_bracket, search_off, search_thresh):
    #But stormfront will always be ma boo.
    #Lets make this into a search function
    #scanner test. Nop typically in a loop.
    #The next question is... does it apt phot. It does but strangeness occurs.
    #Search function. The new protocols don't undergo the search function. Proceed as usual.
    
    #search_off = 150
    x_targ, y_targ, peak_targ  = funcs_star_finder.source_peak_finder(reduced_obj,3,search_bracket, search_thresh , True, True) #what do I return from this?
    apt_pos = funcs_star_finder.plot_spotted(reduced_obj, search_off, search_bracket,x_targ, y_targ, peak_targ, True, False)

    apt_positions = []
    for h in range(0 , len(x_targ)):
        apt_positions.append((x_targ[h], y_targ[h])) #This must undergo some translation    

    print("Aperture positions are:", apt_positions, "from the", reduced_obj[0].header['NAXIS1'], "x", reduced_obj[0].header['NAXIS2'], "image")
    print("Aperture positions are:", apt_pos, "from the", str(search_off), "x", str(search_off) ,"image")

    