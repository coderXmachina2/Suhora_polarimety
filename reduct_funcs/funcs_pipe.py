import gzip
import shutil
import glob
import astropy

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