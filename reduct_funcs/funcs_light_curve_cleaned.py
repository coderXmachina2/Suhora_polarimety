import astropy
import gzip
import glob
import numpy as np
import math
import xlrd
import matplotlib.pyplot as plt
import pickle

from astropy.time import Time,  TimeJD
from datetime import datetime
from copy import deepcopy as cp
from scipy import stats
from photutils.aperture import aperture_photometry
from astropy.visualization import simple_norm
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.visualization.mpl_normalize import ImageNormalize
#from astropy.visualization import SqrtStftch

from reduct_funcs import funcs_utils

#Functions that plot pol data with light curve?
#EECep_light_curve_based_pol    10 args
#EECep_stacked_based_pol        11 args
#EECep_stacked_based            10 args

#EECep_light_curve_loader_n_cleaner
#EECep_light_curve_split_A
#EECep_light_curve_split_B

#destroy outliers
#plot things

#Misc Funtions

def EECep_stacked_based(target_data, #PD, PDerror, MJD
                        light_curve_data, #colour = [0], MJD = [1]
                        MJD_cuts, #Tuple... Left and Right
                        lc=False,
                        verb_t=False,
                        true_pa=False, 
                        title_t=False):
    """
    I think we need to just redoo
    """
    
    #t_b = Time(B_JD , format='jd', scale='utc')
    #t_v = Time(V_JD , format='jd', scale='utc')
    #t_r = Time(R_JD , format='jd', scale='utc')
    #t_i = Time(I_JD , format='jd', scale='utc')
        
    color_arr = [B_arr, V_arr, R_arr, I_arr]
    mjd_arr_k=light_curve_data.append(target_data[2])
    
    
    #originals
    c_arr = ['blue',  'green',  'red', 'purple']
    c_labs = ['B', 'V', 'R', 'I', 'pol_d']
        
    index_cutt = []
    index_L_cutt = []
    
    for jk in range(0, len(mjd_arr_k)):
        idx = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts[0]).argmin() #Left_side
        idxL = np.abs(mjd_arr_k[jk].to_value('mjd', 'float')   - MJD_cuts[1]).argmin() #Right_Side
        index_cutt.append(int(idx))
        index_L_cutt.append(int(idxL))
    
    print("Index cutt:", index_cutt)
    print("Index L cutt:", index_L_cutt)
    
    #fig = plt.figure(figsize=(36, 12))
    #gs = fig.add_gridspec(3, hspace=0)
    #axs = gs.subplots(sharex=True, sharey=False)
    
    
    #if(true_pa):
    #    fig.suptitle('Polarization degree (PD), true Position angle (PA) and EE CEP light curve', 
    #                     fontsize=32)    
    #else:
    #    fig.suptitle('Polarization degree (PD), Position angle (PA) and EE CEP light curve', fontsize=32)
    
    