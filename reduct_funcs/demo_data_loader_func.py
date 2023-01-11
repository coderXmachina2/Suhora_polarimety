#Import dependencies
import re #
import pandas as pd #
import xlrd #

import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime, timedelta
from astropy.time import Time

def load_pol_data_demo_exr(in_str, 
                           verbose_file=False, 
                           verbose_arrays = False, 
                           verbose_midpoint=False,
                           verbose_val_checks=False):
    """
    #Start demonstrating concepts here
    """
    
    path_s = re.search('(.*)master_', in_str)
    f_name_s = re.search('master_(.*)', in_str)
    
    f_name = 'master_' + f_name_s.group(1)
    
    result_date = re.search('master_(.*)_', f_name[:19])
    result_obj = re.search('_(.*)_P', f_name[15:])
    
    df1 = pd.read_excel(in_str,engine='openpyxl') #1 Since we use this 
    
    if(verbose_file):
        print("File path:", path_s)
        print("File name:", f_name)
        print("Date:", result_date.group(1))
        print("obj:", result_obj.group(1))
       
    t_obs = np.array(df1['time obs'])
    q_arr = np.array(df1['q'])
    q_err_arr = np.array(df1['q error'])
    u_arr = np.array(df1['u'])
    u_err_arr = np.array(df1['u error'])
    
    #Introduce the problem... we need time. the time data is in the excel.
    #Suppose we want to return the time value but we only need one timestamp for the entire observation 
    
    #If this was not the case it could be as easy as this:
    #ret_dict = {result_date.group(1)+"_"+ result_obj.group(1): (q_arr,  q_err_arr, u_arr, u_err_arr, t_obs)} 
    #Or it could be as easy as this:
    #ret_dict = {result_date.group(1)+"_"+ result_obj.group(1): [(q_arr,  q_err_arr, u_arr, u_err_arr), t_obs]}
    
    #In this case I just want it returned separately from the dictionary result in a tuple or array.
    #How would I do that?
    
    #We could return the midpoint of the t_obs
    #Ok, how do we calculate midpoint of t_obs?
    #Use date time. We can start by collecting an array of date times. First import date time
    
    date_list = [datetime.strptime(ts, "%Y-%m-%d_%H:%M:%S.%f") for ts in [result_date.group(1)+"_"+x for x in t_obs]]   
    
    #we can verify with a print prints
    if(verbose_arrays):
        print("Array of dates unsorted:", date_list )
        tarr = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in date_list], format='isot', scale='utc')
        plt.plot(tarr.mjd)
        plt.show()
    
    #we can naively calculate a midpoint
    midpoint = date_list[0] + (date_list[-1] - date_list[0])/2 
    
    if(verbose_midpoint):
        print("Time Midpoint Naive:", midpoint)
        
    #However there is a problem. The first and last are not necessarily the beginning and end of the time series
    #So we need a sort. Naively we could do .sort(). If that doesnt work try numpy
    
    date_list = np.sort(date_list)

    if(verbose_arrays):
        #print("Array of dates sorted:", date_list.sort()) #does not work
        print("Array of dates first sort:", date_list)
        tarr = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in date_list], format='isot', scale='utc')
        plt.plot(tarr.mjd)
        plt.show()
        
    if(verbose_midpoint):
        print("Time Midpoint First Sorted:", date_list[0] + (date_list[-1] - date_list[0])/2 )
    
    #Now can we return it? Not quite. If you think about it some of the time is measured before midnight. other time is measured after midnight but referring to the day before. we need to correct for those data after midnight by adding one day
    
    before_midnight = []
    after_midnight = []
    
    for tstamps in t_obs:
        #print(tstamps[:] , tstamps[:2], int(tstamps[:2])  ) #why does int(tstamps[:]) give you error?
        if(int(tstamps[:2])>= 0 and int(tstamps[:2])< 5):
            after_midnight.append(tstamps[:])
        else:
            before_midnight.append(tstamps[:])
            
    #We can recycle the same code get list of date times
    before_midnight = [datetime.strptime(ts, "%Y-%m-%d_%H:%M:%S.%f") for ts in [result_date.group(1)+"_"+x for x in np.array(before_midnight)]]
    
    after_midnight = [datetime.strptime(ts, "%Y-%m-%d_%H:%M:%S.%f") for ts in [result_date.group(1)+"_"+x for x in np.array(after_midnight)]]
    #correct by one day. Import time delta for the correction
    #after_midnight_correction = [ty + timedelta(days=1) for ty in after_midnight]

    #The next intuitive step would be to join the two lists. We can do that like this:
    #tot_dates = before_midnight + after_midnight
    
    #However, again if you think about it what if there are data taken purely after midnight. So then we consider 3 possible cases 
    #all data before midnight 
        #len(before_midnight) != 0, len(after_midnight) = 0
    #data before midnight lasting until after midnight 
        #len(before_midnight) != 0, len(after_midnight) != 0
    #all data after midnight
        #len(before_midnight) = 0, len(after_midnight) != 0
        
    if(len(before_midnight) != 0 and len(after_midnight) != 0):
        #Import time delta for the correction
        after_midnight = [ty + timedelta(days=1) for ty in after_midnight]
        tot_dates = before_midnight + after_midnight
    elif(len(before_midnight) == 0 and len(after_midnight) != 0):
        tot_dates = after_midnight
    elif(len(before_midnight) != 0 and len(after_midnight) == 0):
        tot_dates = before_midnight 
        
    if(verbose_arrays):
        #print("Array of dates sorted:", date_list.sort()) #does not work
        print("Array of dates midnight sorted:", tot_dates)
        tarr = Time([x.strftime("%Y-%m-%dT%H:%M:%S.%f") for x in tot_dates], format='isot', scale='utc')
        plt.plot(tarr.mjd)
        plt.show()
        
    #Do a final sort
    tot_dates = np.sort(tot_dates)
    
    if(verbose_midpoint):
        print("Time Binned Midpoint Sort:", tot_dates[0] + (tot_dates[-1] - tot_dates[0])/2 )
        
    #We can do a valuidation check
    if(verbose_val_checks):
        val_verify = len(before_midnight) + len(after_midnight)
        print(str(len(before_midnight) ), "+", str(len(after_midnight)), 
              val_verify, len(tot_dates ), 
              'Error' if val_verify != len(tot_dates ) else '')
        
    #Now we can return it    
    ret_dict = {result_date.group(1)+"_"+ result_obj.group(1): (q_arr,  q_err_arr, u_arr, u_err_arr)} #4 Finally can remove

    return [ret_dict, midpoint]
    

def data_loader_sample():
    """
    Loads all the data. Returns a tuple
    """
    
    ret_list_data = []
    target_tstamps = []
    
    print("Load all excel data") #One string the path is everything before master. the filename is after.
    
    #Loading single data without the changes

    #Load single data loader.

    xlfile = "master_2020-03-05_eecep_P1-P3R0-102_combined_corr.xlsx"
    filedir = "./stats/2020-03-05/target/EE_Cep/"
    filename = filedir + xlfile

    return_data = load_pol_data_demo_exr(filename, verbose_val_checks=True)
    
    ret_list_data.append(return_data )

    xlfile = "master_2020-04-01_EE_Cep_P1-P3R307-425_man_comb_corr.xlsx"
    filedir = "./stats/2020-04-01/target/EE_Cep/"
    filename = filedir + xlfile

    return_data = load_pol_data_demo_exr(filename, verbose_val_checks=True)
    
    ret_list_data.append(return_data )

    xlfile = "master_2020-03-31_eecep_P1-P3R140-194_mac_comb_corr.xlsx"
    filedir = "./stats/2020-03-31/target/EE_Cep/"
    filename = filedir + xlfile

    return_data = load_pol_data_demo_exr(filename, verbose_val_checks=True)
    
    ret_list_data.append(return_data )

    return (ret_list_data) 
    

def load_pol_data_ori(in_str, verbose=False):
    """
    #This should return a dictioanary with the key being a date.
    #
    """
    
    path_s = re.search('(.*)master_', in_str)
    f_name_s = re.search('master_(.*)', in_str)
    
    f_name = 'master_' + f_name_s.group(1)
    
    result_MJD = re.search('master_(.*)_', f_name[:19])
    result_obj = re.search('_(.*)_P', f_name[15:])
    
    if(verbose):
        print("File path:", path)
        print("File name:", f_name)
        print("MJD:", result_MJD.group(1))
        print("obj:", result_obj.group(1))
    
    file_dir = path_s.group(1)+f_name
    loc = (file_dir)

    wb = xlrd.open_workbook(loc)
    sheet = wb.sheet_by_index(0)

    q_arr = []
    q_err_arr = []
    u_arr = []
    u_err_arr = []

    for j in range(0, sheet.nrows ):
        q_arr.append(sheet.cell_value(j, 26))
        q_err_arr.append(sheet.cell_value(j, 27))
        u_arr.append(sheet.cell_value(j, 28))
        u_err_arr.append(sheet.cell_value(j, 29))
        
    ret_dict = {result_MJD.group(1)+"_"+ result_obj.group(1): (q_arr,  q_err_arr, u_arr, u_err_arr)} 
    
    return ret_dict