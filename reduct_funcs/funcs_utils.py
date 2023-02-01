import astropy
import gzip
import glob
import re
import shutil
import importlib

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy.time import Time
from copy import deepcopy as cp
from datetime import datetime, timedelta

from reduct_funcs import funcs_polarimetry
importlib.reload(funcs_polarimetry)

def FITS_CHECKER(data_str, param_check, verb_arg=False):   
    """
    A function that that checks a FITs file for quick checking. THis function is wrapped around astropy.io.fits.open(). Try using iteratively in a loop. Returns parameter that is checked.

    Parameters
    ----------
    data_str : str
        Directory filename path.
    param_check : str
         parameter to be checked
    verb_arg : bool, optional
         Print checked parameter. True by default
    """
    
    load_data = astropy.io.fits.open(data_str)
    if(verb_arg):
        print(data_str, param_check, load_data[0].header[param_check])
    
    return (load_data[0].header[param_check])

def load_pol_data(in_str,
                  verbose_fname=False,  
                  verbose_dframe=False, 
                  verbose_midnight=False,
                  verbose_arrays = False,
                  verbose_val_checks = False,
                  verbose_midnight_calc=False):
    """
    A function that takes in a relative path filename and directory string and returns a list containing a dictionary of data in the standard format for polarimetric analysis ({'<date>_<obj_name>':(q, q_data, u, u_data)}) and the midpoint of all time obs in the file. Used recursively in util funcs to load data. Returns a dictionary of list data.

    Parameters
    ----------
    in_str : str
        Single excel q and u stats file name with their relative directory path
    verbose_fname : bool, optional
        Prints some file name things for verbosity. By default false.  
    verbose_dframe : bool, optional
        Prints some pandas data frame for verbosity. By default false.  
    verbose_midnight : bool, optional
        Prints some t_obs for either before or after midnight for verbosity. By default false.
    verbose_arrays : bool, optional
        Prints some arrays for verbosity. By default false.  
    verbose_val_checks : bool, optional
        Prints some array validation checks for verbosity. By default false.  
    verbose_midnight_calc : bool, optional
        Prints some measurements of t_obs midpoint for verbosity. By default false.
    """
    
    path_s = re.search('(.*)master_', in_str)
    f_name_s = re.search('master_(.*)', in_str)
    
    f_name = 'master_' + f_name_s.group(1)
    
    result_MJD = re.search('master_(.*)_', f_name[:19])
    result_obj = re.search('_(.*)_P', f_name[15:])
    
    df1 = pd.read_excel(in_str,engine='openpyxl')
    
    if(verbose_fname):
        print("File path:", path_s)
        print("File name:", f_name)
        print("MJD:", result_MJD.group(1)) #
        print("obj:", result_obj.group(1))
        print("\n")
    
    #TODO: if the last row has nans in the first 34 columns drop.
    #TODO: if their are nans in the nrowsx34 column drop the row.
    df1.drop(df1.tail(1).index,inplace=True) #Currently drops last row

    if(verbose_dframe):
        print("DataFrame:", df1)
        print("Shape:", df1.shape)
    
    t_obs = np.array(df1['time obs'])
    q_arr = np.array(df1['q'])
    q_err_arr = np.array(df1['q error'])
    u_arr = np.array(df1['u'])
    u_err_arr = np.array(df1['u error'])
       
    before_midnight = []
    after_midnight = []
    
    for tstamps in t_obs:
        if(verbose_midnight):
            print(tstamps[:] , tstamps[:2], int(tstamps[:2])  )

        if(int(tstamps[:2])>= 0 and int(tstamps[:2])< 6):
            if(verbose_midnight):
                print("After midnight")
            after_midnight.append(tstamps[:]) 
        else:
            if(verbose_midnight):
                print("Before midnight")
            before_midnight.append(tstamps[:])

    #all obs before midnight 
        #len(before_midnight) != 0, len(after_midnight) = 0
    #obs before midnight lasting until after midnight 
        #len(before_midnight) != 0, len(after_midnight) != 0
    #all obs after midnight
        #len(before_midnight) = 0, len(after_midnight) != 0
            
    dates_before_midnight = [datetime.strptime(ts, "%Y-%m-%d_%H:%M:%S.%f") for ts in [result_MJD.group(1)+"_"+x for x in np.array(before_midnight)]]
    
    dates_after_midnight = [datetime.strptime(ts, "%Y-%m-%d_%H:%M:%S.%f") for ts in [result_MJD.group(1)+"_"+x for x in np.array(after_midnight)]]
    
    #you only need to do this if the obs last from before midnight lasting to after midnight
    if(len(dates_before_midnight) != 0 and len(dates_after_midnight) != 0):
        dates_after_midnight_plus_one = [ty + timedelta(days=1) for ty in dates_after_midnight]
        tot_dates = dates_before_midnight + dates_after_midnight_plus_one
    elif(len(dates_before_midnight) == 0 and len(dates_after_midnight) != 0):
        tot_dates = dates_after_midnight
    elif(len(dates_before_midnight) != 0 and len(dates_after_midnight) == 0):
        tot_dates = dates_before_midnight    

    if(verbose_arrays):
        print("Before Midnight:", dates_before_midnight, len(dates_before_midnight))
        print("After Midnight:", dates_after_midnight, len(dates_after_midnight))    
        print("Combined unsorted:", tot_dates , len(tot_dates ))
        
    if(verbose_val_checks):
        val_verify = len(dates_before_midnight) + len(dates_after_midnight)
        print(str(len(dates_before_midnight) ), "+", str(len(dates_after_midnight)), 
              val_verify, len(tot_dates ), 
              'FUBAR' if val_verify != len(tot_dates ) else '')
    
    dates = sorted(tot_dates)    
    
    midpoint = dates[0] + (dates[-1] - dates[0])/2
    if(verbose_midnight_calc):
        print("Time Midpoint:", midpoint)

    ret_dict = {result_MJD.group(1)+"_"+ result_obj.group(1): (np.array(df1['q']),  np.array(df1['q error']), np.array(df1['u']), np.array(df1['u error']))}

    return [ret_dict,  midpoint]

def data_splitter(list_data):
    """
    A function that takes in a list of data from data loader. Returns a tuple with a list of data in [0] and a list of MJDs in [1].

    Parameters
    ----------
    list_data : list
        Single excel q and u stats file name with their relative directory path
    """
    height, width = np.array(list_data).shape

    # Cut the image in half
    width_cutoff = width // 2
    s1 = np.array(list_data)[:, :width_cutoff] #data array
    s2 = np.array(list_data)[:, width_cutoff:] #MJD array
    
    return ([s1.flatten(), s2.flatten()])

def make_set_list(s_query, param_check, verb_arg, verb_arg_MJD=False):
    """
    A function that takes in a list of data from data loader. Returns a tuple with a list of data in [0] and a list of MJDs in [1].

    Parameters
    ----------
    list_data : list
        Single excel q and u stats file name with their relative directory path
    """
    master_list_device = []
    data_dict = {}
    
    list_MJDs = [word for word in glob.glob('./files_sorted/*') if len(word) >= 25 ]

    for MJDs in list_MJDs[:]:  #For each MJD
        if(verb_arg_MJD):
            print(MJDs, MJDs[15:])
        mjd_list = []

        target = glob.glob(MJDs +'/' + s_query  + '/*')#Get all the data
        for targs in target:
            y = FITS_CHECKER(target[0], 'DEVICE', verb_arg) #Get the device. You cant call it y
            master_list_device.append(y) #and put into master list
            mjd_list.append(y)
            
        #Append here:
        data_dict.update({MJDs[15:]:mjd_list})

    return(master_list_device, data_dict)

def date_split(return_d):
    """
    Date split
    """
    
    c_z = 0
    dnew = {}
    for key in return_d[1]:

        print("Date:", key) #
        unique_labels, unique_counts = np.unique(return_d[1][key], return_counts=True)
        for l in range(0 , len(unique_labels)):
            print("Camera:", unique_labels[l])
            print("Count:" , unique_counts[l])
        c_z += unique_counts
        print("\n")   

    print(c_z)

def correct_time(Starg, 
                 verbose_s = False, 
                 verbose_t = False, 
                 verbose_u = False):
    """
    A function that corrects time based on FITS files

    Parameters
    ----------
    Starg : str
        Single excel q and u stats file name with their relative directory path
    verbose_s : bool, optional
        Prints some file name things for verbosity. By default false.  
    verbose_t : bool, optional
        Prints some pandas data frame for verbosity. By default false.  
    verbose_u : bool, optional
        Prints some t_obs for either before or after midnight for verbosity. By default false.
    """
    
    long_list = []
    dict_my = {}
    dict_mjd = {}
    c = 0
    
    list_MJDs = [word for word in glob.glob('./files_sorted/*') if len(word) >= 25 ] #Looks through files sorted
    
    for MJDs in list_MJDs[:]: #The big loop. Does for all.
        mjd_list = []    
        tot_exp_time = [ ]
        if(Starg == 'EECEP'):
            target = glob.glob(MJDs +'/' + 'ee_cep'  + '/*') #Get all the EECEP FITS data            
        elif(Starg == 'G191'):
            target = glob.glob(MJDs +'/' + 'pol_std'  + '/*') #Get all the 191 FITS data
            if(verbose_u):
                print('Dis target:', target)
            c_t = target
            nwtarget = []
            for things in target:
                if '191' in things:
                    nwtarget.append(things)
            target = nwtarget
            if(verbose_t):
                print("Total pol_stds:",len(c_t), "\n", 
                      "Total G191B2B:", len(nwtarget),len(target), "Percentage is 191:", target/len(c_t)  )
        elif(Starg == '2158'):
            target = glob.glob(MJDs +'/' + 'pol_std'  + '/*') #Get all the 191 FITS data
            
            if(verbose_u):
                print('Dis target:', target)
            
            c_t = target
            nwtarget = []
            for things in target:
                if '2158' in things:
                    nwtarget.append(things)
            target = nwtarget
            if(verbose_t):
                print("Total pol_stds:",len(c_t), "\n", 
                      "Total G191B2B:", len(nwtarget),len(target), "Percentage is 191:", target/len(c_t))                
        elif(Starg == '64106'):
            target = glob.glob(MJDs +'/' + 'pol_std'  + '/*') #Get all the 191 FITS data
            
            if(verbose_u):
                print('Dis target:', target)
            
            c_t = target
            nwtarget = []
            for things in target:
                if '64106' in things:
                    nwtarget.append(things)
            target = nwtarget
            if(verbose_t):
                print("Total pol_stds:",len(c_t), "\n", 
                      "Total G191B2B:", len(nwtarget),len(target), "Percentage is 191:", target/len(c_t)  )
        
        if(verbose_s):
            print(c, "Day:" , MJDs[len('./files_sorted\\'):], len(target), "FITS files = ", len(target)/2 , "q n u points")
        #IT searches through all the files sorted

        for k in range(0, len(target)):
            #Openning will take time ofcourse
            met = astropy.io.fits.open(target[k])[0].header['DATE-OBS']+'T'+astropy.io.fits.open(target[k])[0].header['TIME-OBS']
            exptime = astropy.io.fits.open(target[k])[0].header['EXPTIME']
            t = Time(met, format='isot', scale='utc')

            if(verbose_t):
                print("Exposure time:", astropy.io.fits.open(target[k])[0].header['EXPTIME']) 
                print("FITS JD TIME:", astropy.io.fits.open(target[k])[0].header['JD'])                   
                print("Calculated JD TIME:", t.jd)
                print("Calculated MJD TIME:", t.mjd, "\n")

            long_list.append(t.mjd)
            tot_exp_time.append(exptime)
            mjd_list.append(t.mjd)

        dict_my.update({MJDs[15:]:mjd_list})
        dict_mjd.update({MJDs[15:]:tot_exp_time}) #it doent have to be 
        c+=1

    corr_dict = {} #This gets returned
    mjd_coo = {}

    c = 0
    
    for keyes in dict_my.keys():
        corr_dict.update({keyes: [np.median(np.sort(dict_my[keyes])) , np.std(dict_my[keyes])]  })
        if(verbose_t): #this is kind of like a final summary
            print(c, "Date:", keyes, "Time center:" , np.median(np.sort(dict_my[keyes])) ,u"\u00B1", np.std(dict_my[keyes]))
        c+=1
        
    return (corr_dict, dict_mjd)
    
def combine_excels(excel_1, excel_2, sv_out , MJD, targ_obj, obs_filt ,strt_ind, end_ind):
    """
    Function that takes in two excel scripts and combines them into the usual 1. Does not work since 
    
    One last chance to do this bro. Today 23/10/2021. Do it.
    """
    #create new workbookq
    
    print("Sheet 1:",excel_1)
    print("Sheet 2:",excel_2)
    
    filename= targ_obj+"_P1-P3"+obs_filt
    workbook = xlsxwriter.Workbook(sv_out+'master_'+MJD+"_"+ filename+str(strt_ind)+"-" +str(end_ind)+'_mac_comb.xlsx')
    worksheet = workbook.add_worksheet()
    
    print('output:', sv_out+'master_'+MJD+"_"+filename+str(strt_ind)+"-" +str(end_ind)+'_mac_comb.xlsx')
    
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
    
    worksheet.write('N1', 'int')
    worksheet.write('O1', 'time obs')
    worksheet.write('P1', 'filter')
    worksheet.write('Q1', 'exptime')
    worksheet.write('R1', 'aperture radius')
    worksheet.write('S1', 'target 1 x center')
    worksheet.write('T1', 'target 1 y center')
    worksheet.write('U1', 'target 1 counts')
    worksheet.write('V1', 'target 1 error')
    worksheet.write('W1', 'target 2 x center')
    worksheet.write('X1', 'target 2 y center')
    worksheet.write('Y1', 'target 2 counts')
    worksheet.write('Z1', 'target 2 error')

    worksheet.write('AA1', 'q')
    worksheet.write('AB1', 'q error')
    worksheet.write('AC1', 'u')
    worksheet.write('AD1', 'u error')
    
    #Just calc. Previously not calculated but we should now calculate anyway
    worksheet.write('AE1', 'PD') 
    worksheet.write('AF1', 'PD error')
    worksheet.write('AG1', 'PA')
    worksheet.write('AH1', 'PA error')
    
    wb_1 = xlrd.open_workbook((excel_1))
    wb_2 = xlrd.open_workbook((excel_2))
    
    sheet_1 = wb_1.sheet_by_index(0)
    sheet_2 = wb_2.sheet_by_index(0)
    
    #This reokucates
    for u in range(1, np.min([sheet_1.nrows,   sheet_2.nrows      ])):
        for x in range(0, 13):
            worksheet.write(u, x, sheet_1.cell_value(rowx=u, colx=x))
        for x in range(0, 13):
            worksheet.write(u, x+13, sheet_2.cell_value(rowx=u, colx=x))
            
    for u in range(1, np.min([sheet_1.nrows,   sheet_2.nrows      ])):
        
        #First do q and q error
        p_a = 4*(sheet_2.cell_value(rowx=u, colx=7)**2)*(sheet_2.cell_value(rowx=u, colx=12)**2) + (sheet_2.cell_value(rowx=u, colx=11)**2)*(sheet_2.cell_value(rowx=u, colx=8)**2)
        p_b = (sheet_2.cell_value(rowx=u, colx=7) + sheet_2.cell_value(rowx=u, colx=11))**4
        p_ans = np.sqrt(p_a/p_b)
        
        worksheet.write(u, 26, (sheet_2.cell_value(rowx=u, colx=7)-sheet_2.cell_value(rowx=u, colx=11))/(sheet_2.cell_value(rowx=u, colx=7)+sheet_2.cell_value(rowx=u, colx=11))   )
        worksheet.write(u, 27, p_ans)       

        #Then do u and u error
        p_a = 4*(sheet_1.cell_value(rowx=u, colx=7)**2)*(sheet_1.cell_value(rowx=u, colx=12)**2) + (sheet_1.cell_value(rowx=u, colx=11)**2)*(sheet_1.cell_value(rowx=u, colx=8)**2)
        p_b = (sheet_1.cell_value(rowx=u, colx=7) + sheet_1.cell_value(rowx=u, colx=11))**4
        p_ans = np.sqrt(p_a/p_b)
        
        #Then do PD
        
        #Then do PA
        
        worksheet.write(u, 28, (sheet_1.cell_value(rowx=u, colx=7)-sheet_1.cell_value(rowx=u, colx=11))/(sheet_1.cell_value(rowx=u, colx=7)+sheet_1.cell_value(rowx=u, colx=11))) #write q, q_err
        worksheet.write(u, 29, p_ans  ) #write q, q_err
    
    workbook.close() #I guess this writes the thing

def list_autoloader(input_string, 
                    verbose_file=False, 
                    verbose_filename=False,
                    verbose_process=False):
    """
    #A function that takes in a list
    """
    list_return = []
    
    if(verbose_file):
        print("Load all the files in this text file:", input_string)
    
    f = open(input_string, "r")
    
    for things in f.readlines():
        if(verbose_filename):
            print("Loading:", things)
            
        pol_data = load_pol_data(things.strip()) #this returns a tuple  
        list_return.append(pol_data)

    if(verbose_process):
        print(len(list_return), "lists processed.")
        
    splitted_data = data_splitter(list_return)

    return (splitted_data)
    
def autoloader(verb_arg):
    """
    #Not working. 31/07/2021. Still not working. 04/01/2023.
    #Autloader using glob and a bunch of exceptions
    """
    
    if(verb_arg):
        print("That old way of doing things sucked")
        print("This is the way\n")
    
    x = glob.glob('./stats/*')
    
    ret_list_zero_pol = []
    ret_list_high_pol = []
    ret_list_target = []

    for things in x:
        path_s = re.search('./stats\\\(.*)', things)
        if(verb_arg):
            print(path_s.group(1))
        y = glob.glob(things+'/*')

        for blah in y:
            z = glob.glob(blah+'/*')
            for jah in z:
                xl = glob.glob(jah+'/*') #get all xls
                load_this = ''
                sub_xl = [ x for x in xl if "comb" in x ]
                if(verb_arg):
                    print(sub_xl)
                if len(sub_xl) > 1:
                    load_this = [j for j in sub_xl if "corr" in j][0]
                else:
                    load_this = sub_xl[0]

                if(verb_arg):
                    print("Loading this:", load_this)
                    
                if('ee' in load_this or 'EE' in load_this or 'Cep' in load_this or 'Cet' in load_this or 'cep' in load_this or 'Cep' in load_this):
                    target = funcs_polarimetry.load_pol_data(load_this, False)
                    ret_list_target.append(target)
                elif('191' in load_this or '212311' in load_this):
                    low_pol = funcs_polarimetry.load_pol_data(load_this, False)       
                    ret_list_zero_pol.append(low_pol)
                elif('215806' in load_this or '287' in load_this or '204827' in load_this or '251204' in load_this or '64106' in load_this):
                    high_pol = funcs_polarimetry.load_pol_data(load_this, False)
                    ret_list_high_pol.append(high_pol)
            
                #print("combtainer:", combtainer, len(combtainer)) #now go throgh container 
                #If the list contains an item with the word corr. Drop everything else except the item with the word corr
                #else if there is not item with the word corr. Just take the thing with comb (else do nothing)
                #if 'g191b2b' in sha or 'hd212311' in sha:
                #print("Low Polarization standard:", sha)
               
            
            
    #Check this out!
    #Suggestion: Create a text file that holds a list of text files for different data types and use preexisting
    #developed load pol data code to iterate through this list and the same product as data_loader
    #This would make good practice. Write and demo your own data loader that has to handle different data types 
    #and this will be a step towards data mastery.

    
    
    return (ret_list_target, ret_list_zero_pol, ret_list_high_pol)
    
def sample_data_loader():
    """
    #Load small sample
    """
       
    ret_list_zero_pol = []
    ret_list_high_pol = []
    ret_list_target = []
    
    print("Load a sample of excel data")
    
    #2020_03_05
    zero_pol_std_03_05 = funcs_polarimetry.load_pol_data("./stats/GRN_03_method_check/master_2020-03-05_g191b2b_P1-P3R0-46_mac_comb.xlsx")
    
    ret_list_zero_pol.append(zero_pol_std_03_05)
    
    #2020_03_14
    zero_pol_std_03_14 = funcs_polarimetry.load_pol_data("./stats/GRN_03_method_check/master_2020-03-14_g191b2b_P1-P3R0-60_mac_comb.xlsx")
    
    ret_list_zero_pol.append(zero_pol_std_03_14)
    
    #2020_03_15
    zero_pol_std_03_15 = funcs_polarimetry.load_pol_data("./stats/GRN_03_method_check/master_2020-03-15_G191B2B_P1-P3R0-40_mac_comb.xlsx")
    
    ret_list_zero_pol.append(zero_pol_std_03_15)
    
    #2020_03_26
    zero_pol_std_03_26 = funcs_polarimetry.load_pol_data("./stats/GRN_03_method_check/master_2020-03-26_g191b2b_P1-P3R0-30_mac_comb.xlsx")
    
    ret_list_zero_pol.append(zero_pol_std_03_26)
    
    #2020_04_16
    zero_pol_std_04_16 = funcs_polarimetry.load_pol_data("./stats/GRN_03_method_check/master_2020-04-16_g191b2b_P1-P3R62-134_mac_comb.xlsx")
    
    ret_list_zero_pol.append(zero_pol_std_04_16)
    
    #2020_05_01
    zero_pol_std_05_01 = funcs_polarimetry.load_pol_data("./stats/GRN_03_method_check/master_2020-05-01_g191b2b_P1-P3R102-126_mac_comb.xlsx")
    
    ret_list_zero_pol.append(zero_pol_std_05_01)
    
    #2020_05_04
    zero_pol_std_05_04 = funcs_polarimetry.load_pol_data("./stats/GRN_03_method_check/master_2020-05-04_g191b2b_P1-P3R24-44_mac_comb.xlsx")
    
    ret_list_zero_pol.append(zero_pol_std_05_04)
    
    #2020_07_29
    target_data_07_29 = funcs_polarimetry.load_pol_data("./stats/GRN_EE_Cep_Off_Eclipse/master/master_2020-07-29_eecep_P1-P3R360-484_mac_comb.xlsx")
    
    ret_list_target.append(target_data_07_29)
    
    #2020_07_25
    target_data_08_25 = funcs_polarimetry.load_pol_data("./stats/GRN_EE_Cep_Off_Eclipse/master/master_2020-08-25_eecep_P1-P3R388-534_mac_comb.xlsx")
    
    ret_list_target.append(target_data_07_29)
    
    #2020_09_08
    target_data_09_08 = funcs_polarimetry.load_pol_data("./stats/GRN_EE_Cep_Off_Eclipse/master/master_2020-09-08_eecep_P1-P3R212-332_mac_comb.xlsx")
    
    ret_list_target.append(target_data_07_29)
    
    #2020_09_22
    target_data_09_22 = funcs_polarimetry.load_pol_data("./stats/GRN_EE_Cep_Off_Eclipse/master/master_2020-09-22_eecep_P1-P3R200-320_mac_comb.xlsx")
    
    ret_list_target.append(target_data_07_29)
    
    return (ret_list_target, ret_list_zero_pol, ret_list_high_pol)
    
def data_loader():
    """
    Loads all the data. Returns a tuple
    """
    
    ret_list_zero_pol = []
    ret_list_high_pol = []
    ret_list_target = []
    
    target_tstamps = []
    zero_pol_tstamps = []
    high_pol_tstamps = []
    
    print("Load all excel data") #One string the path is everything before master. the filename is after.
    
    ##########
    #2020-03-05
    
    target_data_03_05 = load_pol_data( "./stats/2020-03-05/target/EE_Cep/master_2020-03-05_eecep_P1-P3R0-102_combined_corr.xlsx")
    
    zero_pol_std_03_05 = load_pol_data("./stats/2020-03-05/pol_std/G191B2B/master_2020-03-05_g191b2b_P1-P3R0-46_combined_corr.xlsx")
    
    high_pol_std_03_05 = load_pol_data("./stats/2020-03-05/pol_std/HD215806/master_2020-03-05_hd215806_P1-P3R46-106_combined_corr.xlsx")
    
    ret_list_target.append(target_data_03_05)
    ret_list_zero_pol.append(zero_pol_std_03_05)
    ret_list_high_pol.append(high_pol_std_03_05)
    
    ##########
    #2020-03-13
    #I am making a fiducial date for trouble shooting. 2020-03-13
    #There is shit in here that should not even be loaded. This is fake as sheet.
    
    #zero_pol_std_03_13 = funcs_polarimetry.load_pol_data( "./stats/2020-03-13/pol_std/G191B2B/master_2020-03-13_G191B2B_P1-P3R0-60_mac_comb.xlsx", False) #This is unstable. Try Relading it
    
    #ret_list_zero_pol.append(zero_pol_std_03_13) #worse
    
    ##########
    #2020-03-14
    
    target_data_03_14 = load_pol_data( "./stats/2020-03-14/target/EE_Cep/master_2020-03-14_EE Cep_P1-P3R0-69_combined.xlsx")
    
    zero_pol_std_03_14 = load_pol_data( "./stats/2020-03-14/pol_std/G191B2B/master_2020-03-14_G191B2B_P1-P3R0-60_combined_corr.xlsx") #This is unstable. Try Relading it
    
    high_pol_std_03_14 = load_pol_data("./stats/2020-03-14/pol_std/HD215806/master_2020-03-14_HD215806_P1-P3R60-301_combined_corr.xlsx")

    ret_list_target.append(target_data_03_14)
    ret_list_zero_pol.append(zero_pol_std_03_14) #bad
    ret_list_high_pol.append(high_pol_std_03_14)
    
    ##########
    #2020-03-15
    
    target_data_03_15 = load_pol_data("./stats/2020-03-15/target/EE_Cep/master_2020-03-15_EE Cep_P1-P3R19-59_combined.xlsx")
    
    zero_pol_std_03_15 = load_pol_data("./stats/2020-03-15/pol_std/G191B2B/master_2020-03-15_G191B2B_P1-P3R0-40_combined_corr.xlsx")
    
    high_pol_std_03_15 = load_pol_data("./stats/2020-03-15/pol_std/HD215806/master_2020-03-15_HD215806_P1-P3R40-109_combined.xlsx")
    
    ret_list_target.append(target_data_03_15)
    ret_list_zero_pol.append(zero_pol_std_03_15)
    ret_list_high_pol.append(high_pol_std_03_15)
    
    ##########
    #2020-03-24

    target_data_03_24 = load_pol_data("./stats/2020-03-24/target/EE_Cep/master_2020-03-24_eecep_P1-P3R600-770_mac_comb_corr.xlsx")
    
    ret_list_target.append(target_data_03_24)
    
    ##########
    #2020-03-25
    
    target_data_03_25 = load_pol_data("./stats/2020-03-25/target/EE_Cep/master_2020-03-25_eecep_P1-P3R464-540_combined_corr.xlsx")
    
    zero_pol_std_03_25 = load_pol_data("./stats/2020-03-25/pol_std/HD212311/master_2020-03-25_hd212311_P1-P3R0-44_combined_corr.xlsx")
    
    ret_list_target.append(target_data_03_25)
    ret_list_zero_pol.append(zero_pol_std_03_25)
    
    ##########
    #2020-03-26

    target_data_03_26 = load_pol_data("./stats/2020-03-26/target/EE_Cep/master_2020-03-26_eecep_P1-P3R40-109_mac_comb_corr.xlsx")
    
    zero_pol_std_03_26 = load_pol_data("./stats/2020-03-26/pol_std/G191B2B/master_2020-03-26_g191b2b_P1-P3R0-30_combined_corr.xlsx")
    zero_pol_std_03_26_02 = load_pol_data("./stats/2020-03-26/pol_std/HD212311/master_2020-03-26_hd212311_P1-R90-130_combined_corr.xlsx")
    
    #Broken XL... it was err not error
    high_pol_std_03_26 = load_pol_data("./stats/2020-03-26/pol_std/oj287/master_2020-03-26_oj287_P1-P3R258-296_comb_corr.xlsx")
    high_pol_std_03_26_02 = load_pol_data("./stats/2020-03-26/pol_std/HD204827/master_2020-03-26_hd204827_P1-P3R60-90_combined_corr.xlsx")
    high_pol_std_03_26_03 = load_pol_data("./stats/2020-03-26/pol_std/HD215806/master_2020-03-26_HD215806_P1-P3R40-109_mac_comb_corr.xlsx")
    high_pol_std_03_26_04 = load_pol_data("./stats/2020-03-26/pol_std/HD251204/master_2020-03-26_hd251204_P1-P3R40-109_mac_comb.xlsx")
    
    ret_list_target.append(target_data_03_26)
    
    ret_list_zero_pol.append(zero_pol_std_03_26)
    ret_list_zero_pol.append(zero_pol_std_03_26_02)
    
    #ret_list_high_pol.append(high_pol_std_03_26)
    ret_list_high_pol.append(high_pol_std_03_26_02)
    ret_list_high_pol.append(high_pol_std_03_26_03)
    ret_list_high_pol.append(high_pol_std_03_26_04)
    
    ##########
    #2020-03-31

    target_data_03_31 = load_pol_data("./stats/2020-03-31/target/EE_Cep/master_2020-03-31_eecep_P1-P3R140-194_mac_comb_corr.xlsx")
    
    zero_pol_std_03_31 = load_pol_data("./stats/2020-03-31/pol_std/HD212311/master_2020-03-31_hd212311_P1-P3R0-40_mac_comb_corr.xlsx")
    
    ret_list_target.append(target_data_03_31)
    ret_list_zero_pol.append(zero_pol_std_03_31)
    
    ##########
    #2020-04-01

    target_data_04_01 = load_pol_data("./stats/2020-04-01/target/EE_Cep/master_2020-04-01_EE_Cep_P1-P3R307-425_man_comb_corr.xlsx")
    
    ret_list_target.append(target_data_04_01)
    
    ##########
    #2020-04-02
    
    target_data_04_02 = load_pol_data("./stats/2020-04-02/target/EE_Cep/master_2020-04-02_eecep_P1-P3R180-230_mac_comb.xlsx")
    
    ret_list_target.append(target_data_04_02)
    
    ##########
    #2020-04-04

    target_data_04_04 = load_pol_data("./stats/2020-04-04/target/EE_Cep_test_x/master_2020-04-04_EE Cep_P1-P3R0-46_mac_comb_corr.xlsx")
    
    ret_list_target.append(target_data_04_04)
    
    ##########
    #2020-04-05

    target_data_04_05 = load_pol_data("./stats/2020-04-05/target/EE_Cep/master_2020-04-05_EE Cep_P1-P3R0-43_comb_corr.xlsx")
    
    zero_pol_std_04_05 = load_pol_data("./stats/2020-04-05/pol_std/HD212311/master_2020-04-05_HD212311_P1-P3R0-40_mac_comb_corr.xlsx")
    
    high_pol_std_04_05 = load_pol_data("./stats/2020-04-05/pol_std/HD215806/master_2020-04-05_hd215806_P1-P3R40-160_mac_comb.xlsx")
    
    ret_list_target.append(target_data_04_05)
    
    ret_list_zero_pol.append(zero_pol_std_04_05)
    
    ret_list_high_pol.append(high_pol_std_04_05)
    
    ##########
    #2020-04-06

    target_data_04_06 = load_pol_data("./stats/2020-04-06/target/EE_Cep/master_2020-04-06_eecep_P1-P3R428-618_mac_comb_corr.xlsx")
    
    zero_pol_std_04_06 = load_pol_data("./stats/2020-04-06/pol_std/HD212311/master_2020-04-06_hd212311_P1-P3R0-180_mac_comb_corr.xlsx")
    
    high_pol_std_04_06 = load_pol_data("./stats/2020-04-06/pol_std/HD215806/master_2020-04-06_hd215806_P1-P3R180-388_mac_comb_corr.xlsx")
    high_pol_std_04_06_02 = load_pol_data("./stats/2020-04-06/pol_std/oj287/master_2020-04-06_oj287_P1-P3R388-402_mac_comb.xlsx")


    ret_list_target.append(target_data_04_06)
    
    ret_list_zero_pol.append(zero_pol_std_04_06)

    ret_list_high_pol.append(high_pol_std_04_06)
    ret_list_high_pol.append(high_pol_std_04_06_02)
    
    ##########
    #2020-04-07
    
    target_data_04_07 = load_pol_data("./stats/2020-04-07/target/EE_Cep/master_2020-04-07_eecep_P1-P3R504-712_mac_comb_corr.xlsx")
    
    zero_pol_std_04_07 = load_pol_data("./stats/2020-04-07/pol_std/HD212311/master_2020-04-07_hd212311_P1-P3R0-180_mac_comb.xlsx")
    
    high_pol_std_04_07 = load_pol_data("./stats/2020-04-07/pol_std/HD215806/master_2020-04-07_hd215806_P1-P3R180-300_mac_comb_corr.xlsx")
    high_pol_std_04_07_02 = load_pol_data("./stats/2020-04-07/pol_std/oj287/master_2020-04-07_oj287_P1-P3R300-312_mac_comb.xlsx")
    
    ret_list_target.append(target_data_04_07)
    
    ret_list_zero_pol.append(zero_pol_std_04_07)
    
    ret_list_high_pol.append(high_pol_std_04_07)
    ret_list_high_pol.append(high_pol_std_04_07_02)
    
    ##########
    #2020-04-08
    
    
    target_data_04_08 = load_pol_data("./stats/2020-04-08/target/EE_Cep/master_2020-04-08_eecep_P1-P3R288-508_mac_comb_corr.xlsx")
    
    zero_pol_std_04_08 = load_pol_data("./stats/2020-04-08/pol_std/HD212311/master_2020-04-08_hd212311_P1-P3R344-538_mac_comb_corr.xlsx")
    
    high_pol_std_04_08 = load_pol_data("./stats/2020-04-08/pol_std/HD204827/master_2020-04-08_hd204827_P1-P3R42-276_mac_comb.xlsx")
    high_pol_std_04_08_02 = load_pol_data("./stats/2020-04-08/pol_std/HD215806/master_2020-04-08_hd215806_P1-P3R589-741_mac_comb_corr.xlsx")
    
    ret_list_target.append(target_data_04_08)
    
    ret_list_zero_pol.append(zero_pol_std_04_08)
    
    ret_list_high_pol.append(high_pol_std_04_08)
    ret_list_high_pol.append(high_pol_std_04_08_02)
    
    ##########
    #2020-04-15
    
    target_data_04_15 = load_pol_data("./stats/2020-04-15/target/EE_Cep/master_2020-04-15_eecep_P1-P3R320-426_mac_comb.xlsx")
    
    zero_pol_std_04_15 = load_pol_data("./stats/2020-04-15/pol_std/HD212311/master_2020-04-15_hd212311_P1-P3R120-340_mac_comb.xlsx")
    
    high_pol_std_04_15 = load_pol_data("./stats/2020-04-15/pol_std/HD204827/master_2020-04-15_hd204827_P1-P3R0-120_mac_comb_corr.xlsx")
    high_pol_std_04_15_02 = load_pol_data("./stats/2020-04-15/pol_std/HD215806/master_2020-04-15_hd215806_P1-P3R340-442_mac_comb.xlsx")
    high_pol_std_04_15_03 = load_pol_data("./stats/2020-04-15/pol_std/OJ287/master_2020-04-15_oj287_P1-P3R503-541_mac_comb.xlsx")
    
    ret_list_target.append(target_data_04_15)
    
    ret_list_zero_pol.append(zero_pol_std_04_15)
    
    ret_list_high_pol.append(high_pol_std_04_15)
    ret_list_high_pol.append(high_pol_std_04_15_02)
    ret_list_high_pol.append(high_pol_std_04_15_03)
   
    ##########
    #2020-04-16
    
    target_data_04_16 = load_pol_data("./stats/2020-04-16/target/EE_Cep/master_2020-04-16_eecep_P1-P3R300-500_mac_comb_corr.xlsx")
    
    zero_pol_std_04_16 = load_pol_data("./stats/2020-04-16/pol_std/G191B2B/master_2020-04-16_g191b2b_P1-P3R62-134_mac_comb_corr.xlsx")
    zero_pol_std_04_16_02 = load_pol_data("./stats/2020-04-16/pol_std/HD212311/master_2020-04-16_hd212311_P1-P3R194-234_mac_comb.xlsx")
    
    high_pol_std_04_16 = load_pol_data("./stats/2020-04-16/pol_std/HD204827/master_2020-04-16_hd204827_P1-P3R134-193_mac_comb_corr.xlsx")
    high_pol_std_04_16_02 = load_pol_data("./stats/2020-04-16/pol_std/BD64106/master_2020-04-16_bd64106_P1-P3R0-62_mac_comb.xlsx")
    high_pol_std_04_16_03 = load_pol_data("./stats/2020-04-16/pol_std/HD215806/master_2020-04-16_hd215806_P1-P3R234-286_mac_comb_corr.xlsx")
    high_pol_std_04_16_04 = load_pol_data("./stats/2020-04-16/pol_std/OJ287/master_2020-04-16_oj287_P1-P3R324-344_mac_comb.xlsx")


    ret_list_target.append(target_data_04_16)
    
    ret_list_zero_pol.append(zero_pol_std_04_16)
    ret_list_zero_pol.append(zero_pol_std_04_16_02)
    
    ret_list_high_pol.append(high_pol_std_04_16)
    ret_list_high_pol.append(high_pol_std_04_16_02)
    ret_list_high_pol.append(high_pol_std_04_16_03)
    ret_list_high_pol.append(high_pol_std_04_16_04)
    
    ##########
    #2020-04-18
    
    target_data_04_18 = load_pol_data("./stats/2020-04-18/target/EE_Cep/master_2020-04-18_EE Cet_P1-P3R0-100_mac_comb_corr.xlsx")
    
    high_pol_std_04_18 = load_pol_data("./stats/2020-04-18/pol_std/oj287/master_2020-04-18_OJ287_P1-P3R0-20_mac_comb.xlsx")

    ret_list_high_pol.append(high_pol_std_04_18)
    
    ret_list_target.append(target_data_04_18)
    
    ##########
    #2020-04-19
    
    target_data_04_19 = load_pol_data("./stats/2020-04-19/target/EE_Cep/master_2020-04-19_EE Cep_P1-P3R0-50_mac_comb_corr.xlsx")
    
    zero_pol_std_04_19 = load_pol_data("./stats/2020-04-19/pol_std/HD212311/master_2020-04-19_HD212311_P1-P3R0-165_mac_comb_corr.xlsx")
    
    high_pol_std_04_19 = load_pol_data("./stats/2020-04-19/pol_std/HD215806/master_2020-04-19_hd215806_P1-P3R165-265_mac_comb_corr.xlsx")
    high_pol_std_04_19_02 = load_pol_data("./stats/2020-04-19/pol_std/OJ287/master_2020-04-19_oj287_P1-P3R265-285_mac_comb.xlsx")
        
    ret_list_zero_pol.append(zero_pol_std_04_19)
    
    ret_list_high_pol.append(high_pol_std_04_19)
    ret_list_high_pol.append(high_pol_std_04_19_02)
    
    ret_list_target.append(target_data_04_19)
    
    ##########
    #2020-04-21
    
    target_data_04_21 = load_pol_data("./stats/2020-04-21/target/EE_Cep/master_2020-04-21_eecep_P1-P3R524-586_mac_comb.xlsx")
    
    zero_pol_std_04_21 = load_pol_data("./stats/2020-04-21/pol_std/HD212311/master_2020-04-21_hd212311_P1-P3R0-60_mac_comb_corr.xlsx")
    
    high_pol_std_04_21 = load_pol_data("./stats/2020-04-21/pol_std/HD215806/master_2020-04-21_hd215806_P1-P3R60-120_mac_comb.xlsx")
    high_pol_std_04_21_02 = load_pol_data("./stats/2020-04-21/pol_std/OJ287/master_2020-04-21_oj287_P1-P3R167-209_mac_comb_corr.xlsx")
    
    ret_list_zero_pol.append(zero_pol_std_04_21)
    ret_list_high_pol.append(high_pol_std_04_21)
    ret_list_high_pol.append(high_pol_std_04_21_02)
    
    ret_list_target.append(target_data_04_21)
    
    ##########
    #2020-04-22
    
    target_data_04_22 = load_pol_data("./stats/2020-04-22/target/EE_Cep/master_2020-04-22_eecep_P1-P3R240-342_mac_comb.xlsx")
    
    zero_pol_std_04_22 = load_pol_data("./stats/2020-04-22/pol_std/HD212311/master_2020-04-22_hd212311_P1-P3R74-154_mac_comb_corr.xlsx")
    
    high_pol_std_04_22 = load_pol_data("./stats/2020-04-22/pol_std/HD215806/master_2020-04-22_hd215806_P1-P3R154-226_mac_comb_corr.xlsx")
    high_pol_std_04_22_02 = load_pol_data("./stats/2020-04-22/pol_std/BD64106/master_2020-04-22_bd64106_P1-P3R0-74_mac_comb_corr.xlsx")
    high_pol_std_04_22_03 = load_pol_data("./stats/2020-04-22/pol_std/OJ287/master_2020-04-22_oj287_P1-P3R266-288_mac_comb.xlsx")
 
    ret_list_zero_pol.append(zero_pol_std_04_22)
    
    ret_list_high_pol.append(high_pol_std_04_22)
    ret_list_high_pol.append(high_pol_std_04_22_02)
    ret_list_high_pol.append(high_pol_std_04_22_03)
    
    ret_list_target.append(target_data_04_22)
    
    ##########
    #2020-04-23

    target_data_04_23 = load_pol_data("./stats/2020-04-23/target/EE_Cep/master_2020-04-23_eecep_P1-P3R320-410_mac_comb.xlsx")
    
    zero_pol_std_04_23 = load_pol_data("./stats/2020-04-23/pol_std/HD212311/master_2020-04-23_hd212311_P1-P3R60-120_mac_comb.xlsx")
    
    high_pol_std_04_23 = load_pol_data("./stats/2020-04-23/pol_std/BD64106/master_2020-04-23_bd64106_P1-P3R0-60_mac_comb.xlsx")
    high_pol_std_04_23_02 = load_pol_data("./stats/2020-04-23/pol_std/HD215806/master_2020-04-23_hd215806_P1-P3R120-180_mac_comb_corr.xlsx")
    high_pol_std_04_23_03 = load_pol_data("./stats/2020-04-23/pol_std/OJ287/master_2020-04-23_oj287_P1-P3R229-249_mac_comb.xlsx")
        
    ret_list_zero_pol.append(zero_pol_std_04_23)
    
    ret_list_target.append(target_data_04_23)
    
    ret_list_high_pol.append(high_pol_std_04_23)
    ret_list_high_pol.append(high_pol_std_04_23_02)
    ret_list_high_pol.append(high_pol_std_04_23_03)
    
    ##########
    #2020-04-27
    #The bad batch
    
    target_data_04_27 = load_pol_data("./stats/2020-04-27/target/EE_Cep/master_2020-04-27_EECep_P1-P3R0-60_mac_comb_corr.xlsx")    
    
    zero_pol_std_04_27 = load_pol_data("./stats/2020-04-27/pol_std/HD212311/master_2020-04-27_HD212311_P1-P3R0-80_mac_comb_corr.xlsx")
    high_pol_std_04_27 = load_pol_data("./stats/2020-04-27/pol_std/HD215806/master_2020-04-27_HD215806_P1-P3R80-160_mac_comb.xlsx")
    high_pol_std_04_27_02 = load_pol_data("./stats/2020-04-27/pol_std/OJ287/master_2020-04-27_OJ287_P1-P3R160-172_mac_comb_corr.xlsx")

    ret_list_target.append(target_data_04_27)
    
    ret_list_zero_pol.append(zero_pol_std_04_27)
    
    ret_list_high_pol.append(high_pol_std_04_27)
    ret_list_high_pol.append(high_pol_std_04_27_02)
        
    ##########
    #2020-04-28
    #The OK batch
    
    target_data_04_28 = load_pol_data("./stats/2020-04-28/target/EE_Cep/master_2020-04-28_eecep_P1-P3R200-270_mac_comb_corr.xlsx")
    
    zero_pol_std_04_28 = load_pol_data("./stats/2020-04-28/pol_std/HD212311/master_2020-04-28_hd212311_P1-P3R60-120_mac_comb.xlsx")

    high_pol_std_04_28 = load_pol_data("./stats/2020-04-28/pol_std/BD64106/master_2020-04-28_bd64106_P1-P3R0-60_mac_comb.xlsx") 
    high_pol_std_04_28_02 = load_pol_data("./stats/2020-04-28/pol_std/HD215806/master_2020-04-28_hd215806_P1-P3R120-180_mac_comb_corr.xlsx")
    high_pol_std_04_28_03 = load_pol_data("./stats/2020-04-28/pol_std/OJ287/master_2020-04-28_oj287_P1-P3R215-243_mac_comb_corr.xlsx")
    
    ret_list_target.append(target_data_04_28)
    
    ret_list_zero_pol.append(zero_pol_std_04_28)
    
    ret_list_high_pol.append(high_pol_std_04_28)
    ret_list_high_pol.append(high_pol_std_04_28_02)
    ret_list_high_pol.append(high_pol_std_04_28_03)
    
    ##########
    #2020-04-30

    target_data_04_30 = load_pol_data("./stats/2020-04-30/target/EE_Cep/master_2020-04-30_eecep_P1-P3R80-110_mac_comb.xlsx")
    
    zero_pol_std_04_30 = load_pol_data("./stats/2020-04-30/pol_std/HD212311/master_2020-04-30_hd212311_P1-P3R0-48_mac_comb.xlsx")

    high_pol_std_04_30 = load_pol_data("./stats/2020-04-30/pol_std/HD215806/master_2020-04-30_hd215806_P1-P3R48-78_mac_comb_corr.xlsx") 

    ret_list_target.append(target_data_04_30)
    
    ret_list_zero_pol.append(zero_pol_std_04_30)
    
    ret_list_high_pol.append(high_pol_std_04_30)
    
    ##########
    #2020-05-01

    target_data_05_01 = load_pol_data("./stats/2020-05-01/target/EE_Cep/master_2020-05-01_eecep_P1-P3R152-240_mac_comb_corr.xlsx")
    
    zero_pol_std_05_01 = load_pol_data("./stats/2020-05-01/pol_std/G191B2B/master_2020-05-01_g191b2b_P1-P3R102-126_mac_comb_corr.xlsx")
    zero_pol_std_05_01_02 = load_pol_data("./stats/2020-05-01/pol_std/HD212311/master_2020-05-01_hd212311_P1-P3R126-172_mac_comb.xlsx")

    high_pol_std_05_01 = load_pol_data("./stats/2020-05-01/pol_std/BD64106/master_2020-05-01_bd64106_P1-P3R0-102_mac_comb_corr.xlsx")
    high_pol_std_05_01_02 = load_pol_data("./stats/2020-05-01/pol_std/HD215806/master_2020-05-01_hd215806_P1-P3R172-270_mac_comb.xlsx") 

    
    ret_list_zero_pol.append(zero_pol_std_05_01)
    ret_list_zero_pol.append(zero_pol_std_05_01_02)
    
    ret_list_high_pol.append(high_pol_std_05_01)
    ret_list_high_pol.append(high_pol_std_05_01_02)
    
    ret_list_target.append(target_data_05_01)
    
    ##########
    #2020-05-04
    
    target_data_05_04 = load_pol_data("./stats/2020-05-04/target/EE_Cep/master_2020-05-04_eecep_P1-P3R140-204_mac_comb.xlsx")
    
    zero_pol_std_05_04 = load_pol_data("./stats/2020-05-04/pol_std/G191B2B/master_2020-05-04_g191b2b_P1-P3R24-44_mac_comb_x.xlsx")
    zero_pol_std_05_04_02 = load_pol_data("./stats/2020-05-04/pol_std/HD212311/master_2020-05-04_hd212311_P1-P3R44-104_mac_comb_corr.xlsx")
    
    
    high_pol_std_05_04 = load_pol_data("./stats/2020-05-04/pol_std/BD64106/master_2020-05-04_bd64106_P1-P3R0-24_mac_comb.xlsx")
    high_pol_std_05_04_02 = load_pol_data("./stats/2020-05-04/pol_std/HD215806/master_2020-05-04_hd215806_P1-P3R104-164_mac_comb_corr.xlsx")
    high_pol_std_05_04_03 = load_pol_data("./stats/2020-05-04/pol_std/OJ287/master_2020-05-04_oj287_P1-P3R243-265_mac_comb.xlsx")
    
    ret_list_high_pol.append(high_pol_std_05_04)
    ret_list_high_pol.append(high_pol_std_05_04_02)
    ret_list_high_pol.append(high_pol_std_05_04_03)
    
    ret_list_zero_pol.append(zero_pol_std_05_04)
    ret_list_zero_pol.append(zero_pol_std_05_04_02)
    
    ret_list_target.append(target_data_05_04)
    
    ##########
    #2020-05-06
    
    target_data_05_06 = load_pol_data("./stats/2020-05-06/target/EE_Cep/master_2020-05-06_eecep_P1-P3R100-168_mac_comb_corr_x.xlsx")
    
    zero_pol_std_05_06 = load_pol_data("./stats/2020-05-06/pol_std/HD212311/master_2020-05-06_hd212311_P1-P3R0-50_mac_comb_corr.xlsx")
    
    high_pol_std_05_06 = load_pol_data("./stats/2020-05-06/pol_std/HD215806/master_2020-05-06_hd215806_P1-P3R50-134_mac_comb_corr.xlsx")
    
    
    ret_list_high_pol.append(high_pol_std_05_06)
    
    ret_list_zero_pol.append(zero_pol_std_05_06)
    
    ret_list_target.append(target_data_05_06)
    
    ##########
    #2020-05-12
    
    target_data_05_12 = load_pol_data("./stats/2020-05-12/target/EE_Cep/master_2020-05-12_eecep_P1-P3R108-154_mac_comb_corr.xlsx")
    
    zero_pol_std_05_12 = load_pol_data("./stats/2020-05-12/pol_std/HD212311/master_2020-05-12_hd212311_P1-P3R86-216_mac_comb_corr.xlsx")
    
    high_pol_std_05_12 = load_pol_data("./stats/2020-05-12/pol_std/BD64106/master_2020-05-12_bd64106_P1-P3R0-85_mac_comb_corr.xlsx")
    high_pol_std_05_12_02 = load_pol_data("./stats/2020-05-12/pol_std/HD215806/master_2020-05-12_hd215806_P1-P3R216-350_mac_comb_corr.xlsx")
    high_pol_std_05_12_03 = load_pol_data("./stats/2020-05-12/pol_std/oj287/master_2020-05-12_oj287_P1-P3R409-429_mac_comb_corr.xlsx")
    
    
    ret_list_high_pol.append(high_pol_std_05_12)
    ret_list_high_pol.append(high_pol_std_05_12_02)
    ret_list_high_pol.append(high_pol_std_05_12_03)
    
    ret_list_zero_pol.append(zero_pol_std_05_12)
    
    ret_list_target.append(target_data_05_12)
    
    ##########
    #2020-05-14
    #Target needs to be corrected
    
    target_data_05_14 = load_pol_data("./stats/2020-05-14/target/EE_Cep/master_2020-05-14_eecep_P1-P3R108-156_mac_comb_corr.xlsx")
    
    zero_pol_std_05_14 = load_pol_data("./stats/2020-05-14/pol_std/HD212311/master_2020-05-14_hd212311_P1-P3R80-252_mac_comb_corr.xlsx")
    
    high_pol_std_05_14 = load_pol_data("./stats/2020-05-14/pol_std/BD64106/master_2020-05-14_bd64106_P1-P3R0-80_mac_comb.xlsx")
    high_pol_std_05_14_02 = load_pol_data("./stats/2020-05-14/pol_std/HD215806/master_2020-05-14_hd215806_P1-P3R252-348_mac_comb_corr.xlsx")
        
    ret_list_target.append(target_data_05_14)
    
    ret_list_high_pol.append(high_pol_std_05_14)
    ret_list_high_pol.append(high_pol_std_05_14_02)
    
    ret_list_zero_pol.append(zero_pol_std_05_14)
    
    ##########
    #2020-05-20
    
    target_data_05_20 = load_pol_data("./stats/2020-05-20/target/EE_Cep/master_2020-05-20_eecep_P1-P3R112-160_mac_comb_corr.xlsx")
    
    ret_list_target.append(target_data_05_20)
    
    ##########
    #2020-06-03
    
    target_data_06_03 = load_pol_data("./stats/2020-06-03/target/EE_Cep/master_2020-06-03_eecep_P1-P3R104-154_mac_comb_corr.xlsx")
    
    ret_list_target.append(target_data_06_03)

    ##########
    #2020-07-29
    #For now let 07_29 be the last batch of data for now
    target_data_07_29 = load_pol_data("./stats/2020-07-29/target/EE_Cep/master_2020-07-29_eecep_P1-P3R360-484_mac_comb.xlsx")
    
    zero_pol_std_07_29 = load_pol_data("./stats/2020-07-29/pol_std/HD212311/master_2020-07-29_hd212311_P1-P3R0-112_mac_comb.xlsx")
    
    high_pol_std_07_29 = load_pol_data( "./stats/2020-07-29/pol_std/HD215806/master_2020-07-29_hd215806_P1-P3R112-202_mac_comb.xlsx")
    #a moment in time i'm not crying you left me on my own
    
    ret_list_target.append(target_data_07_29)
    ret_list_zero_pol.append(zero_pol_std_07_29)
    ret_list_high_pol.append(high_pol_std_07_29)
    
    ##########
    #2020-08-25
    
    target_data_08_25 = load_pol_data("./stats/2020-08-25/target/EE_Cep/master_2020-08-25_eecep_P1-P3R388-534_mac_comb.xlsx")
    
    ret_list_target.append(target_data_08_25)
    
    ##########
    #2020-09-08
    
    target_data_09_08 = load_pol_data("./stats/2020-09-08/target/EE_Cep/master_2020-09-08_eecep_P1-P3R212-332_mac_comb.xlsx")
    
    ret_list_target.append(target_data_09_08)
    
    ##########
    #2020-09-22
    
    target_data_09_22 = load_pol_data("./stats/2020-09-22/target/EE_Cep/master_2020-09-22_eecep_P1-P3R200-317_mac_comb.xlsx")
    
    ret_list_target.append(target_data_09_22)
    
    ##########
    #2020-10-08
    
    target_data_10_08 = load_pol_data("./stats/2020-10-08/target/EE_Cep/master_2020-10-08_eecep_P1-P3R72-124_mac_comb.xlsx")
    
    ret_list_target.append(target_data_10_08)
    
    ##########
    #2020-10-28
    
    target_data_10_28 = load_pol_data("./stats/2020-10-28/target/EE_Cep/master_2020-10-28_eecep_P1-P3R216-352_mac_comb.xlsx")
    
    ret_list_target.append(target_data_10_28)
    
    #2021 = 0.2 phase where NIR shenanigans were thought to happen.
    #######################################################################################
    #######################################################################################
    #######################################################################################
    #2021-04-06
    
    target_data_04_06 = load_pol_data("./stats/2021-04-06/target/EE_Cep/master_2021-04-06_eecep_P1-P3R100-168_mac_comb_corr.xlsx")
    
    ret_list_target.append(target_data_04_06)
    
    ##########
    #2021-05-05
    
    target_data_05_05 = load_pol_data("./stats/2021-05-05/target/EE_Cep/master_2021-05-05_eecep_P1-P3R212-284_mac_comb_corr.xlsx")
    
    ret_list_target.append(target_data_05_05)
    
    ##########
    #2021-05-20
    
    target_data_05_20 = load_pol_data("./stats/2021-05-20/target/EE_Cep/master_2021-05-20_eecep_P1-P3R212-284_mac_comb_corr.xlsx")
    
    ret_list_target.append(target_data_05_20)
    
    ##########
    #2021-05-25
    
    target_data_05_25 = load_pol_data("./stats/2021-05-25/target/EE_Cep/master_2021-05-25_eecep_P1-P3R260-340_mac_comb_corr.xlsx")
    
    ret_list_target.append(target_data_05_25)
    

    return (ret_list_target, ret_list_zero_pol, ret_list_high_pol)
    
def filter_data(pol_data, filter_strings, verb_arg):
    """
    #A look for this funtion
    """
    
    TOI_list = []

    for targs in pol_data:
        for filt_strings in filter_strings:
            if filt_strings in list(targs.keys())[0]:# or '191b2b' in list(targs.keys())[0]:
                TOI_list.append(targs)

    if(verb_arg):
        print("Returned", len(TOI_list), "results")
                
    return(TOI_list)
    
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
                
def outlier_removal(inp_data, q_lim, u_lim):
    print("Specify qlim and ulim to reject scatter outliers")
    
    outlier_rejected_data = cp(inp_data)
    ###
    #Do some outlier rejection using q and u plane
    ###
    
    return outlier_rejected_data

def combine_latest_flates(input_flat_files):
    #print("Flat dir")
    #print(input_flat_files)
    
    p1_filters = []
    p3_filters = []
    
    if input_flat_files == []:
        print("Empty List")
    else:
        for shit in input_flat_files:
            if('p3-r'in shit):
                print(shit)
                path_s = re.search('/flat(.*)', shit)
                print(path_s.group(0))
            elif('p1-r'in shit):
                print(shit)
                path_s = re.search('/flat(.*)', shit)
                print(path_s.group(0))

def my_quadrature_compute(list_of_errors):
    """
    #Takes a list of errors and computes the square root of the sume of squares "quadrature"
    """
    
    squares = []
    for things in list_of_errors:
        squares.append(things**2)

    return(np.sqrt(np.sum(squares)))

def ringo_error_prop( test_list ):
    denom_err = []
    for hi in range(0, len(test_list)):
        de_nom= 1/(test_list[hi]**2)
        denom_err.append(de_nom)

    sigma_q = np.sqrt(1/sum(denom_err))

    return (sigma_q)

def print_list(MJD, f_list):
    print(MJD)
    for i in range(0, len(f_list)):
        print(i, 
        astropy.io.fits.open(f_list[i])[0].header['TIME-OBS'],
        astropy.io.fits.open(f_list[i])[0].header['FILTER'],
        astropy.io.fits.open(f_list[i])[0].header['EXPTIME'],
        astropy.io.fits.open(f_list[i])[0].header['OBJECT'])
        
def data_to_excel(data_struct, obj_name ,who_am_i):
    """
    #It takes in the val array and the the error array
    """
    #print("Translate this data to excel:", who_am_i)
    #print("Data structure:",data_struct)
    #for o in range(0, len(data_struct[0])):
    #    print(data_struct[0][o], data_struct[1][o], data_struct[2][o])
    #write_this = [    ]
    
    workbook = xlsxwriter.Workbook('./tstats/'+ who_am_i  +'.xlsx')
    worksheet = workbook.add_worksheet()
    worksheet.write('A1', 'MJD')
    worksheet.write('B1', 'Adam '+who_am_i+' val')
    worksheet.write('C1', 'Adam '+who_am_i+' error')
    
    #Look through that shit here:
    row = 1
    column = 0
    
    for l in range(0, len(data_struct[0])):
        write_this = [data_struct[0][l], data_struct[1][l], data_struct[2][l]]
    
        for item in write_this:
            worksheet.write(row, column, item)
            column += 1
        
        column = 0
        row += 1 #This row 
    
    workbook.close()