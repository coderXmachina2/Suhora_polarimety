import gzip
import shutil
import glob
import astropy
import xlrd
import xlwt
import xlsxwriter
import importlib
import re
import numpy as np

from astropy.time import Time
from xlwt import Workbook
from reduct_funcs import funcs_polarimetry
importlib.reload(funcs_polarimetry)

def combine_excels(excel_1, excel_2, sv_out , MJD, targ_obj, obs_filt ,strt_ind, end_ind):
    """
    Function that takes in two excel scripts and combines them into the usual 1
    
    One last chance to do this bro. Today 23/10/2021. Do it.
    """
    #create new workbookq
    
    print("Sheet 1:",excel_1)
    print("Sheet 2:",excel_2)
    
    filename= targ_obj+"_P1-P3"+obs_filt
    workbook = xlsxwriter.Workbook(sv_out+'master_'+MJD+"_"+ filename+str(strt_ind)+"-" +str(end_ind)+'_mac_comb.xlsx') #This is what. a new workbook?
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
    worksheet.write('AE1', 'PD')
    worksheet.write('AF1', 'PD error')
    worksheet.write('AG1', 'PA')
    worksheet.write('AH1', 'PA error')
    
    wb_1 = xlrd.open_workbook((excel_1))
    wb_2 = xlrd.open_workbook((excel_2))
    
    sheet_1 = wb_1.sheet_by_index(0)
    sheet_2 = wb_2.sheet_by_index(0)
       
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
        
        worksheet.write(u, 28, (sheet_1.cell_value(rowx=u, colx=7)-sheet_1.cell_value(rowx=u, colx=11))/(sheet_1.cell_value(rowx=u, colx=7)+sheet_1.cell_value(rowx=u, colx=11))) #write q, q_err
        worksheet.write(u, 29, p_ans  ) #write q, q_err
    
    #Its truly only saving time if its calculating the fluexs and erros 
    
    workbook.close() #I guess this writes the thing

def data_loader():
    ret_list_zero_pol = []
    ret_list_high_pol = []
    ret_list_target = []
    
    print("Load all excel data") #One string the path is everything before master. the filename is after.

    target_data_03_05 = funcs_polarimetry.load_pol_data( "./stats/2020-03-05/target/EE_Cep/master_2020-03-05_eecep_P1-P3R0-102_combined.xlsx", False )
    zero_pol_std_03_05 = funcs_polarimetry.load_pol_data("./stats/2020-03-05/pol_std/G191B2B/master_2020-03-05_g191b2b_P1-P3R0-46_combined.xlsx", False )
    high_pol_std_03_05 = funcs_polarimetry.load_pol_data("./stats/2020-03-05/pol_std/HD215806/master_2020-03-05_hd215806_P1-P3R46-106_combined.xlsx", False)
    
    ret_list_zero_pol.append(zero_pol_std_03_05)
    ret_list_high_pol.append(high_pol_std_03_05)
    ret_list_target.append(target_data_03_05)
    
    target_data_03_14 = funcs_polarimetry.load_pol_data( "./stats/2020-03-14/target/EE_Cep/master_2020-03-14_EE Cep_P1-P3R0-69_combined.xlsx", False )
    zero_pol_std_03_14 = funcs_polarimetry.load_pol_data( "./stats/2020-03-14/pol_std/G191B2B/master_2020-03-14_G191B2B_P1-P3R0-60_combined.xlsx", False)
    high_pol_std_03_14 = funcs_polarimetry.load_pol_data("./stats/2020-03-14/pol_std/HD215806/master_2020-03-14_HD215806_P1-P3R60-301_combined.xlsx", False)

    ret_list_zero_pol.append(zero_pol_std_03_14)
    ret_list_high_pol.append(high_pol_std_03_14)
    ret_list_target.append(target_data_03_14)
    
    target_data_03_15 = funcs_polarimetry.load_pol_data("./stats/2020-03-15/target/EE_Cep/master_2020-03-15_EE Cep_P1-P3R19-59_combined.xlsx", False )
    zero_pol_std_03_15 = funcs_polarimetry.load_pol_data("./stats/2020-03-15/pol_std/G191B2B/master_2020-03-15_G191B2B_P1-P3R0-40_combined.xlsx", False )
    high_pol_std_03_15 = funcs_polarimetry.load_pol_data("./stats/2020-03-15/pol_std/HD215806/master_2020-03-15_HD215806_P1-P3R40-109_combined.xlsx", False)
    
    ret_list_zero_pol.append(zero_pol_std_03_15)
    ret_list_high_pol.append(high_pol_std_03_15)
    ret_list_target.append(target_data_03_15)

    target_data_03_24 = funcs_polarimetry.load_pol_data("./stats/2020-03-24/target/EE_Cep/master_2020-03-24_eecep_P1-P3R600-770_mac_comb.xlsx", False)
    
    ret_list_target.append(target_data_03_24)
    
    target_data_03_25 = funcs_polarimetry.load_pol_data("./stats/2020-03-25/target/EE_Cep/master_2020-03-25_eecep_P1-P3R464-540_combined.xlsx", False)
    zero_pol_std_03_25 = funcs_polarimetry.load_pol_data("./stats/2020-03-25/pol_std/HD212311/master_2020-03-25_hd212311_P1-P3R0-44_combined.xlsx", False)
    
    ret_list_zero_pol.append(zero_pol_std_03_25)
    ret_list_target.append(target_data_03_25)

    target_data_03_26 = funcs_polarimetry.load_pol_data("./stats/2020-03-26/target/EE_Cep/master_2020-03-26_eecep_P1-P3R40-109_mac_comb.xlsx", False)
    zero_pol_std_03_26 = funcs_polarimetry.load_pol_data("./stats/2020-03-26/pol_std/G191B2B/master_2020-03-26_g191b2b_P1-P3R0-30_combined.xlsx", False )
    zero_pol_std_03_26_02 = funcs_polarimetry.load_pol_data("./stats/2020-03-26/pol_std/HD212311/master_2020-03-26_hd212311_P1-R90-130_combined.xlsx", False)
    high_pol_std_03_26 = funcs_polarimetry.load_pol_data("./stats/2020-03-26/pol_std/oj287/master_2020-03-26_oj287_P1-P3R258-296_combined_man_samp.xlsx",False)
    high_pol_std_03_26_02 = funcs_polarimetry.load_pol_data("./stats/2020-03-26/pol_std/HD204827/master_2020-03-26_hd204827_P1-P3R60-90_combined.xlsx" , False)
    high_pol_std_03_26_03 = funcs_polarimetry.load_pol_data("./stats/2020-03-26/pol_std/HD215806/master_2020-03-26_HD215806_P1-P3R40-109_mac_comb.xlsx", False)
    high_pol_std_03_26_04 = funcs_polarimetry.load_pol_data("./stats/2020-03-26/pol_std/HD251204/master_2020-03-26_hd251204_P1-P3R40-109_mac_comb.xlsx", False)
    
    ret_list_target.append(target_data_03_26)
    ret_list_zero_pol.append(zero_pol_std_03_26)
    ret_list_zero_pol.append(zero_pol_std_03_26_02)
    ret_list_high_pol.append(high_pol_std_03_26)
    ret_list_high_pol.append(high_pol_std_03_26_02)
    ret_list_high_pol.append(high_pol_std_03_26_03)
    ret_list_high_pol.append(high_pol_std_03_26_04)

    target_data_03_31 = funcs_polarimetry.load_pol_data("./stats/2020-03-31/target/EE_Cep/master_2020-03-31_eecep_P1-P3R140-194_mac_comb_x.xlsx", False)
    zero_pol_std_03_31 = funcs_polarimetry.load_pol_data("./stats/2020-03-31/pol_std/HD212311/master_2020-03-31_hd212311_P1-P3R0-40_mac_comb_outlier_zapped.xlsx", False )
    
    ret_list_target.append(target_data_03_31)
    ret_list_zero_pol.append(zero_pol_std_03_31)

    target_data_04_01 = funcs_polarimetry.load_pol_data("./stats/2020-04-01/target/EE_Cep/master_2020-04-01_EE_Cep_P1-P3R307-425.xlsx", False)
    
    ret_list_target.append(target_data_04_01)

    target_data_04_04 = funcs_polarimetry.load_pol_data("./stats/2020-04-04/target/EE_Cep/master_2020-04-04_EE Cep_P3-R26-46_mancombined_corr.xlsx", False)
    
    ret_list_target.append(target_data_04_04)

    target_data_04_05 = funcs_polarimetry.load_pol_data("./stats/2020-04-05/target/EE_Cep/master_2020-04-05_EE Cep_P1-P3R0-43_corr.xlsx", False)
    
    ret_list_target.append(target_data_04_05)

    zero_pol_std_04_06 = funcs_polarimetry.load_pol_data("./stats/2020-04-06/pol_std/HD212311/master_2020-04-06_hd212311_P1-P3R0-180_mac_comb.xlsx", False)
    high_pol_std_04_06 = funcs_polarimetry.load_pol_data("./stats/2020-04-06/pol_std/HD215806/master_2020-04-06_hd215806_P1-P3R180-388_mac_comb.xlsx", False)
    high_pol_std_04_06_02 = funcs_polarimetry.load_pol_data("./stats/2020-04-06/pol_std/oj287/master_2020-04-06_oj287_P1-P3R388-402_mac_comb.xlsx", False)
    target_data_04_06 = funcs_polarimetry.load_pol_data("./stats/2020-04-06/target/EE_Cep/master_2020-04-06_eecep_P1-P3R428-618_mac_comb.xlsx", False)

    ret_list_target.append(target_data_04_06)
    ret_list_zero_pol.append(zero_pol_std_04_06)
    ret_list_high_pol.append(high_pol_std_04_06)
    ret_list_high_pol.append(high_pol_std_04_06_02)
    
    zero_pol_std_04_07 = funcs_polarimetry.load_pol_data("./stats/2020-04-07/pol_std/HD212311/master_2020-04-07_hd212311_P1-P3R0-180_mac_comb.xlsx", False)
    high_pol_std_04_07 = funcs_polarimetry.load_pol_data("./stats/2020-04-07/pol_std/HD215806/master_2020-04-07_hd215806_P1-P3R180-300_mac_comb.xlsx", False)
    high_pol_std_04_07_02 = funcs_polarimetry.load_pol_data("./stats/2020-04-07/pol_std/oj287/master_2020-04-07_oj287_P1-P3R300-312_mac_comb.xlsx", False)
    target_data_04_07 = funcs_polarimetry.load_pol_data("./stats/2020-04-07/target/EE_Cep/master_2020-04-07_eecep_P1-P3R504-712_mac_comb_corr.xlsx", False)
    
    ret_list_target.append(target_data_04_07)
    ret_list_zero_pol.append(zero_pol_std_04_07)
    ret_list_high_pol.append(high_pol_std_04_07)
    ret_list_high_pol.append(high_pol_std_04_07_02)
    
    zero_pol_std_04_08 = funcs_polarimetry.load_pol_data("./stats/2020-04-08/pol_std/HD212311/master_2020-04-08_hd212311_P1-P3R344-538_mac_comb.xlsx", False)
    high_pol_std_04_08 = funcs_polarimetry.load_pol_data("./stats/2020-04-08/pol_std/HD204827/master_2020-04-08_hd204827_P1-P3R42-276_mac_comb.xlsx", False)
    high_pol_std_04_08_02 = funcs_polarimetry.load_pol_data("./stats/2020-04-08/pol_std/HD215806/master_2020-04-08_hd215806_P1-P3R589-741_mac_comb.xlsx", False)
    target_data_04_08 = funcs_polarimetry.load_pol_data("./stats/2020-04-08/target/EE_Cep/master_2020-04-08_eecep_P1-P3R288-508_mac_comb_corr.xlsx", False)
    
    ret_list_target.append(target_data_04_08)
    ret_list_zero_pol.append(zero_pol_std_04_08)
    ret_list_high_pol.append(high_pol_std_04_08)
    ret_list_high_pol.append(high_pol_std_04_08_02)
    
    zero_pol_std_04_15 = funcs_polarimetry.load_pol_data("./stats/2020-04-15/pol_std/HD212311/master_2020-04-15_hd212311_P1-P3R120-340_mac_comb.xlsx", False)
    high_pol_std_04_15 = funcs_polarimetry.load_pol_data("./stats/2020-04-15/pol_std/HD204827/master_2020-04-15_hd204827_P1-P3R0-120_mac_comb.xlsx", False)
    high_pol_std_04_15_02 = funcs_polarimetry.load_pol_data("./stats/2020-04-15/pol_std/HD215806/master_2020-04-15_hd215806_P1-P3R340-442_mac_comb.xlsx", False)
    high_pol_std_04_15_03 = funcs_polarimetry.load_pol_data("./stats/2020-04-15/pol_std/OJ287/master_2020-04-15_oj287_P1-P3R503-541_mac_comb.xlsx", False)
    target_data_04_15 = funcs_polarimetry.load_pol_data("./stats/2020-04-15/target/EE_Cep/master_2020-04-15_eecep_P1-P3R320-426_mac_comb.xlsx", False)
    
    ret_list_target.append(target_data_04_15)
    ret_list_zero_pol.append(zero_pol_std_04_15)
    ret_list_high_pol.append(high_pol_std_04_15)
    ret_list_high_pol.append(high_pol_std_04_15_02)
    ret_list_high_pol.append(high_pol_std_04_15_03)
   
    #That is a duplicate
    zero_pol_std_04_16 = funcs_polarimetry.load_pol_data("./stats/2020-04-16/pol_std/G191B2B/master_2020-04-16_g191b2b_P1-P3R62-134_mac_comb_corr.xlsx", False)
    zero_pol_std_04_16_02 = funcs_polarimetry.load_pol_data("./stats/2020-04-16/pol_std/HD212311/master_2020-04-16_hd212311_P1-P3R194-234_mac_comb.xlsx", False)
    
    high_pol_std_04_16 = funcs_polarimetry.load_pol_data("./stats/2020-04-16/pol_std/HD204827/master_2020-04-16_hd204827_P1-P3R134-193_mac_comb_corr.xlsx", False)
    high_pol_std_04_16_02 = funcs_polarimetry.load_pol_data("./stats/2020-04-16/pol_std/BD64106/master_2020-04-16_bd64106_P1-P3R0-62_mac_comb.xlsx", False)
    high_pol_std_04_16_03 = funcs_polarimetry.load_pol_data("./stats/2020-04-16/pol_std/HD215806/master_2020-04-16_hd215806_P1-P3R234-286_mac_comb.xlsx", False)
    high_pol_std_04_16_04 = funcs_polarimetry.load_pol_data("./stats/2020-04-16/pol_std/OJ287/master_2020-04-16_oj287_P1-P3R324-344_mac_comb.xlsx", False)
    
    target_data_04_16 = funcs_polarimetry.load_pol_data("./stats/2020-04-16/target/EE_Cep/master_2020-04-16_eecep_P1-P3R300-500_mac_comb_corr.xlsx", False)
    #missing a target here ma dude. Please perform photometry.
    #target_data_04_16 = funcs_polarimetry.load_pol_data("somefile", False)
    ret_list_zero_pol.append(zero_pol_std_04_16)
    ret_list_zero_pol.append(zero_pol_std_04_16_02)
    ret_list_high_pol.append(high_pol_std_04_16)
    ret_list_high_pol.append(high_pol_std_04_16_02)
    ret_list_high_pol.append(high_pol_std_04_16_03)
    ret_list_high_pol.append(high_pol_std_04_16_04)
    #append target!
    ret_list_target.append(target_data_04_16)
    
    #Ma dude
    high_pol_std_04_18 = funcs_polarimetry.load_pol_data("./stats/2020-04-18/pol_std/oj287/master_2020-04-18_OJ287_P1-P3R0-20_mac_comb.xlsx", False)
    ret_list_high_pol.append(high_pol_std_04_18)
    
    #target_data_04_15 = funcs_polarimetry.load_pol_data("./stats/2020-04-15/target/EE_Cep/master_2020-04-15_eecep_P1-P3R320-426_mac_comb.xlsx", False)
    
    #For now let 07_29 be the last batch of data for now
    target_data_07_29 = funcs_polarimetry.load_pol_data("./stats/2020-07-29/target/EE_Cep/master_2020-07-29_eecep_P1-P3R360-484_mac_comb.xlsx", False)
    zero_pol_std_07_29 = funcs_polarimetry.load_pol_data( "./stats/2020-07-29/pol_std/HD212311/master_2020-07-29_hd212311_P1-P3R0-112_mac_comb.xlsx", False)
    high_pol_std_07_29 = funcs_polarimetry.load_pol_data( "./stats/2020-07-29/pol_std/HD215806/master_2020-07-29_hd215806_P1-P3R112-202_mac_comb.xlsx", False)
    #a moment in time i'm not crying you left me on my own
    
    ret_list_target.append(target_data_07_29)
    ret_list_zero_pol.append(zero_pol_std_07_29)
    ret_list_high_pol.append(high_pol_std_07_29)

    return (ret_list_target, ret_list_zero_pol, ret_list_high_pol)
    
def filter_data(pol_data, filter_strings):
    TOI_list = []

    for targs in pol_data:
        for filt_strings in filter_strings:
            if filt_strings in list(targs.keys())[0]:# or '191b2b' in list(targs.keys())[0]:
                TOI_list.append(targs)
                
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