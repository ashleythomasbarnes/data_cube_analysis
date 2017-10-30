import numpy as np
import math
import matplotlib.pyplot as plt
import os
import sys
# sys.path.insert(0, '/Library/Python/2.7/site-packages')
# import scipy
from scipy.optimize import curve_fit
from astropy import units as u
from astropy.table import Column
from astropy.table import Table
from astropy.io import ascii
import errors
import gradients

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size'   : 10} )
plt.rc('text', usetex=True)

#
ntwoh_dir = '/Users/ashleybarnes/Dropbox/work/IRAM-30m/Barnes_kinematics/leodis_output_ash/cloudf/n2h+/final/output_ascii/'
ntwoh_files = os.listdir(ntwoh_dir)
number_comps_ntwoh = len(ntwoh_files)
comps_ntwoh = np.empty([number_comps_ntwoh, 50, 15]) * np.nan
all_fit_count = 0; ind_fit_counts_ntwoh = np.empty([number_comps_ntwoh])
for i in range(number_comps_ntwoh):
    ind_fit_counts_ntwoh[i] = len(np.genfromtxt(ntwoh_dir+ntwoh_files[i]))
    print 'N2H Fil '+str(i)+': ' + str(ind_fit_counts_ntwoh[i])
    comps_ntwoh[i, 0:len(np.genfromtxt(ntwoh_dir+ntwoh_files[i])), :] =  np.genfromtxt(ntwoh_dir+ntwoh_files[i])
    all_fit_count += len(np.genfromtxt(ntwoh_dir+ntwoh_files[i]))
print 'Total number of fits in file: '+str(all_fit_count)
#

#
c18o_dir = '/Users/ashleybarnes/Dropbox/work/IRAM-30m/Barnes_kinematics/leodis_output_ash/cloudf/c18o/final/output_ascii/'
c18o_files = os.listdir(c18o_dir)
number_comps_c18o = len(c18o_files)
comps_c18o = np.empty([number_comps_c18o, 200, 15]) * np.nan
all_fit_count = 0; ind_fit_counts_c18o = np.empty([number_comps_c18o])
for i in range(number_comps_c18o):
    ind_fit_counts_c18o[i] = len(np.genfromtxt(c18o_dir+c18o_files[i]))
    print 'C18O Fil '+str(i)+': ' + str(ind_fit_counts_c18o[i])
    comps_c18o[i, 0:len(np.genfromtxt(c18o_dir+c18o_files[i])), :] =  np.genfromtxt(c18o_dir+c18o_files[i])
    all_fit_count += len(np.genfromtxt(c18o_dir+c18o_files[i]))
print 'Total number of fits in file: '+str(all_fit_count)
#

source_distance_pc = 3.7e3
distance_err = source_distance_pc * 0.15
##

table_arr_ntwoh = []
table_arr_c18o = []

comp_names_n2h = [r'F$_{\rm PPV4}$']
comp_names_c18o = [r'F$_{\rm PPV1}$', r'F$_{\rm PPV2}$', r'F$_{\rm PPV3}$', r'F$_{\rm PPV4}$']

for i in range(number_comps_ntwoh):

    comp_ntwoh = comps_ntwoh[i]

    comp_ntwoh_velo = comp_ntwoh[:,5]
    comp_ntwoh_width = comp_ntwoh[:,7]
    comp_ntwoh_ave_velo = np.nanmean(comp_ntwoh_velo); comp_ntwoh_ave_velo_err = np.nanstd(comp_ntwoh_velo); comp_ntwoh_ave_velo_eom = comp_ntwoh_ave_velo_err / np.sqrt(ind_fit_counts_ntwoh[i])
    comp_ntwoh_width_ave = np.nanmean(comp_ntwoh_width); comp_ntwoh_width_err = np.nanstd(comp_ntwoh_width)

    comp_ntwoh_ave_velo_array = np.empty([len(comp_ntwoh_velo)]) * np.nan; comp_ntwoh_ave_velo_array.fill(comp_ntwoh_ave_velo)
    comp_ntwoh_ra = comp_ntwoh[:,1] / 3600.; comp_ntwoh_ra_radian = np.empty([len(comp_ntwoh_ra)]) * np.nan
    comp_ntwoh_dec = comp_ntwoh[:,2] / 3600.; comp_ntwoh_dec_radian = np.empty([len(comp_ntwoh_dec)]) * np.nan
    for i in range(len(comp_ntwoh_ra)):
        comp_ntwoh_ra_radian[i] = math.radians(comp_ntwoh_ra[i])
        comp_ntwoh_dec_radian[i] = math.radians(comp_ntwoh_dec[i])

    p0 = (1000, 1000)
    nan_ID = ~np.isnan(comp_ntwoh_ra_radian)
    comp_ntwoh_popt, comp_ntwoh_pcov = curve_fit(gradients.velo_gradients, (comp_ntwoh_ave_velo_array[nan_ID], comp_ntwoh_ra_radian[nan_ID], comp_ntwoh_dec_radian[nan_ID]), comp_ntwoh_velo[nan_ID], p0)
    comp_ntwoh_perr = np.sqrt(np.diag(comp_ntwoh_pcov))
    comp_ntwoh_delta_v_lsr, comp_ntwoh_delta_v_lsr_perr = errors.mag_error(comp_ntwoh_popt, comp_ntwoh_perr, source_distance_pc, distance_err)
    comp_ntwoh_theta_delta_v_lsr, comp_ntwoh_theta_delta_v_lsr_perr = errors.ang_error(comp_ntwoh_popt, comp_ntwoh_perr)

    print 'Cloud F N2H+ - Velocity: ', np.round(comp_ntwoh_ave_velo, decimals = 1), ' pm ', np.round(comp_ntwoh_ave_velo_err, decimals = 1), ' , ', \
                        'Width: ', np.round(comp_ntwoh_width_ave, decimals = 1), ' pm ', np.round(comp_ntwoh_width_err, decimals = 1), ' , ', \
                        'Gradient mag: ', np.round(comp_ntwoh_delta_v_lsr, decimals = 1), ' pm ', np.round(comp_ntwoh_delta_v_lsr_perr, decimals = 1), ' , ', \
                        'Gradient angle: ', np.round(math.degrees(comp_ntwoh_theta_delta_v_lsr), decimals = 1), ' pm ', np.round(math.degrees(comp_ntwoh_theta_delta_v_lsr_perr), decimals = 1)

    table_arr_temp = [comp_ntwoh_ave_velo, comp_ntwoh_ave_velo_err, comp_ntwoh_width_ave, comp_ntwoh_width_err, comp_ntwoh_delta_v_lsr, comp_ntwoh_delta_v_lsr_perr, math.degrees(comp_ntwoh_theta_delta_v_lsr), math.degrees(comp_ntwoh_theta_delta_v_lsr_perr)]
    table_arr_ntwoh = table_arr_ntwoh + [table_arr_temp]


for j in range(number_comps_c18o):

    comp_c18o = comps_c18o[j]
    comp_c18o_velo = comp_c18o[:,5]
    comp_c18o_width = comp_c18o[:,7]

    comp_c18o_ave_velo = np.nanmean(comp_c18o_velo); comp_c18o_ave_velo_err = np.nanstd(comp_c18o_velo); comp_c18o_ave_velo_eom = comp_c18o_ave_velo_err / np.sqrt(ind_fit_counts_c18o[j])
    comp_c18o_width_ave = np.nanmean(comp_c18o_width); comp_c18o_width_err = np.nanstd(comp_c18o_width)

    comp_c18o_ave_velo_array = np.empty([len(comp_c18o_velo)]) * np.nan; comp_c18o_ave_velo_array.fill(comp_c18o_ave_velo)
    comp_c18o_ra = comp_c18o[:,1] / 3600.; comp_c18o_ra_radian = np.empty([len(comp_c18o_ra)]) * np.nan
    comp_c18o_dec = comp_c18o[:,2] / 3600.; comp_c18o_dec_radian = np.empty([len(comp_c18o_dec)]) * np.nan
    for i in range(len(comp_c18o_ra)):
        comp_c18o_ra_radian[i] = math.radians(comp_c18o_ra[i])
        comp_c18o_dec_radian[i] = math.radians(comp_c18o_dec[i])
    nan_ID = ~np.isnan(comp_c18o_dec_radian)
    comp_c18o_popt, comp_c18o_pcov = curve_fit(gradients.velo_gradients, (comp_c18o_ave_velo_array[nan_ID], comp_c18o_ra_radian[nan_ID], comp_c18o_dec_radian[nan_ID]), comp_c18o_velo[nan_ID], p0)
    comp_c18o_perr = np.sqrt(np.diag(comp_c18o_pcov))
    comp_c18o_delta_v_lsr, comp_c18o_delta_v_lsr_perr = errors.mag_error(comp_c18o_popt, comp_c18o_perr, source_distance_pc, distance_err)
    comp_c18o_theta_delta_v_lsr, comp_c18o_theta_delta_v_lsr_perr = errors.ang_error(comp_c18o_popt, comp_c18o_perr)

    print 'Cloud F C18O+ - Velocity: ', np.round(comp_c18o_ave_velo, decimals = 1), ' pm ', np.round(comp_c18o_ave_velo_err, decimals = 1), ' , ', \
                        'Width: ', np.round(comp_c18o_width_ave, decimals = 1), ' pm ', np.round(comp_c18o_width_err, decimals = 1), ' , ', \
                        'Gradient mag: ', np.round(comp_c18o_delta_v_lsr, decimals = 1), ' pm ', np.round(comp_c18o_delta_v_lsr_perr, decimals = 1), ' , ', \
                        'Gradient angle: ', np.round(math.degrees(comp_c18o_theta_delta_v_lsr), decimals = 1), ' pm ', np.round(math.degrees(comp_c18o_theta_delta_v_lsr_perr), decimals = 1)


    table_arr_temp = [comp_c18o_ave_velo, comp_c18o_ave_velo_err, comp_c18o_width_ave, comp_c18o_width_err, comp_c18o_delta_v_lsr, comp_c18o_delta_v_lsr_perr, math.degrees(comp_c18o_theta_delta_v_lsr), math.degrees(comp_c18o_theta_delta_v_lsr_perr)]
    table_arr_c18o = table_arr_c18o + [table_arr_temp]


table_arr_c18o = np.asarray(np.around(table_arr_c18o, decimals = 2))
table_arr_ntwoh = np.asarray(np.around(table_arr_ntwoh, decimals = 2))

table_ntwoh = Table(table_arr_ntwoh)
table_c18o = Table(table_arr_c18o)

col_ntwoh = Column(name = 'name', data = comp_names_n2h)
col_c18o = Column(name = 'name', data = comp_names_c18o)

table_ntwoh.add_column(col_ntwoh, index = 0)
table_c18o.add_column(col_c18o, index = 0)

col_c18o = Column(name = 'number',  data = ind_fit_counts_c18o)
table_c18o.add_column(col_c18o, index = 1)

col_ntwoh = Column(name = 'number',  data = ind_fit_counts_ntwoh)
table_ntwoh.add_column(col_ntwoh, index = 1)

print '\n PRINTING LAXTEX TABLES'
ascii.write(table_ntwoh, format='latex')
ascii.write(table_c18o, format='latex')
