import numpy as np
import math
import matplotlib.pyplot as plt
import os
import sys
sys.path.insert(0, '/Library/Python/2.7/site-packages')
import scipy
from scipy.optimize import curve_fit
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size'   : 10} )
plt.rc('text', usetex=True)


ntwoh_dir = '/Users/ashleybarnes/Dropbox/work/IRAM-30m/Barnes_kinematics/leodis_output/cloudf/n2h+/output_ascii/'
ntwoh_files = ['leodis_nonhierarchical_tree_0.dat', 'leodis_nonhierarchical_tree_1.dat']
c18o_dir = '/Users/ashleybarnes/Documents/work/IRAM-30m/cloudf/line_linking/VESPAOdp20120727_eff_freq_baseline_smoothed_header_cloudf_c18o10/'
c18o_files = ['component_1.dat', 'component_2.dat', 'component_3.dat', 'component_4.dat', 'component_5.dat']
##

comp1_ntwoh = np.genfromtxt(str(ntwoh_dir)+str(ntwoh_files[0])); 
comp2_ntwoh = np.genfromtxt(str(ntwoh_dir)+str(ntwoh_files[1]))
comp1_c18o = np.genfromtxt(str(c18o_dir)+str(c18o_files[0])); 
comp2_c18o = np.genfromtxt(str(c18o_dir)+str(c18o_files[1])); 
comp3_c18o = np.genfromtxt(str(c18o_dir)+str(c18o_files[2])); 
comp4_c18o = np.genfromtxt(str(c18o_dir)+str(c18o_files[3])); 
comp5_c18o = np.genfromtxt(str(c18o_dir)+str(c18o_files[4]))

source_distance_pc = 3.7e3 

def velo_gradients_1d(X, b):
    ave_velo, dec = X
    v_lsr = ave_velo + ( b * dec )
    return v_lsr

def velo_gradients_number_1d(b , distance):
    delta_v_lsr = b / distance
    return delta_v_lsr

def velo_gradients(X, a, b):
    ave_velo, ra, dec = X
    v_lsr = ave_velo + ( a * ra ) + ( b * dec )
    return v_lsr
    
def velo_gradients_number(a, b, distance):
    delta_v_lsr = (a**2 + b**2)**0.5 / distance
    theta_delta_v_lsr = math.atan(a / b)
    return delta_v_lsr, theta_delta_v_lsr

# Velocities  
comp1_ntwoh_velo = comp1_ntwoh[:,4] 
comp2_ntwoh_velo = comp2_ntwoh[:,4] 
comp1_c18o_velo = comp1_c18o[:,3]
comp2_c18o_velo = comp2_c18o[:,3]
comp3_c18o_velo = comp3_c18o[:,3]
comp4_c18o_velo = comp4_c18o[:,3]
comp5_c18o_velo = comp5_c18o[:,3]  
    
# Average velocities    
comp1_ntwoh_ave_velo = np.mean(comp1_ntwoh_velo); comp1_ntwoh_ave_velo_err = np.std(comp1_ntwoh_velo)
comp2_ntwoh_ave_velo = np.mean(comp2_ntwoh_velo); comp2_ntwoh_ave_velo_err = np.std(comp2_ntwoh_velo)
comp1_c18o_ave_velo = np.mean(comp1_c18o_velo); comp1_c18o_ave_velo_err = np.std(comp1_c18o_velo)
comp2_c18o_ave_velo = np.mean(comp2_c18o_velo); comp2_c18o_ave_velo_err = np.std(comp2_c18o_velo)
comp3_c18o_ave_velo = np.mean(comp3_c18o_velo); comp3_c18o_ave_velo_err = np.std(comp3_c18o_velo)
comp4_c18o_ave_velo = np.mean(comp4_c18o_velo); comp4_c18o_ave_velo_err = np.std(comp4_c18o_velo)
comp5_c18o_ave_velo = np.mean(comp5_c18o_velo); comp5_c18o_ave_velo_err = np.std(comp5_c18o_velo)



comp2_ntwoh_ave_velo_array = np.empty([len(comp2_ntwoh_velo)]) * np.nan; comp2_ntwoh_ave_velo_array.fill(comp2_ntwoh_ave_velo)
comp1_c18o_ave_velo_array = np.empty([len(comp1_c18o_velo)]) * np.nan; comp1_c18o_ave_velo_array.fill(comp1_c18o_ave_velo)

# RA/DEC offsets in radians
comp2_ntwoh_ra = comp2_ntwoh[:,1] / 3600.; comp2_ntwoh_ra_radian = np.empty([len(comp2_ntwoh_ra)]) * np.nan 
comp1_c18o_ra = comp1_c18o[:,0] / 3600.; comp1_c18o_ra_radian = np.empty([len(comp1_c18o_ra)]) * np.nan
comp2_ntwoh_dec = comp2_ntwoh[:,2] / 3600.; comp2_ntwoh_dec_radian = np.empty([len(comp2_ntwoh_dec)]) * np.nan 
comp1_c18o_dec = comp1_c18o[:,1] / 3600.; comp1_c18o_dec_radian = np.empty([len(comp1_c18o_dec)]) * np.nan 

for i in range(len(comp2_ntwoh_ra)): 
    comp2_ntwoh_ra_radian[i] = math.radians(comp2_ntwoh_ra[i])
    comp2_ntwoh_dec_radian[i] = math.radians(comp2_ntwoh_dec[i])

for i in range(len(comp1_c18o_ra)): 
    comp1_c18o_ra_radian[i] = math.radians(comp1_c18o_ra[i])
    comp1_c18o_dec_radian[i] = math.radians(comp1_c18o_dec[i])

##### 1D 
p0 = 1
comp2_ntwoh_popt, comp2_ntwoh_pcov = curve_fit(velo_gradients_1d, (comp2_ntwoh_ave_velo_array, comp2_ntwoh_dec_radian), comp2_ntwoh_velo, p0)
comp1_c18o_popt, comp1_c18o_pcov = curve_fit(velo_gradients_1d, (comp1_c18o_ave_velo_array, comp1_c18o_dec_radian), comp1_c18o_velo, p0)

comp2_ntwoh_perr = np.sqrt(np.diag(comp2_ntwoh_pcov)); comp1_c18o_perr = np.sqrt(np.diag(comp1_c18o_pcov))

comp2_ntwoh_delta_v_lsr = velo_gradients_number_1d(comp2_ntwoh_popt, source_distance_pc)
comp1_c18o_delta_v_ls = velo_gradients_number_1d(comp1_c18o_popt, source_distance_pc)

comp2_ntwoh_delta_v_lsr_perr = velo_gradients_number_1d(comp2_ntwoh_perr, source_distance_pc)
comp1_c18o_delta_v_ls_perr = velo_gradients_number_1d(comp1_c18o_perr, source_distance_pc)

print 'Cloud F 1D N2H+: ', comp2_ntwoh_delta_v_lsr, ' pm ', comp2_ntwoh_delta_v_lsr_perr
print 'Cloud F 1D C18O: ', comp1_c18o_delta_v_ls, ' pm ', comp1_c18o_delta_v_ls_perr

##### 3D
p0 = (1, 1)
comp2_ntwoh_popt, comp2_ntwoh_pcov = curve_fit(velo_gradients, (comp2_ntwoh_ave_velo_array, comp2_ntwoh_ra_radian, comp2_ntwoh_dec_radian), comp2_ntwoh_velo, p0)
comp1_c18o_popt, comp1_c18o_pcov = curve_fit(velo_gradients, (comp1_c18o_ave_velo_array, comp1_c18o_ra_radian, comp1_c18o_dec_radian), comp1_c18o_velo, p0)

comp2_ntwoh_perr = np.sqrt(np.diag(comp2_ntwoh_pcov)); comp1_c18o_perr = np.sqrt(np.diag(comp1_c18o_pcov))

comp2_ntwoh_delta_v_lsr, comp2_ntwoh_theta_delta_v_lsr = velo_gradients_number(comp2_ntwoh_popt[0], comp2_ntwoh_popt[1], source_distance_pc)
comp1_c18o_delta_v_ls,  comp1_c18o_theta_delta_v_lsr = velo_gradients_number(comp1_c18o_popt[0], comp1_c18o_popt[1], source_distance_pc)

comp2_ntwoh_delta_v_lsr_perr, comp2_ntwoh_theta_delta_v_lsr_perr = velo_gradients_number(comp2_ntwoh_perr[0], comp2_ntwoh_perr[1], source_distance_pc)
comp1_c18o_delta_v_lsr_perr,  comp1_c18o_theta_delta_v_lsr_perr = velo_gradients_number(comp1_c18o_perr[0], comp1_c18o_perr[1], source_distance_pc)

print 'Cloud H N2H+: ', comp2_ntwoh_delta_v_lsr, ' pm ', comp2_ntwoh_delta_v_lsr_perr, ' , ', math.degrees(comp2_ntwoh_theta_delta_v_lsr), ' pm ', comp2_ntwoh_theta_delta_v_lsr_perr
print 'Cloud H C18O: ', comp1_c18o_delta_v_ls, ' pm ', comp1_c18o_delta_v_lsr_perr, ' , ', math.degrees(comp1_c18o_theta_delta_v_lsr), ' pm ', comp1_c18o_theta_delta_v_lsr_perr



























