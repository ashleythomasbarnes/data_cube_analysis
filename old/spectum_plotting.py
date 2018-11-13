import numpy as np
import os
import math
import aplpy
import pyfits
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import gauss
matplotlib.style.use('classic')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size'   : 10} )
plt.rc('text', usetex=True)

prefix = '/data01/atbarnes/'
# prefix = '/Users/ashleybarnes/'

msd_map = prefix+'Dropbox/work/KT_extinction_maps/J2000_msd/F_J2000_msd.fits'
co_dir =  prefix+'Dropbox/work/IRAM-30m/Barnes_kinematics/leodis_output_ash/cloudf/13co/final/output_ascii/'
co_file_fits =  prefix+'Dropbox/work/IRAM-30m/Barnes_kinematics/cloudf_grs/cloudf_13co10/STAGE_7/final_solns_updated_7_4sig_manual_JDH.round.dat'
raw_co_file =  prefix+'Dropbox/work/IRAM-30m/Barnes_kinematics/cloudf_grs/cloudf_13co10.fits'
co_file_momentdat =  prefix+'Dropbox/work/IRAM-30m/Barnes_kinematics/cloudf_grs/cloudf_13co10/MISC/coords.dat'

list_co_dir = os.listdir(co_dir)
list_co_dir = np.sort(list_co_dir)

count1 = 0
all_fit_count = 0
co_files = list_co_dir

coms_co = np.empty([len(list_co_dir)+1, 1000, 15]) * np.nan

count1 = 0
for file in list_co_dir:
    if file.endswith(".dat"):
        for i in range(np.genfromtxt(str(co_dir)+str(file)).shape[0]):
            coms_co[count1, i, :] = np.genfromtxt(str(co_dir)+str(file))[i]

        print 'Fil '+str(count1)+': '+str(len(np.genfromtxt(str(co_dir)+str(file))) )
        all_fit_count += len(np.genfromtxt(str(co_dir)+str(file)))
        count1 = count1 + 1

print 'Total number of fits in file: '+str(all_fit_count)

fits_co = np.genfromtxt(co_file_fits)
moment_co = np.genfromtxt(co_file_momentdat)

x = np.unique(moment_co[:, 0])
x = x[np.logical_not(np.isnan(x))]
ra = x[~np.isnan(x)]
ra = ra[::-1]

y = np.unique(moment_co[:, 1])
y = y[np.logical_not(np.isnan(y))]
dec = y[~np.isnan(y)]

do_stage1 = False
blanked_fits = [[]]
if do_stage1:
    print 'STAGE 1'
    counts1 = 0
    counts2 = 0

    for m in [3,4,5]:#range(len(co_files)):
        counts3 = 0
        arraylength = sum(~np.isnan(coms_co[m,:,0]),2)
        for k in range(arraylength):
            for i in range(len(fits_co)):

                if ((( (fits_co[i, 1] > (coms_co[m, k, 1] - 1)) and (fits_co[i, 1] < (coms_co[m, k, 1] + 1)) ) and \
                    ( (fits_co[i, 2] > (coms_co[m, k, 2] - 1)) and (fits_co[i, 2] < (coms_co[m, k, 2] + 1)) ) and \
                    ( (fits_co[i, 5] > (coms_co[m, k, 5] - 0.07)) and (fits_co[i, 5] < (coms_co[m, k, 5] + 0.07)) )) == False).all():

            #if ((fits_co[i, 1] in coms_co[m, k, :]) & (fits_co[i,2] in coms_co[m, k, :]) & (fits_co[i,3] in coms_co[m, k, :]) & (fits_co[i,5] in coms_co[m, k, :]) & (fits_co[i,7] in coms_co[m, k, :])) == False:
                    counts1 = counts1 + 1
                else:
                    #print fits_co[i,:]
                    fits_co[i,:] = np.empty(fits_co[i,:].shape) * np.nan
                    counts2 = counts2 + 1
                    counts3 = counts3 + 1
        print counts3

    print 'Fits blanked '+str(counts2)
    print 'Remaining Fits ' +str(counts1)
    np.savetxt('blanked_arr.txt', fits_co)
else:
    print 'LOADING STAGE 1'
    fits_co = np.genfromtxt('blanked_arr.txt')


ticks_msd = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
fig = plt.figure(figsize=(13,13))
gc = aplpy.FITSFigure(msd_map, figure=fig)
gc.show_colorscale(vmin = 0.0, vmax = 0.5, interpolation='nearest', cmap='Greys')
gc.ticks.set_xspacing(0.05)
gc.ticks.set_yspacing(0.05)
gc.tick_labels.set_xformat('hh:mm:ss'); gc.tick_labels.set_yformat('dd:mm:ss')
gc.add_colorbar()
gc.colorbar.set_width(0.2)
gc.colorbar.set_location('right')
gc.colorbar.set_axis_label_text(r'$\Sigma$ (g cm$^{-2}$)')
gc.recenter(283.3094358, 1.49597906 , width=0.0858175, height=0.3632936)
gc.colorbar.set_ticks(ticks_msd)
gc.ticks.set_color('white')
gc.ticks.hide_x()
gc.ticks.hide_y()

data_co, header_co = pyfits.getdata(raw_co_file, 0, header=True)
data_co = np.squeeze(data_co)
naxis3 = header_co['NAXIS3']
crpix3 = header_co['CRPIX3']
cdelt3 = np.round(header_co['CDELT3'])
crval3 = header_co['CRVAL3']
velo_axis_co = ((np.linspace(1,naxis3,naxis3)-crpix3)*cdelt3+crval3) / 1000.0

grid_Fppv = np.empty([ len(ra) , len(dec) ]) * np.nan

colors = np.array(['green','b','r','y','purple','magenta','cyan','dodgerblue', 'orange',
'r','orange','g','y','red','magenta','cyan','dodgerblue', 'blue',
'r','orange','g','y','purple','magenta','cyan','dodgerblue', 'blue'
,'r','orange','g','y','purple','magenta','cyan','dodgerblue', 'blue'])

lw = 0.3
lw1 = 0.2
size = 0.013542
alpha = 0.5
counter = 0

fit_count = 0

print 'STAGE 2'

for i in range(header_co['NAXIS1']):
    for j in range(header_co['NAXIS2']):

        # if i > 0:
        #     continue

        # if ((i != 0 & j != 0) or (i != 13 & j != 58)):
        #     continue
        # print i, j

        counter = counter + 1
        ax = fig.add_subplot(header_co['NAXIS2'], header_co['NAXIS1'], counter)
        ax.plot(velo_axis_co, data_co[:,j,i], c='black', drawstyle='steps-mid', linewidth=lw1, alpha = alpha)

        [spine.set_linewidth(0.3) for spine in ax.spines.itervalues()]

        for n in range(len(fits_co)):
            if ( (fits_co[n, 1] > (ra[i] - 2)) and (fits_co[n, 1] < (ra[i] + 2)) ) and ( (fits_co[n, 2] > (dec[j] - 2)) and (fits_co[n, 2] < (dec[j] + 2)) ):
                coeff = np.array([fits_co[n, 3], fits_co[n, 5], fits_co[n, 7]/2.35482])
                gauss.gauss_fit = gauss.gauss(velo_axis_co, *coeff)
                ax.plot(velo_axis_co, gauss.gauss_fit, c = 'grey', linewidth=lw, alpha = alpha)

        for m in range(len(co_files)):
            arraylength = sum(~np.isnan(coms_co[m,:,0]),2)
            for k in range(arraylength):
                if ( (coms_co[m, k, 1] > (ra[i] - 2)) and (coms_co[m, k, 1] < (ra[i] + 3)) ) and ( (coms_co[m, k, 2] > (dec[j] - 3)) and (coms_co[m, k, 2] < (dec[j] + 3)) ):
                    fit_count += 1
                    coeff = np.array([coms_co[m, k, 3], coms_co[m, k, 5], coms_co[m, k, 7]/2.35482])
                    coms_co[m, k, 1], ra[i], coms_co[m, k, 2], dec[j]
                    gauss.gauss_fit = gauss.gauss(velo_axis_co, *coeff)
                    ax.plot(velo_axis_co, gauss.gauss_fit, c = colors[m], linewidth=lw, alpha = alpha)

        ax.set_position([0.4081 + (i* size), 0.10075 + (j* size), size, size])
        ax.set_xlim(45, 70)
        ax.set_ylim(-0.5, 5)
        plt.tick_params(axis='both', which='both', bottom='off', top='off', left = 'off', right = 'off', labelleft = 'off', labelbottom='off', direction = 'in')
        ax.patch.set_facecolor('none')
        ax.patch.set_alpha(0)

fig.axes[0].tick_params(color = 'black', direction = 'in', which = 'both', right = True, top = True)
fig.axes[1].tick_params(color = 'black', direction = 'in', which = 'both',  bottom = True, top = True, labeltop = True, labelbottom = False)

print 'For figure plotting: '
size = fig.axes[0].get_position()
print size

fig.set_rasterized(False); plt.savefig('final_figures/big_spec_plot_grs_cloudf.fullres.pdf')
plt.close
