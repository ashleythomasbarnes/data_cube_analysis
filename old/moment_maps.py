import numpy as np
import os
import math
import matplotlib as plt
import matplotlib.patheffects as path_effects
# import matplotlib
# matplotlib.style.use('classic')

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size'   : 10} )
plt.rc('text', usetex=True)
#-----#
do_loads = True
do_casa_moments = False
do_aplpy_moments = True
#-----#

if do_loads:
	moments_dir = '/Users/ashleybarnes/Dropbox/work/IRAM-30m/cloudf/moment_maps/'
	maps_dir = '/Users/ashleybarnes/Dropbox/work/IRAM-30m/cloudf/line_fitting/'
	maps_files = ['VESPAOdp20120727_VESPAOdp20120728_combined_eff_baseline_iso_smoothed_header_cloudf_ntwoh10.fits', 'VESPAOdp20120727_eff_freq_baseline_smoothed_header_cloudf_c18o10.fits', 'VESPAOdp20120727_eff_freq_baseline_smoothed_specsmoothed4_header_cloudf_c18o21.fits']
	msd_map = '/Users/ashleybarnes/Dropbox/work/KT_extinction_maps/J2000_msd/F_J2000_msd.fits'
	moments_dir_fits = str(moments_dir)+'fits/'
	moments_dir_fits_n2h = [str(moments_dir_fits)+'VESPAOdp20120727_VESPAOdp20120728_combined_eff_baseline_iso_smoothed_header_cloudf_ntwoh10.fits.image.integrated.fits', str(moments_dir_fits)+'VESPAOdp20120727_VESPAOdp20120728_combined_eff_baseline_iso_smoothed_header_cloudf_ntwoh10.fits.image.weighted_coord.fits', str(moments_dir_fits)+'VESPAOdp20120727_VESPAOdp20120728_combined_eff_baseline_iso_smoothed_header_cloudf_ntwoh10.fits.image.weighted_dispersion_coord.fits']
	moments_dir_fits_c18o = [str(moments_dir_fits)+'VESPAOdp20120727_eff_freq_baseline_smoothed_header_cloudf_c18o10.fits.image.integrated.fits', str(moments_dir_fits)+'VESPAOdp20120727_eff_freq_baseline_smoothed_header_cloudf_c18o10.fits.image.weighted_coord.fits', str(moments_dir_fits)+'VESPAOdp20120727_eff_freq_baseline_smoothed_header_cloudf_c18o10.fits.image.weighted_dispersion_coord.fits']

	mom_chan_range = ['237~332', '258 ~ 388', '86 ~ 152']
	mom_chan_num = np.array([95., 130., 66.])
	rms = np.array([0.045833333, 0.065461538, 0.185969697])
	delta_v = np.array([0.063, 0.054, 0.1067])
	sigma = (mom_chan_num)**0.5 * delta_v * rms
	sigma_threshold = 3.

if do_casa_moments:

	i=2									#file counter
	importfits(fitsimage=str(maps_dir)+str(maps_files[i]), imagename=str(maps_dir)+str(maps_files[i])+'.image', overwrite=True)
	print '\n### Creating moment maps for file: '+str(maps_files[i])+', in channel range: '+str(mom_chan_range[i])+'###\n'
	immoments(imagename = str(maps_dir)+str(maps_files[i])+'.image', moments = [0,1,2,8], chans = mom_chan_range[i], outfile = str(moments_dir)+str(maps_files[i])+'.image')


	for file in os.listdir(moments_dir): 	#masking loop
		mom_files = [str(maps_files[i])+'.image.integrated', str(maps_files[i])+'.image.maximum', str(maps_files[i])+'.image.weighted_coord', str(maps_files[i])+'.image.weighted_dispersion_coord']
		mom_files_withdir = [str(moments_dir)+str(maps_files[i])+'.image.integrated', str(moments_dir)+str(maps_files[i])+'.image.maximum', str(moments_dir)+str(maps_files[i])+'.image.weighted_coord', str(moments_dir)+str(maps_files[i])+'.image.weighted_dispersion_coord']

		if file == mom_files[0]:
			ia.open(mom_files[0])
			ia.calcmask(mask = str(mom_files[0])+' > '+str(sigma[i] * sigma_threshold), name ='masked_int')
			ia.maskhandler('set','masked_int')
			exportfits(imagename=str(mom_files[0]), fitsimage=str(mom_files[0])+'.fits', dropdeg=True, overwrite=True)
			for j in range(1,4):
				ia.open(str(mom_files[j]))
				ia.maskhandler('copy', [str(mom_files[0])+':masked_int', 'masked_int'])
				ia.maskhandler('set','masked_int')
				exportfits(imagename=str(mom_files[j]), fitsimage=str(mom_files[j])+'.fits', dropdeg=True, overwrite=True)

	# for file in os.listdir(moments_dir):
	#     if file.startswith("VESPAOdp2012072"):
	#         exportfits(imagename=str(moments_dir)+str(file), fitsimage=str(moments_dir)+str(file)+'.fits', dropdeg=True, overwrite=True)

if do_aplpy_moments:
	import aplpy
	from matplotlib import pyplot as plt
	from astropy.io import fits
	from astropy.wcs import WCS
	from astropy.utils.data import download_file
	from astropy import units as u
	from astropy.coordinates import SkyCoord
	import pyfits
	import astropy
	from astropy.io import ascii

	## Butler and Tan core loading
	bt12_cores = np.genfromtxt('bt12_cores_radec.txt')
	bt12_cores_names = ['F1', 'F2', 'F3', 'F4']
	file_name_rathborne = 'r_06_cores.txt'
	table = ascii.read(file_name_rathborne, delimiter = ';', guess = False)
	core_ID = np.where(table['MSXDC'] == 'G034.43+00.24')
	table_cores = table[core_ID]
	##

	## Figure paramters
	#v_min_0 = -0.1; v_max_0 = 7
	v_min_0 = 0; v_max_0 = np.array([6, 9])
	v_min_1 = 56; v_max_1 = 59
	v_min_2 = 0.2; v_max_2 = 1.7
	ntwo_ii_levs = np.array([sigma[0]*10., sigma[0]*20, sigma[0]*30, sigma[0]*40, sigma[0]*50, sigma[0]*70, sigma[0]*100])
	c18o_ii_levs = np.array([sigma[0]*120, sigma[0]*160, sigma[0]*180, sigma[0]*200, sigma[0]*220])
	ticks_msd = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.7])
	#ticks_ii = np.array([-0.1, 1.5, 3, 5, 7])
	ticks_ii = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
	ticks_vf = np.array([56, 57, 58 , 59])
	ticks_lw = np.array([0.2, 1, 1.7])
	lw = 0.5
	lw1 = 0.8
	alpha = 0.7
	##

	fig = plt.figure(figsize=(10, 8))
	gc = aplpy.FITSFigure(msd_map, figure=fig, subplot = (2,4,1))
	gc.show_colorscale(vmin = 0.0, vmax = 0.3, interpolation='nearest', cmap='Greys')
	gc.ticks.set_xspacing(0.02)
	gc.ticks.set_yspacing(0.025)
	gc.add_scalebar(1./60.); gc.scalebar.show(1./60., corner = 'top left', alpha=0.9, color='black')
	gc.scalebar.set_label('1 pc'); gc.scalebar.set_font_size(10)
	gc.tick_labels.set_xformat('hh:mm:ss'); gc.tick_labels.set_yformat('dd:mm:ss')
	gc.add_colorbar()
	gc.colorbar.set_width(0.2)
	#gc.colorbar.set_location('top')
	gc.colorbar.set_axis_label_text(r'$\Sigma$ (g cm$^{-2}$)')
	gc.recenter(283.32748, 1.4562265, width=0.031200702, height=0.073794918)
	gc.show_contour(moments_dir_fits_n2h[0], levels= ntwo_ii_levs, linewidths=lw1, alpha = alpha, colors='red')

	pos = np.genfromtxt('pos.txt').T
	path_effects2 = [path_effects.Stroke(linewidth=0.5, foreground='black'), path_effects.Normal()]
	gc.show_markers(bt12_cores[:,0], bt12_cores[:,1], s = (bt12_cores[:,2] / np.nanmax(bt12_cores[:,2])) * 120 , marker='+', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	for i in range(len(bt12_cores)):
		gc.add_label(bt12_cores[i,0]+(10./3600.), bt12_cores[i,1]+(6./3600.), bt12_cores_names[i], color='black', size = 13, path_effects=path_effects2)
	mass_temp = np.array(table_cores['Mass'], dtype = 'float')
	gc.show_markers(np.array(table_cores['_RAJ2000']), np.array(table_cores['_DEJ2000']), s = (mass_temp / np.nanmax(mass_temp)) * 100 , marker='x', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)

	for i in range(len(table_cores)):
		if table_cores['m_MSXDC'][i] == 'MM2' or table_cores['m_MSXDC'][i] == 'MM4':
			continue
		gc.add_label(table_cores['_RAJ2000'][i]+(pos[0, i]/3600.), table_cores['_DEJ2000'][i]+(pos[1, i]/3600.), table_cores['m_MSXDC'][i], color='black', size = 13, path_effects=path_effects2)
	#
	gc.colorbar.set_ticks(ticks_msd)
	gc.ticks.set_color('white')
	gc.hide_xaxis_label()

	color_map = ['Reds', 'Blues']
	gc = aplpy.FITSFigure(moments_dir_fits_n2h[0], figure=fig, subplot = (2,4,2))
	# gc.show_colorscale(pmin= v_min_0, vmax= v_max_0, interpolation='nearest', cmap=color_map[0])
	gc.show_colorscale(vmin = 0, vmax = v_max_0[0], interpolation='nearest', cmap=color_map[0])
	gc.add_beam()
	gc.beam.set_color('black')
	gc.beam.set_hatch('///')
	gc.beam.set_edgecolor('black')
	gc.beam.set_facecolor('none')
	gc.ticks.set_xspacing(0.02)
	gc.ticks.set_yspacing(0.025)
	gc.tick_labels.set_xformat('hh:mm:ss'); gc.tick_labels.set_yformat('dd:mm:ss')
	gc.add_colorbar()
	gc.colorbar.set_width(0.2)
	#gc.colorbar.set_location('top')
	gc.add_label(283.32731, 1.4889172, 'N$_2$H$^{+}$ (1$-$0)', color='black', size = 13)
	gc.colorbar.set_axis_label_text(r'Intergrated intensity (K km\,s$^{-1}$)')
	gc.show_contour(moments_dir_fits_n2h[0], colors='black', levels= ntwo_ii_levs, linewidths=lw, alpha = alpha)
	gc.show_markers(bt12_cores[:,0], bt12_cores[:,1], s = (bt12_cores[:,2] / np.nanmax(bt12_cores[:,2])) * 120 , marker='+', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.show_markers(np.array(table_cores['_RAJ2000']), np.array(table_cores['_DEJ2000']), s = (mass_temp / np.nanmax(mass_temp)) * 100 , marker='x', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.hide_xaxis_label()
	gc.hide_yaxis_label()
	gc.ticks.set_color('white')
	gc.colorbar.set_ticks(ticks_ii)


	path_effects2 = [path_effects.Stroke(linewidth=0.5, foreground='black'), path_effects.Normal()]
	gc.show_markers(bt12_cores[:,0], bt12_cores[:,1], s = (bt12_cores[:,2] / np.nanmax(bt12_cores[:,2])) * 120 , marker='+', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	for i in range(len(bt12_cores)):
		gc.add_label(bt12_cores[i,0]+(10./3600.), bt12_cores[i,1]+(6./3600.), bt12_cores_names[i], color='black', size = 13, path_effects=path_effects2)
	mass_temp = np.array(table_cores['Mass'], dtype = 'float')
	gc.show_markers(np.array(table_cores['_RAJ2000']), np.array(table_cores['_DEJ2000']), s = (mass_temp / np.nanmax(mass_temp)) * 100 , marker='x', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	for i in range(len(table_cores)):
		if table_cores['m_MSXDC'][i] == 'MM2' or table_cores['m_MSXDC'][i] == 'MM4':
			continue
		gc.add_label(table_cores['_RAJ2000'][i]+(pos[0, i]/3600.), table_cores['_DEJ2000'][i]+(pos[1, i]/3600.), table_cores['m_MSXDC'][i], color='black', size = 13, path_effects=path_effects2)
	#

	gc = aplpy.FITSFigure(moments_dir_fits_n2h[1], figure=fig, subplot = (2,4,3))
	gc.show_colorscale(vmin= v_min_1, vmax= v_max_1, interpolation='nearest', cmap=color_map[0])
	gc.add_colorbar()
	gc.colorbar.set_width(0.2)
	#gc.colorbar.set_location('top')
	gc.colorbar.set_axis_label_text(r'Intensity weighted velocity (km\,s$^{-1}$)')
	gc.ticks.set_xspacing(0.02)
	gc.ticks.set_yspacing(0.025)
	gc.tick_labels.set_xformat('hh:mm:ss'); gc.tick_labels.set_yformat('dd:mm:ss')
	gc.show_contour(moments_dir_fits_n2h[0], colors='black', levels= ntwo_ii_levs, linewidths=lw, alpha = alpha)
	gc.show_markers(bt12_cores[:,0], bt12_cores[:,1], s = (bt12_cores[:,2] / np.nanmax(bt12_cores[:,2])) * 120 , marker='+', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.show_markers(np.array(table_cores['_RAJ2000']), np.array(table_cores['_DEJ2000']), s = (mass_temp / np.nanmax(mass_temp)) * 100 , marker='x', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.hide_xaxis_label()
	gc.hide_yaxis_label()
	gc.ticks.set_color('white')
	gc.colorbar.set_ticks(ticks_vf)

	gc = aplpy.FITSFigure(moments_dir_fits_n2h[2], figure=fig, subplot = (2,4,4))
	gc.show_colorscale(vmin= v_min_2, vmax= v_max_2, interpolation='nearest', cmap=color_map[0])
	gc.add_colorbar()
	gc.colorbar.set_width(0.2)
	#gc.colorbar.set_location('top')
	gc.colorbar.set_axis_label_text(r'Intensity weighted width (km\,s$^{-1}$)')
	gc.ticks.set_xspacing(0.02)
	gc.ticks.set_yspacing(0.025)
	gc.tick_labels.set_xformat('hh:mm:ss'); gc.tick_labels.set_yformat('dd:mm:ss')
	gc.show_contour(moments_dir_fits_n2h[0], colors='black', levels= ntwo_ii_levs, linewidths=lw, alpha = alpha)
	gc.show_markers(bt12_cores[:,0], bt12_cores[:,1], s = (bt12_cores[:,2] / np.nanmax(bt12_cores[:,2])) * 120 , marker='+', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.show_markers(np.array(table_cores['_RAJ2000']), np.array(table_cores['_DEJ2000']), s = (mass_temp / np.nanmax(mass_temp)) * 100 , marker='x', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.hide_xaxis_label()
	gc.hide_yaxis_label()
	gc.ticks.set_color('white')
	gc.colorbar.set_ticks(ticks_lw)

	gc = aplpy.FITSFigure(msd_map, figure=fig, subplot = (2,4,5))
	gc.show_colorscale(vmin = 0.0, vmax = 0.3, interpolation='nearest', cmap='Greys')
	gc.ticks.set_xspacing(0.02)
	gc.ticks.set_yspacing(0.025)
	gc.tick_labels.set_xformat('hh:mm:ss'); gc.tick_labels.set_yformat('dd:mm:ss')
	gc.add_colorbar()
	gc.colorbar.set_width(0.2)
	#gc.colorbar.set_location('top')
	gc.colorbar.set_axis_label_text(r'$\Sigma$ (g cm$^{-2}$)')
	gc.recenter(283.32748, 1.4562265, width=0.031200702, height=0.073794918)
	gc.show_contour(moments_dir_fits_c18o[0], levels= c18o_ii_levs, linewidths=lw1, alpha = alpha, colors='blue')
	gc.show_markers(bt12_cores[:,0], bt12_cores[:,1], s = (bt12_cores[:,2] / np.nanmax(bt12_cores[:,2])) * 120 , marker='+', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.show_markers(np.array(table_cores['_RAJ2000']), np.array(table_cores['_DEJ2000']), s = (mass_temp / np.nanmax(mass_temp)) * 100 , marker='x', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.ticks.set_color('white')
	gc.colorbar.set_ticks(ticks_msd)

	gc = aplpy.FITSFigure(moments_dir_fits_c18o[0], figure=fig, subplot = (2,4,6))
	gc.show_colorscale(vmin = 0, vmax = v_max_0[1], interpolation='nearest', cmap=color_map[1])
	gc.ticks.set_xspacing(0.02)
	gc.ticks.set_yspacing(0.025)
	gc.tick_labels.set_xformat('hh:mm:ss'); gc.tick_labels.set_yformat('dd:mm:ss')
	gc.add_colorbar()
	gc.colorbar.set_width(0.2)
	#gc.colorbar.set_location('top')
	gc.colorbar.set_axis_label_text(r'Intergrated intensity (K km\,s$^{-1}$)')
	gc.show_contour(moments_dir_fits_c18o[0], colors='black', levels= c18o_ii_levs, linewidths=lw, alpha = alpha)
	gc.show_markers(bt12_cores[:,0], bt12_cores[:,1], s = (bt12_cores[:,2] / np.nanmax(bt12_cores[:,2])) * 120 , marker='+', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.show_markers(np.array(table_cores['_RAJ2000']), np.array(table_cores['_DEJ2000']), s = (mass_temp / np.nanmax(mass_temp)) * 100 , marker='x', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.colorbar.set_ticks(ticks_ii)
	gc.ticks.set_color('white')
	gc.add_label(283.32731, 1.4889172, 'C$^{18}$O (1$-$0)', color='black', size = 13)
	gc.hide_yaxis_label()

	gc = aplpy.FITSFigure(moments_dir_fits_c18o[1], figure=fig, subplot = (2,4,7))
	gc.show_colorscale(vmin= v_min_1, vmax= v_max_1, interpolation='nearest', cmap=color_map[1])
	gc.ticks.set_xspacing(0.02)
	gc.ticks.set_yspacing(0.025)
	gc.tick_labels.set_xformat('hh:mm:ss'); gc.tick_labels.set_yformat('dd:mm:ss')
	gc.add_colorbar()
	gc.colorbar.set_width(0.2)
	#gc.colorbar.set_location('top')
	gc.colorbar.set_axis_label_text(r'Intensity weighted velocity (km\,s$^{-1}$)')
	gc.show_contour(moments_dir_fits_c18o[0], colors='black', levels= c18o_ii_levs, linewidths=lw, alpha = alpha)
	gc.show_markers(bt12_cores[:,0], bt12_cores[:,1], s = (bt12_cores[:,2] / np.nanmax(bt12_cores[:,2])) * 120 , marker='+', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.show_markers(np.array(table_cores['_RAJ2000']), np.array(table_cores['_DEJ2000']), s = (mass_temp / np.nanmax(mass_temp)) * 100 , marker='x', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.ticks.set_color('white')
	gc.colorbar.set_ticks(ticks_vf)
	gc.hide_yaxis_label()

	gc = aplpy.FITSFigure(moments_dir_fits_c18o[2], figure=fig, subplot = (2,4,8))
	gc.show_colorscale(vmin= v_min_2, vmax= v_max_2, interpolation='nearest', cmap=color_map[1])
	gc.ticks.set_xspacing(0.02)
	gc.ticks.set_yspacing(0.025)
	gc.tick_labels.set_xformat('hh:mm:ss'); gc.tick_labels.set_yformat('dd:mm:ss')
	gc.add_colorbar()
	gc.colorbar.set_width(0.2)
	#gc.colorbar.set_location('top')
	gc.colorbar.set_axis_label_text(r'Intensity weighted width (km\,s$^{-1}$)')
	gc.show_contour(moments_dir_fits_c18o[0], colors='black', levels= c18o_ii_levs, linewidths=lw, alpha = alpha)
	gc.show_markers(bt12_cores[:,0], bt12_cores[:,1], s = (bt12_cores[:,2] / np.nanmax(bt12_cores[:,2])) * 120 , marker='+', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.show_markers(np.array(table_cores['_RAJ2000']), np.array(table_cores['_DEJ2000']), s = (mass_temp / np.nanmax(mass_temp)) * 100 , marker='x', edgecolor='black', zorder = 32, facecolor='black', linewidth = 1.25, path_effects=path_effects2)
	gc.ticks.set_color('white')
	gc.colorbar.set_ticks(ticks_lw)
	gc.hide_yaxis_label()

	ax = fig.axes
	for ax1 in ax[::2]:
		ax1.tick_params(color = 'black', direction = 'in', which = 'both', right = True, top = True)
	for ax1 in ax[1::2]:
		ax1.tick_params(color = 'black', direction = 'in', which = 'both',  bottom = True, top = True, labeltop = True, labelbottom = False)

	## Tick rotation
	lab = [0, 2, 4, 6, 8, 10, 12, 14]
	for x in lab:
	    ax = fig.axes[x]
	    for tick in ax.get_yticklabels():
	        tick.set_rotation(90)
	##


	plt.tight_layout()
	fig.set_rasterized(False)
	plt.savefig('moment_maps_cloudf.pdf', dpi = 200)
	plt.close()
