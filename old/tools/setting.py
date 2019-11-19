import numpy as np
import astropy.units as u
from astropy.io import fits
from spectral_cube import SpectralCube

import smooth


def get_cube(inputfile, unit = u.K):
    
    img, hdr = fits.getdata(inputfile, header = True)
    cube = SpectralCube.read(inputfile, format = 'fits', unit = unit, beam_threshold = 0.05)
    cube = cube.with_spectral_unit(u.km / u.s)
    
    cube.beam_threshold = 0.5
    
    return img, hdr, cube


def get_rms_velocities(cube, velocity_range):
    
    spectral_extrema = cube.spectral_extrema
    rms_velocities = [[spectral_extrema[0].value, velocity_range[0]],
                        [velocity_range[1], spectral_extrema[1].value]]

    return rms_velocities


def getrms(cube, rms_velocities, outputfile = ''):

    count = 0
    rms_arr = np.empty([len(rms_velocities), cube.shape[1], cube.shape[2]]) * np.nan
    for rms_velocity in rms_velocities:

        rms_velocity = rms_velocity * (u.km / u.s)
        rms_cube = cube.spectral_slab(rms_velocity[0], rms_velocity[1])
        rms_arr[count, :, :] = np.nanstd(rms_cube, axis = 0)
        
        count += 1

    rms_map = np.nanmean(rms_arr, axis = 0)
    
    if outputfile != '':
        
        header = cube.header
        del header['*3']
        del header['PV1_*']
        header['WCSAXES'] = 2
        
        rms_map_hdu = fits.PrimaryHDU(data = rms_map, header = header)
        rms_map_hdu.writeto(outputfile, overwrite = True)
    
    return rms_map
    
    
def set_mask(cube, rms_map, velocity = '', spectal_units = (u.km / u.s), 
            val_thresh = 0.0, sigma_thresh = 0,
            ):
                
    "Masks array with the set rms file or value and returns the masked data and the velocity axis"
            
    if velocity != '':
        cube = cube.spectral_slab(velocity[0]*spectal_units, velocity[1]*spectal_units)
    else: 
        cube = cube
    
    if val_thresh != 0:  
        print '\nMasking with a value of %0.2fK \n' %(val_thresh)
        mask = cube > val_thresh
        
    if sigma_thresh != 0:
        print '\nMasking with an array with an average %i rms value of %0.2fK \n' %(sigma_thresh, (np.nanmean(rms_map) * sigma_thresh))
        mask = cube > ((rms_map * u.K) * sigma_thresh)
        
    cube = cube.with_mask(mask)   
                        
    return mask, cube
    
    
def set_mask_to_other(cube, mask):
    
    cube = cube.with_mask(mask)    
                    
    return mask, cube    
    
    
def getstatistics(map):

    rms_statistics = np.array([np.nanmin(map.data),
                                np.nanmean(map.data),
                                np.nanmedian(map.data),
                                np.nanmax(map.data),
                                np.nanstd(map.data),
                                np.nanvar(map.data)])
                                
    return rms_statistics
    