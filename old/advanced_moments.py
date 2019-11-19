import numpy as np  
import astropy.units as u
from astropy.io import fits

from tools import *
  
def make_momentmaps(inputfile, outputfile, velocity_range = '',
                    spectralsmoothfactor = 2, spatialsmoothfactor = 2,
                    moments = ['mom0', 'max', 'mom1', 'sigma', 'mom9', 'mom9-mom1'], 
                    singlefile = True, sigma_thresh=5, unit=u.K,
                    ):
    
    __version__ = '0.1'
    __author__  = 'Ashley Barnes'

    """
    Module to produce moment maps from input .fits data cube. 

    Parameters
    ----------
        inputfile : string
            Input file name.
        outputfile: string
            Output file name (without file extension)
        velocity_range: list
            List of velocity range for emission () - i.e. that will not be taken when determining the rms.
        spectralsmoothfactor: float  
            Factor the spectra axis will be smoothed to produce lower resolution cube for mask.
        spatialsmoothfactor: float
            Factor the spatial axis will be smoothed to produce lower resolution cube for mask.
        moments: list
            List of required output moments.
        singlefile: bool
            If all the momoents should be put into a single .fits HDU, otherwise the moments will be output to multiple .fits files. 
        sigma_thresh: float
            Factor of rms used to mask the lower resolution cube mask [looking like 5sigma is the best, but needs testing].
        unit: astropy.unit
            Unit for input datacube.
    Returns
    -------
        momentmaps: .fits
            Will create one or multiple (depending on singlefile) .fits files that contain the moment maps.

    """

    img, hdr, cube = setting.get_cube(inputfile, unit=unit)
    
    cube_sm = smooth.spectral(cube, spectralsmoothfactor = spectralsmoothfactor, just_smooth = True)
    cube_sm = smooth.spatial(cube_sm, spatialsmoothfactor = spatialsmoothfactor, just_smooth = True)

    rms_velocities = setting.get_rms_velocities(cube, velocity_range)                        
    
    cube_rms = setting.getrms(cube, rms_velocities)
    cube_rms_sm = setting.getrms(cube_sm, rms_velocities)
    
    cube_sm_mask, cube_sm_masked = setting.set_mask(cube_sm, cube_rms_sm, sigma_thresh = 5)
    cube_sm_mask_dilation = smooth.smooth_mask(cube_sm_mask, cube.wcs, diam = 3)
    
    _, cube_masked = setting.set_mask_to_other(cube, cube_sm_mask_dilation)
    _, cube_masked_nozero = setting.set_mask(cube_masked, cube_rms, val_thresh = 0)
    
    if outputfile.endswith('.fits'):
        outputfile = outputfile[:-5]
        print 'Removing output .fits'
    
    momentmaps.get_momentmaps(cube_masked_nozero, outputfile, 
                            moments = moments, singlefile = singlefile)
    
    return 