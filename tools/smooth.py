import math
import datetime
import numpy as np
from radio_beam import Beam

import astropy.units as u
from astropy.io import fits
from astropy.convolution import Gaussian1DKernel

import spectral_cube
from spectral_cube import SpectralCube, BooleanArrayMask
from scipy import ndimage


def spectral(input, target_resolution = '', spectralsmoothfactor = '', outputfile = '', velocity_range = '', 
                    current_resolution = '', target_pixel_scale = '', just_smooth = False):
    """Function for producing spectrally smoothed fits files."""

    if type(input) == spectral_cube.spectral_cube.SpectralCube: 
        
        img, hdr = np.asarray(input), input.header
        cube = input
        
    else: 

        img, hdr = fits.getdata(input, header = True)
        cube = SpectralCube.read(input, format = 'fits', unit = u.K, beam_threshold = 0.05)
        cube = cube.with_spectral_unit(u.km / u.s)

    if spectralsmoothfactor == 0: 
        return spectralsmoothfactor

    if current_resolution == '':
        # Assuming beams are in degrees
        current_resolution = hdr['CDELT3']

    fwhm_factor = np.sqrt(8*np.log(2))
    current_resolution = current_resolution * u.km/u.s
    
    if spectralsmoothfactor != '':
        target_resolution = current_resolution * spectralsmoothfactor
    else: 
        target_resolution = target_resolution * u.km/u.s
        
    current_pixel_scale = current_resolution

    gaussian_width = ((target_resolution**2 - current_resolution**2)**0.5 /
                    current_pixel_scale / fwhm_factor)

    print "Current spectral axis has resolution of %0.2f kms" % (current_pixel_scale.value)
    print "New spectral axis has resolution of %0.2f kms" % (target_resolution.value)
    print "Smoothing spectral axis with a Gaussian kernal of %0.2f kms" % ((gaussian_width*current_pixel_scale).value)
    
    kernel = Gaussian1DKernel(gaussian_width)
    smcube = cube.spectral_smooth(kernel)

    if just_smooth:
        
        print "\nReturning smoothed only cube (i.e. no regridding)."
        
        return smcube

    if velocity_range == '':
        spectral_axis = cube.spectral_axis
        velocity_range = [spectral_axis.min(), spectral_axis.max()]
    else:
        velocity_range = velocity_range * u.km/u.s

    if target_pixel_scale == '':
        target_pixel_scale = target_resolution
    else:
        target_pixel_scale = target_pixel_scale * u.km/u.s

    new_spectral_axis = np.arange(velocity_range[0].value, velocity_range[1].value,
                                    target_pixel_scale.value ) * u.km/u.s
    
    print "Regridding spectral axis to %0.2f kms pixels" % ((target_resolution).value)
    interpcube = smcube.spectral_interpolate(new_spectral_axis, suppress_smooth_warning = True)
    
    if outputfile != '':

        interpcube.write(outputfile, overwrite = True)
    
    return interpcube


def spatial(input, target_resolution = '', target_pixel_scale = '', spatialsmoothfactor = '', 
                            outputfile = '', current_resolution = '', target_map_size = '', 
                            just_smooth = True,
                            ):
                            
    """Function for producing spectrally smoothed fits files."""

    img, hdr = np.asarray(input), input.header
    cube = input
    
    if spatialsmoothfactor == 0: 
        return spectralsmoothfactor
            
    if current_resolution == '':
    
        current_resolution_beam = Beam.from_fits_header(hdr)
        print "Current spatial axis has resolution of:"     
        print current_resolution
    
    else:
    
        current_resolution = Beam(major = current_resolution[0]*u.arcsec, 
                                    minor = current_resolution[1]*u.arcsec, 
                                    pa = current_resolution[2] *u.deg)
                                    
        print "Current spatial axis has resolution of:"     
        print current_resolution
    
    if spatialsmoothfactor != '':
        target_resolution = Beam(major = current_resolution_beam.major * spatialsmoothfactor, 
                                    minor = current_resolution_beam.minor * spatialsmoothfactor, 
                                    pa = current_resolution_beam.pa)
                                
        print "New spatial axis has resolution of:" 
        print target_resolution
    
    else :
        target_resolution = Beam(major = target_resolution*u.arcsec, minor = target_resolution*u.arcsec, pa = 0*u.deg)
        
        print "New spatial axis has resolution of:" 
        print target_resolution
        
    try: 
        smcube = cube.convolve_to(target_resolution)
        
    except ValueError:
        smcube = cube
        
    if just_smooth:
        
        print "\nReturning smoothed only cube (i.e. no regridding)."
        
        if outputfile != '':    
        
            smcube.write(outputfile, overwrite = True)
        
        return smcube
        
    
    target_pixel_scale = target_pixel_scale *u.arcsec
    current_pixel_scale = hdr['CDELT2'] * 3600 *u.arcsec
    
    if target_map_size == '':
        print 'Regridding to map the same size as the input map.'
    
        pixel_scale_factor = (target_pixel_scale / current_pixel_scale).value
    
        hdr['CDELT1'] = hdr['CDELT1'] * pixel_scale_factor
        hdr['CDELT2'] = hdr['CDELT2'] * pixel_scale_factor
        hdr['NAXIS1'] = int(hdr['NAXIS1'] / pixel_scale_factor)
        hdr['NAXIS2'] = int(hdr['NAXIS2'] / pixel_scale_factor)
        hdr['CRPIX1'] = int(hdr['CRPIX1'] / pixel_scale_factor)
        hdr['CRPIX2'] = int(hdr['CRPIX2'] / pixel_scale_factor)
    
    else:
        target_map_size = (target_map_size * 3600.) *u.arcsec
        print 'Regridding to map a map size of: %i x %i arcsec' % (target_map_size[0].value, target_map_size[1].value)
    
        pixel_scale_factor = (target_pixel_scale / current_pixel_scale).value
    
        current_map_size = np.array([ np.absolute(hdr['NAXIS1'] * hdr['CDELT1'] * 3600.) , np.absolute(hdr['NAXIS2'] * hdr['CDELT2'] * 3600.) ])
        map_scale_factor = np.array([ (target_map_size[0].value / current_map_size[0]), (target_map_size[1].value / current_map_size[1]) ])
    
    
        hdr['CDELT1'] = hdr['CDELT1'] * pixel_scale_factor
        hdr['CDELT2'] = hdr['CDELT2'] * pixel_scale_factor
    
        hdr['NAXIS1'] = int(np.absolute( target_map_size[0] / (hdr['CDELT1'] *3600 *u.arcsec) ) )
        hdr['NAXIS2'] = int(np.absolute( target_map_size[1] / (hdr['CDELT2'] *3600 *u.arcsec) ) )
        hdr['CRPIX1'] = math.ceil(hdr['NAXIS1'] / 2.)
        hdr['CRPIX2'] = math.ceil(hdr['NAXIS2'] / 2.)
    
        # hdr['CTYPE1'] = 'RA---CAR'
        # hdr['CTYPE2'] = 'DEC--CAR'
        hdr['CTYPE1'] = 'RA---TAN'
        hdr['CTYPE2'] = 'DEC--TAN'
        hdr['RADESYS'] = 'ICRS'
        hdr['EQUINOX'] = 2000.0
        # del hdr['PV*']
        # del hdr['LONPOLE*']
        # del hdr['LATPOLE*']
    
    interpcube = smcube.reproject(hdr, order = u'bilinear')
    
    data = np.asarray(interpcube)
    header = interpcube.header
    header['EQUINOX'] = 2000.0
    del header['ORIGIN']
    del header['DATE']
    header['HISTORY'] = 'Smoothed and regridded using spectral_cube version '+str(spectral_cube.__version__)+' on '+str(datetime.date.today())
    
    if outputfile != '':    
    
        fits.writeto(outputfile, data, header, overwrite = True)
    
    return interpcube
    
    
def smooth_mask(mask, wcs, diam):
    
    mask_np = np.asarray( mask.include() * 1 )        
    structure = np.ones((diam, diam, diam))
    dist = ((np.indices((diam, diam)) - (diam - 1) / 2.)**2).sum(axis=0)**0.5
    structure[dist > diam / 2.] = 0
    
    mask_dilation = ndimage.binary_dilation(mask_np, structure=structure)
    mask_erosion = ndimage.binary_erosion(mask_np, structure=structure)
    
    mask_dilation_speccube = BooleanArrayMask(mask=mask_erosion, wcs=wcs)
    
    return mask_dilation_speccube