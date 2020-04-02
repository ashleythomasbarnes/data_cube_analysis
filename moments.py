import numpy as np
import astropy.units as au
import astropy.stats as stats
from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.utils.console import ProgressBar

def get_cube(inputfile, x_units=(au.km / au.s), y_units=au.K):

    """Get datacube

    Parameters
    ----------
    inputfile : str
        input filename
    x_units:
        units of x-axis
    y_units:
        units of y-axis

    Returns
    -------
    cube : spectral-cube object
        output cube in xunit and yunit
    """

    cube = SpectralCube.read(inputfile, format = 'fits', unit = y_units, beam_threshold = 0.05)
    cube = cube.with_spectral_unit(x_units, velocity_convention='radio')
    cube.beam_threshold = 0.5
    cube.allow_huge_operations=True

    if cube.unit == "Jy / beam":
        cube = cube.to(au.K)

    return cube

def get_rms(cube, rms_velocities, outputfile=''):

    """Get RMS map from a given cube

    Parameters
    ----------
    cube : spectral cube object
        Input cube
    rms_velocities : list of velocities used to return the rms - unitless
        List of rms values - e.g [[20,30],[50,60]]
    outputfile : string (optional)
        Output filename

    Returns
    -------
    rms_map : fits.PrimaryHDU object
        map of rms values
    """


    count = 0
    rms_arr = np.empty([len(rms_velocities), cube.shape[1], cube.shape[2]]) * np.nan

    for rms_velocity in rms_velocities:

        rms_velocity = rms_velocity * (au.km / au.s)
        rms_cube = cube.spectral_slab(rms_velocity[0], rms_velocity[1])

        rms_arr[count, :, :] = np.nanstd(rms_cube, axis=0)
#         rms_arr[count, :, :] = np.sqrt(np.nansum(rms_cube**2, axis=0) / len(rms_cube)) #RMS CALC
        #MAD -> used incase some emission is incl. in rms
        # rms_arr[count, :, :] = stats.mad_std(rms_cube, axis=0)

        count += 1

    rms_map = np.nanmean(rms_arr, axis = 0)

    header = cube.header
    del header['*3']
    del header['PV1_*']
    header['WCSAXES'] = 2
    rms_map = fits.PrimaryHDU(data=rms_map, header=header)

    if outputfile != '':
        rms_map.writeto(outputfile, overwrite = True)

    return rms_map

def get_threshmask(cube, thresh_map, thresh=0):

    """Set a given threshold map as mask to cube (e.g. for rms)

    Parameters
    ----------
    cube : spectral cube object (in Kelvin)
        Input cube
    thresh_map :
        Input map to provide threshold (e.g. RMS map)
    thresh : float (optional)
        Input threshold to create mask

    Returns
    -------
    mask : fits.PrimaryHDU object
        Output mask
    cube : fits.PrimaryHDU object
        Output masked cube
    """

    mask = cube > ((thresh_map * au.K) * thresh)
    cube_masked = cube.with_mask(mask)

    return mask, cube_masked

def get_momentmaps(cube, mom_velocity='', outputfile='',
            moments=['mom0', 'max', 'mom1', 'sigma', 'mom9', 'mom9-mom1', 'variance', 'fwhm'],
            ):

    """Get moment maps from a given (masked) cube

    Parameters
    ----------
    cube : spectral cube object
        Input cube
    outputfile : str (optional)
        Prefix of output file name - [outputfile]_[mom].fits
    mom_velocity : list (optional)
        list of velocities over which to deteremine moment maps - e.g. [40,50]
    moments : list (optional)
        list of moments to output *to file* - can be one or more of the following
        ['mom0', 'max', 'mom1', 'sigma', 'mom9', 'mom9-mom1', 'variance', 'fwhm']
            'mom0': moment 0 or integrated intensity
            'max': peak intensity moment map
            'mom1': (intensity weighted) centroid velocity moment map
            'sigma': (intensity weighted) line-width moment map
            'mom9': velocity at peak intensity moment map
            'mom9-mom1': "velocity complexity" - measure of asymmetry of velocity profile

    Returns
    -------
    mom_out : dict
        dictionary of output moment maps - will be one or more of the following
        ['mom0', 'max', 'mom1', 'sigma', 'mom9', 'mom9-mom1', 'variance', 'fwhm']

        to reduce run time, if 'mom9' or 'mom9-mom1' not chosen as output to file
        they will not be calculated and None value will be returned in dict.
    """

    if mom_velocity!='':
        velo_min = mom_velocity[0]*au.km / au.s
        velo_max = mom_velocity[1]*au.km / au.s
        cube_slab = cube.spectral_slab(velo_min, velo_max)

    else:
        print('[INFO] No velocity range.')
        cube_slab = cube

    shape = cube_slab.shape

    if 'mom0' in moments:
        moment_0 = cube_slab.moment(order=0, axis=0)
    else:
        moment_0 = None

    if 'max' in moments:
        maximum = cube_slab.max(axis = 0)
    else:
        maximum = None

    if 'mom1' in moments:
        moment_1 = cube_slab.moment(order=1, axis=0)
    else:
        moment_1 = None

    if 'sigma' in moments:
        sigma_map = cube_slab.linewidth_sigma()
    else:
        sigma_map = None

    if 'variance' in moments:
        variance_map = cube_slab.moment(order=2, axis=0)
    else:
        variance_map = None

    if 'fwhm' in moments:
        fwhm_map = cube_slab.linewidth_fwhm()
    else:
        fwhm_map = None

    # moment_0 = cube_slab.moment(order=0, axis=0)
    # maximum = cube_slab.max(axis = 0)
    # moment_1 = cube_slab.moment(order=1, axis=0)
    # sigma_map = cube_slab.linewidth_sigma()
    # variance_map = cube_slab.moment(order=2, axis=0)
    # fwhm_map = cube_slab.linewidth_fwhm()
    # moment_9 = None
    # velo_complexity = None

    if 'mom9' in moments:

        print('[INFO] Creating mom9 and velo_complexity - make take a while.')

        moment_9 = np.empty([shape[1], shape[2]]) * np.nan
        moment_9_ID = cube_slab.argmax(axis=0)
        spectral_axis = cube_slab.spectral_axis

        for velo in range(1, len(spectral_axis)):
            for x in range(shape[1]):
                for y in range(shape[2]):

                    if moment_9_ID[x,y] == 0:
                        moment_9[x,y] = np.nan
                    else:
                        moment_9[x,y] = spectral_axis[moment_9_ID[x,y]].value

        velo_complexity = np.absolute(np.asarray(moment_1) - np.asarray(moment_9))

        moment_9 = moment_9 *au.km/au.s
        velo_complexity = velo_complexity *au.km/au.s

        moment_9_hdu = fits.PrimaryHDU(moment_9.value, moment_1.header)
        velo_complexity_hdu = fits.PrimaryHDU(velo_complexity.value, moment_1.header)

    else:
        moment_9_hdu = None
        velo_complexity_hdu = None

    mom_out = dict.fromkeys(moments)
    mom_out['mom0'] = moment_0
    mom_out['max'] = maximum
    mom_out['mom1'] = moment_1
    mom_out['sigma'] = sigma_map
    mom_out['mom9'] = moment_9_hdu
    mom_out['mom9-mom1'] = velo_complexity_hdu
    mom_out['variance'] = variance_map
    mom_out['fwhm'] = fwhm_map

    if outputfile!='':
        if 'mom0' in moments:
            moment_0.write(outputfile+'_mom0.fits', overwrite=True)
        if 'max' in moments:
            maximum.write(outputfile+'_max.fits', overwrite=True)
        if 'mom1' in moments:
            moment_1.write(outputfile+'_mom1.fits', overwrite=True)
        if 'sigma' in moments:
            sigma_map.write(outputfile+'_sigma.fits', overwrite=True)
        if 'mom9' in moments:
            moment_9_hdu.writeto(outputfile+'_mom9.fits', overwrite=True)
        if 'mom9-mom1' in moments:
            velo_complexity.writeto(outputfile+'_mom9-mom1.fits', overwrite=True)
        if 'variance' in moments:
            variance_map.write(outputfile+'_variance.fits', overwrite=True)
        if 'fwhm' in moments:
            fwhm_map.write(outputfile+'_fwhm.fits', overwrite=True)

    return mom_out

def get_nchan(mask):

    """Get map of number of channels along each pixel within a mask

    Parameters
    ----------
    mask : (2d) np.array
        Input mask

    Returns
    -------
    nchans : (2d) np.array
        Output map of number of channels
    """

    arr1 = mask.include()
    arr2 = arr1 * 1
    nchan = np.sum(arr2, axis=0)

    return nchan

def get_mom0err(nchan, rms_map, velo_res):

    """Get map of err on mom0
    sigma_W = rms * V_res * N_chan**0.5

    Parameters
    ----------
    nchan : (2d) np.array
        Input map of number of channels
    rms_map : fits hdu object
        Input map of rms (inc header for output hdu)
    velo_res : float
        Input velocity resolution (in kms-1 usally)

    Returns
    -------
    mom0err : fits hdu object
        Output map of err on mom0
    """


    cdelt = np.absolute(velo_res)

    nchan_sqrt = np.sqrt(nchan)
    mom0err = rms_map.data * nchan_sqrt * velo_res

    id_nan = np.where(mom0err==0)
    mom0err[id_nan] = np.nan

    header = rms_map.header
    header['BUNIT'] = 'K km s-1'

    mom0err_hdu = fits.PrimaryHDU(mom0err, header)

    return mom0err_hdu
