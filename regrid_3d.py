import numpy as np
import astropy.units as au

def get_spectralregrid(cube, res_t, res_c=''):

    """Spectrally smooth data
    
    Parameters
    ----------
    cube : spectral cube object
    res_t : float
        target *channel width* (in kms-1)
        *not* resolusion - do get_spectralsmooth
    res_c : float (optional) 
        current resolution (in kms-1) if not in header
    
    Returns
    -------
    hdu_out = Output hdu containing the smoothed data and updated header
    """

    if res_c == '':
        # Assuming beams are in kms-1
        res_c = cube.header['CDELT3']
    pix_c = cube.header['CDELT3']

    spectral_axis = cube.spectral_axis
    spec_min = spectral_axis.min()
    spec_max = spectral_axis.max()
    
    new_spectral_axis = np.arange(spec_min.value, spec_max.value, res_t) *au.km/au.s
    interpcube = cube.spectral_interpolate(new_spectral_axis, suppress_smooth_warning=True)

    return interpcube