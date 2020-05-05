import numpy as np
import astropy.units as au
import astropy.stats as stats
from astropy.io import fits
from spectral_cube import SpectralCube
from scipy import ndimage
from scipy.ndimage.morphology import binary_dilation, binary_closing

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

# create circular mask
def get_circmask(h=11, w=11, l=11, center='', radius=5):

    """Creates a sperical mask within 3D cube
        For use in creating structure for e.g. ndimage.dilation

    Parameters
    ----------
    h : int
        Array height
    w : int
        Array width
    l : int
        Array length
    c : int
        Centre of array; where to put circle
    r : int
        Radius of circle mask
        Must be smaller than the size of array
        e.g. set to the beamsize of observations

    Returns
    -------
    mask : array
        Output mask
    """

    if center=='':
        center = int(h/2)

    X, Y, Z = np.ogrid[:h, :w, :l]
    dist = np.sqrt((X - center)**2 + (Y-center)**2 + (Z-center)**2)
    mask = dist <= radius

    return mask *1

def get_prunemask(mask, thresh):

    """Creates a pruned mask such that no regions contain less pixels than
        thresh value; used to remove sub-beam structures within mask

    Parameters
    ----------
    mask : array
        Input mask to be pruned
    thesh : int
        Number of pixel to prune within the mask
        e.g. set the beam area in pixels

    Returns
    -------
    mask : array
        Output prunded mask
    """

    for k in range(mask.shape[0]):

        mask_ = mask[k, :, :]
        l, j = ndimage.label(mask_)
        hist = ndimage.measurements.histogram(l, 0, j+1, j+1)
        os = ndimage.find_objects(l)

        for i in range(j):
            if hist[i+1]<thresh:
                mask_[os[i]] = 0

        mask[k, :, :] = mask_

    return(mask)


def get_expmask(cube, rms, beamarea=2):

    """Get full expanded mask
       Wrapper for getting threshold mask
       binary dilation and closing
       and then pruning with above some beam area
       [under testing]"""

    bmaj = cube.header['BMAJ']
    bmin = cube.header['BMIN']
    pix = np.absolute(cube.header['CDELT1'])
    bmajp = bmaj/pix
    bminp = bmin/pix
    bareap = np.pi*bmajp*bminp
    rad = np.floor(bmajp/2)

    structure = get_circmask(radius=rad)

    mask_l_, _ = get_threshmask(cube, rms.data, thresh=2)
    mask_h_, _ = get_threshmask(cube, rms.data, thresh=5)

    mask_l = mask_l_.include()
    mask_h = mask_h_.include()

    mask_d = binary_dilation(mask_h, iterations=-1, mask=mask_l)
    mask_dc = binary_closing(mask_d, structure=structure, iterations=1)
    mask_dcc = get_prunemask(mask_dc, thresh=bareap*2)

    return(mask_dcc)
